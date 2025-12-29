"""
H5AD Service for handling single-cell h5ad data
Provides lazy loading, caching, and efficient data extraction
"""
import os
import json
import logging
from typing import Optional, List, Dict, Any
from pathlib import Path
import base64
from io import BytesIO
import sqlite3
from contextlib import contextmanager

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import linkage, leaves_list
from PIL import Image

logger = logging.getLogger(__name__)


class H5ADService:
    """Service for managing h5ad file operations"""
    
    def __init__(self, h5ad_base_path: str = "h5ad", precomputed_base_path: str = "h5ad/precomputed"):
        self.h5ad_base_path = Path(h5ad_base_path)
        self.precomputed_base_path = Path(precomputed_base_path)
        self.cache = {}  # In-memory cache for frequently accessed data
        
    @contextmanager
    def get_db_connection(self, dataset_id: str):
        """Context manager for SQLite connections"""
        db_path = self.get_precomputed_dir(dataset_id) / "data.db"
        
        if not db_path.exists():
            yield None
            return
            
        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row  # Enable column access by name
        try:
            yield conn
        finally:
            conn.close()
        
    def get_h5ad_path(self, dataset_id: str) -> Path:
        """Get the path to the h5ad file"""
        return self.h5ad_base_path / "raw" / f"{dataset_id}.h5ad"
    
    def get_precomputed_dir(self, dataset_id: str) -> Path:
        """Get the directory for precomputed data"""
        return self.precomputed_base_path / dataset_id
    
    def _load_json(self, file_path: Path) -> Dict:
        """Load JSON file"""
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"Error loading JSON from {file_path}: {e}")
            raise
    
    def _save_json(self, file_path: Path, data: Dict) -> None:
        """Save data as JSON"""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with open(file_path, 'w') as f:
                json.dump(data, f)
        except Exception as e:
            logger.error(f"Error saving JSON to {file_path}: {e}")
            raise
    
    def _load_h5ad(self, dataset_id: str):
        """Load h5ad file (lazy loading)"""
        try:
            import anndata
            h5ad_path = self.get_h5ad_path(dataset_id)
            if not h5ad_path.exists():
                raise FileNotFoundError(f"H5AD file not found: {h5ad_path}")
            return anndata.read_h5ad(h5ad_path)
        except Exception as e:
            logger.error(f"Error loading h5ad file for {dataset_id}: {e}")
            raise
    
    def get_embedding(self, dataset_id: str, embedding_type: str = "umap") -> Dict:
        """
        Get embedding coordinates (UMAP, tSNE, PCA)
        Returns: {coordinates: [[x, y], ...], cell_ids: [...]}
        """
        cache_key = f"{dataset_id}:{embedding_type}"
        
        # Check in-memory cache
        if cache_key in self.cache:
            logger.info(f"Cache hit for {cache_key}")
            return self.cache[cache_key]
        
        # Check precomputed file
        precomputed_file = self.get_precomputed_dir(dataset_id) / f"{embedding_type}.json"
        if precomputed_file.exists():
            logger.info(f"Loading precomputed {embedding_type} for {dataset_id}")
            data = self._load_json(precomputed_file)
            self.cache[cache_key] = data
            return data
        
        # Fallback: extract from h5ad
        logger.info(f"Extracting {embedding_type} from h5ad for {dataset_id}")
        adata = self._load_h5ad(dataset_id)
        
        # Map embedding type to obsm key
        embedding_key_map = {
            "umap": "X_umap",
            "tsne": "X_tsne",
            "pca": "X_pca"
        }
        
        obsm_key = embedding_key_map.get(embedding_type.lower())
        if not obsm_key or obsm_key not in adata.obsm:
            raise ValueError(f"Embedding {embedding_type} not found in dataset {dataset_id}")
        
        coordinates = adata.obsm[obsm_key].tolist()
        cell_ids = adata.obs_names.tolist()
        
        data = {
            "embedding_type": embedding_type,
            "coordinates": coordinates,
            "cell_ids": cell_ids
        }
        
        # Cache for future use
        self._save_json(precomputed_file, data)
        self.cache[cache_key] = data
        
        return data
    
    def get_gene_expression(self, dataset_id: str, gene: str, use_raw: bool = False) -> Dict:
        """
        Get gene expression for a specific gene
        
        Args:
            dataset_id: Dataset identifier
            gene: Gene symbol
            use_raw: If True, use raw counts from adata.raw.X (if available), otherwise use normalized
        
        Returns: {gene: str, expression: [...], cell_ids: [...], stats: {...}}
        """
        # Check precomputed cache
        cache_suffix = "_raw" if use_raw else ""
        gene_cache_dir = self.get_precomputed_dir(dataset_id) / "gene_expression"
        gene_cache_file = gene_cache_dir / f"{gene}{cache_suffix}.json"
        
        if gene_cache_file.exists():
            logger.info(f"Loading cached gene expression for {gene} in {dataset_id}")
            return self._load_json(gene_cache_file)
        
        # Load from h5ad on-demand
        logger.info(f"Extracting gene expression for {gene} from h5ad")
        adata = self._load_h5ad(dataset_id)
        
        # Check if gene exists
        if gene not in adata.var_names:
            raise ValueError(f"Gene {gene} not found in dataset {dataset_id}")
        
        # Extract expression (use normalized counts from .X, or raw from .raw.X if requested)
        if use_raw and hasattr(adata, 'raw') and adata.raw is not None:
            logger.info(f"Using raw counts for {gene}")
            gene_expr = adata.raw[:, gene].X
        else:
            # Use normalized counts from .X (default, recommended for visualization)
            gene_expr = adata[:, gene].X
        
        # Handle sparse matrix
        if hasattr(gene_expr, 'toarray'):
            expression = gene_expr.toarray().flatten()
        else:
            expression = np.array(gene_expr).flatten()
        
        expression_list = expression.tolist()
        cell_ids = adata.obs_names.tolist()
        
        # Compute statistics
        non_zero = expression[expression > 0]
        stats = {
            "mean": float(np.mean(expression)),
            "max": float(np.max(expression)),
            "min": float(np.min(expression)),
            "expressed_cells": int(np.sum(expression > 0)),
            "mean_expressed": float(np.mean(non_zero)) if len(non_zero) > 0 else 0.0,
            "median": float(np.median(expression)),
            "std": float(np.std(expression))
        }
        
        result = {
            "gene": gene,
            "expression": expression_list,
            "cell_ids": cell_ids,
            "stats": stats,
            "data_type": "raw" if use_raw else "normalized"
        }
        
        # Cache for future use
        self._save_json(gene_cache_file, result)
        
        return result
    
    def get_metadata(self, dataset_id: str, columns: Optional[List[str]] = None) -> Dict:
        """
        Get cell metadata
        Returns: {columns: [...], data: {...}, cell_ids: [...]}
        """
        # Check precomputed file
        metadata_file = self.get_precomputed_dir(dataset_id) / "metadata.json"
        
        if metadata_file.exists():
            logger.info(f"Loading precomputed metadata for {dataset_id}")
            data = self._load_json(metadata_file)
            
            # Filter columns if requested
            if columns:
                filtered_data = {col: data["data"][col] for col in columns if col in data["data"]}
                return {
                    "columns": columns,
                    "data": filtered_data,
                    "cell_ids": data["cell_ids"]
                }
            return data
        
        # Extract from h5ad
        logger.info(f"Extracting metadata from h5ad for {dataset_id}")
        adata = self._load_h5ad(dataset_id)
        
        metadata_dict = adata.obs.to_dict(orient='list')
        cell_ids = adata.obs_names.tolist()
        
        # Convert any non-serializable types
        for key, values in metadata_dict.items():
            metadata_dict[key] = [str(v) if pd.isna(v) else v for v in values]
        
        result = {
            "columns": list(adata.obs.columns),
            "data": metadata_dict,
            "cell_ids": cell_ids
        }
        
        # Cache
        self._save_json(metadata_file, result)
        
        if columns:
            filtered_data = {col: result["data"][col] for col in columns if col in result["data"]}
            return {
                "columns": columns,
                "data": filtered_data,
                "cell_ids": result["cell_ids"]
            }
        
        return result
    
    def get_metadata_column_values(self, dataset_id: str, column: str) -> Dict:
        """
        Get unique values and counts for a specific metadata column
        
        Args:
            dataset_id: Dataset identifier
            column: Metadata column name
            
        Returns:
            {
                "column": column_name,
                "values": [unique_values],
                "counts": {value: count}
            }
        """
        # Load h5ad
        adata = self._load_h5ad(dataset_id)
        
        # Check if column exists
        if column not in adata.obs.columns:
            raise ValueError(f"Column '{column}' not found in dataset {dataset_id}. "
                           f"Available columns: {list(adata.obs.columns)}")
        
        # Get value counts
        value_counts = adata.obs[column].value_counts()
        
        # Convert to dict and handle NaN
        counts_dict = {}
        for value, count in value_counts.items():
            key = str(value) if not pd.isna(value) else "NA"
            counts_dict[key] = int(count)
        
        # Get sorted unique values
        unique_values = sorted(counts_dict.keys())
        
        return {
            "column": column,
            "values": unique_values,
            "counts": counts_dict
        }
    
    def get_plot_data(
        self, 
        dataset_id: str, 
        embedding: str = "umap",
        genes: Optional[List[str]] = None,
        metadata_cols: Optional[List[str]] = None
    ) -> Dict:
        """
        Get combined data for plotting (optimized endpoint)
        Returns everything needed for a complete visualization in one call
        """
        result = {}
        
        # Get embedding coordinates
        embedding_data = self.get_embedding(dataset_id, embedding)
        result["coordinates"] = embedding_data["coordinates"]
        result["cell_ids"] = embedding_data["cell_ids"]
        
        # Get gene expression if requested
        if genes:
            result["genes"] = {}
            for gene in genes:
                try:
                    gene_data = self.get_gene_expression(dataset_id, gene)
                    result["genes"][gene] = gene_data["expression"]
                except ValueError as e:
                    logger.warning(f"Gene {gene} not found: {e}")
                    # Return empty list to satisfy response schema and signal missing gene
                    result["genes"][gene] = []
        
        # Get metadata if requested
        if metadata_cols:
            metadata = self.get_metadata(dataset_id, metadata_cols)
            result["metadata"] = metadata["data"]
        
        return result
    
    def search_genes(self, dataset_id: str, search_term: str, limit: int = 50) -> List[Dict]:
        """
        Search for genes by symbol
        Returns: [{symbol: str, mean_expression: float, expressed_cells: int}, ...]
        """
        # Try SQLite first (fast search)
        with self.get_db_connection(dataset_id) as conn:
            if conn:
                try:
                    logger.info(f"Searching genes in SQLite for {dataset_id}")
                    search_lower = search_term.lower()
                    
                    # Optimize query: use prefix search when possible for better index usage
                    # If search term doesn't contain wildcards, try exact match first, then prefix
                    if '%' not in search_term and '_' not in search_term:
                        # Try exact match first (fastest)
                        cursor = conn.execute("""
                            SELECT gene, stats_json 
                            FROM gene_stats 
                            WHERE gene_lower = ? 
                            ORDER BY mean_expression DESC 
                            LIMIT ?
                        """, (search_lower, limit))
                        
                        results = []
                        for row in cursor:
                            stats = json.loads(row['stats_json'])
                            results.append({
                                "symbol": row['gene'],
                                "mean_expression": stats["mean_expression"],
                                "expressed_cells": stats["expressed_cells"]
                            })
                        
                        # If exact match found enough results, return early
                        if len(results) >= limit:
                            return results
                        
                        # Otherwise, fall through to prefix search
                        cursor = conn.execute("""
                            SELECT gene, stats_json 
                            FROM gene_stats 
                            WHERE gene_lower LIKE ? 
                            ORDER BY mean_expression DESC 
                            LIMIT ?
                        """, (f'{search_lower}%', limit))
                    else:
                        # Use pattern search (contains)
                        cursor = conn.execute("""
                            SELECT gene, stats_json 
                            FROM gene_stats 
                            WHERE gene_lower LIKE ? 
                            ORDER BY mean_expression DESC 
                            LIMIT ?
                        """, (f'%{search_lower}%', limit))
                    
                    results = []
                    for row in cursor:
                        stats = json.loads(row['stats_json'])
                        results.append({
                            "symbol": row['gene'],
                            "mean_expression": stats["mean_expression"],
                            "expressed_cells": stats["expressed_cells"]
                        })
                    
                    return results
                except Exception as e:
                    logger.error(f"Error searching genes in SQLite: {e}")
                    raise
        
        # If no DB connection (should effectively not happen if preprocessed correctly)
        return []
    
    def calculate_module_score(
        self,
        dataset_id: str,
        gene_list: List[str],
        use_raw: bool = False
    ) -> Dict:
        """
        Calculate module score for a list of genes using scanpy.tl.score_genes
        (Seurat-style scoring: average expression minus reference set average)
        
        Args:
            dataset_id: Dataset identifier
            gene_list: List of gene symbols
            use_raw: Use raw counts if True, normalized if False
        
        Returns:
            {
                "module_score": [...],
                "cell_ids": [...],
                "genes_used": [...],
                "genes_not_found": [...],
                "method": "scanpy_score_genes",
                "stats": {...}
            }
        """
        logger.info(f"Calculating module score for {len(gene_list)} genes in {dataset_id} using scanpy.tl.score_genes")
        
        # Load h5ad
        adata = self._load_h5ad(dataset_id)
        
        # Find which genes exist
        # Check in var_names (for normalized) or raw.var_names (for raw)
        if use_raw and hasattr(adata, 'raw') and adata.raw is not None:
            var_names = adata.raw.var_names
        else:
            var_names = adata.var_names
        
        genes_found = [g for g in gene_list if g in var_names]
        genes_not_found = [g for g in gene_list if g not in var_names]
        
        if len(genes_found) == 0:
            raise ValueError(f"None of the genes found in dataset {dataset_id}")
        
        logger.info(f"Found {len(genes_found)}/{len(gene_list)} genes")
        
        # Make a copy to avoid modifying the original
        adata_copy = adata.copy()
        
        # Use scanpy.tl.score_genes
        # This calculates: average expression of gene_list - average expression of reference set
        sc.tl.score_genes(
            adata_copy,
            gene_list=genes_found,
            ctrl_size=50,
            n_bins=25,
            score_name='module_score',
            use_raw=use_raw if hasattr(adata_copy, 'raw') and adata_copy.raw is not None else None
        )
        
        # Extract scores
        module_score = adata_copy.obs['module_score'].values
        module_score_list = module_score.tolist()
        cell_ids = adata_copy.obs_names.tolist()
        
        # Compute statistics
        non_zero = module_score[module_score > 0]
        stats = {
            "mean": float(np.mean(module_score)),
            "median": float(np.median(module_score)),
            "std": float(np.std(module_score)),
            "min": float(np.min(module_score)),
            "max": float(np.max(module_score)),
            "cells_with_expression": int(np.sum(module_score > 0)),
            "mean_expressed": float(np.mean(non_zero)) if len(non_zero) > 0 else 0.0
        }
        
        return {
            "module_score": module_score_list,
            "cell_ids": cell_ids,
            "genes_used": genes_found,
            "genes_not_found": genes_not_found,
            "method": "scanpy_score_genes",
            "n_genes_used": len(genes_found),
            "stats": stats,
            "data_type": "raw" if use_raw else "normalized"
        }
    
    def generate_heatmap(
        self,
        dataset_id: str,
        gene_list: List[str],
        plot_type: str = "heatmap",
        scale_expression: bool = False,
        cluster_rows: bool = False,
        cluster_columns: bool = False,
        group_by: Optional[str] = "cell_type"
    ) -> Dict:
        """
        Generate heatmap or dotplot data for genes across cell types
        
        Args:
            dataset_id: Dataset identifier
            gene_list: List of gene symbols
            plot_type: "heatmap" or "dotplot"
        
        Returns:
            {
                "data": 2D array (genes × cell_types),
                "genes": [...],
                "cell_types": [...],
                "genes_used": [...],
                "genes_not_found": [...],
                "plot_type": "heatmap" or "dotplot",
                "stats": {...}
            }
        """
        if plot_type not in ["heatmap", "dotplot"]:
            raise ValueError(f"plot_type must be 'heatmap' or 'dotplot', got '{plot_type}'")
        
        logger.info(f"Generating {plot_type} for {len(gene_list)} genes in {dataset_id}")
        
        # Load h5ad
        adata = self._load_h5ad(dataset_id)
        
        # Determine grouping column
        group_col = group_by or "cell_type"
        if group_col not in adata.obs.columns:
            available = list(adata.obs.columns)
            raise ValueError(
                f"Dataset {dataset_id} does not have '{group_col}' column in metadata. "
                f"Available columns: {available[:20]}{'...' if len(available) > 20 else ''}"
            )
        
        # Get unique groups in order of appearance (preserve original order)
        cell_types = adata.obs[group_col].unique().tolist()
        
        if len(cell_types) == 0:
            raise ValueError(f"No groups found in column '{group_col}' for dataset {dataset_id}")
        
        logger.info(f"Found {len(cell_types)} groups in '{group_col}'")
        
        # Find which genes exist and preserve the order from gene_list
        genes_found = [g for g in gene_list if g in adata.var_names]
        genes_not_found = [g for g in gene_list if g not in adata.var_names]
        
        if len(genes_found) == 0:
            raise ValueError(f"None of the genes found in dataset {dataset_id}")
        
        logger.info(f"Found {len(genes_found)}/{len(gene_list)} genes")
        
        # Prepare data matrix
        data = []
        all_expressions = []
        all_percentages = []
        
        # Use normalized data (not raw)
        # Process genes in the order they appear in gene_list (preserve user's order)
        for gene in genes_found:
            gene_idx = adata.var_names.get_loc(gene)
            gene_data = []
            
            for cell_type in cell_types:
                # Filter cells by group
                cell_mask = adata.obs[group_col] == cell_type
                gene_expr = adata[cell_mask, gene_idx].X
                
                # Handle sparse matrix
                if hasattr(gene_expr, 'toarray'):
                    gene_expr = gene_expr.toarray().flatten()
                else:
                    gene_expr = np.array(gene_expr).flatten()
                
                # Calculate average expression (log1p transformed)
                # Apply log1p: log(1 + expression)
                log_expr = np.log1p(gene_expr)
                avg_expr = float(np.mean(log_expr))
                all_expressions.append(avg_expr)
                
                if plot_type == "dotplot":
                    # Calculate percentage of cells expressing (> 0)
                    pct_expressing = float(np.mean(gene_expr > 0) * 100)
                    all_percentages.append(pct_expressing)
                    gene_data.append({
                        "expression": avg_expr,
                        "percentage": pct_expressing
                    })
                else:
                    gene_data.append(avg_expr)
            
            data.append(gene_data)
        
        # Convert to numpy array for easier manipulation
        if plot_type == "heatmap":
            data_array = np.array(data)
        else:
            # Extract expression values for scaling/clustering
            expr_array = np.array([[d['expression'] for d in row] for row in data])
        
        # Apply scaling if requested (before clustering)
        if scale_expression:
            if plot_type == "heatmap":
                # Z-score normalize each row (gene)
                row_means = data_array.mean(axis=1, keepdims=True)
                row_stds = data_array.std(axis=1, keepdims=True) + 1e-10  # Add small epsilon to avoid division by zero
                data_array = (data_array - row_means) / row_stds
                data = data_array.tolist()
            else:
                # Z-score normalize each row (gene) for expression values
                row_means = expr_array.mean(axis=1, keepdims=True)
                row_stds = expr_array.std(axis=1, keepdims=True) + 1e-10
                expr_array = (expr_array - row_means) / row_stds
                # Update expression values in data
                for i, row in enumerate(data):
                    for j, point in enumerate(row):
                        point['expression'] = float(expr_array[i, j])
        
        # Apply clustering if requested
        if cluster_rows:
            if plot_type == "heatmap":
                # Cluster genes (rows) using hierarchical clustering
                linkage_matrix = linkage(data_array, method='average', metric='euclidean')
                row_order = leaves_list(linkage_matrix)
                data = [data[i] for i in row_order]
                genes_found = [genes_found[i] for i in row_order]
                # Update data_array after reordering
                data_array = np.array(data)
            else:
                # Cluster genes based on expression array
                linkage_matrix = linkage(expr_array, method='average', metric='euclidean')
                row_order = leaves_list(linkage_matrix)
                data = [data[i] for i in row_order]
                genes_found = [genes_found[i] for i in row_order]
                # Re-extract expr_array after reordering
                expr_array = np.array([[d['expression'] for d in row] for row in data])
        
        if cluster_columns:
            if plot_type == "heatmap":
                # Cluster cell types (columns)
                linkage_matrix = linkage(data_array.T, method='average', metric='euclidean')
                col_order = leaves_list(linkage_matrix)
                data = [[row[i] for i in col_order] for row in data]
                cell_types = [cell_types[i] for i in col_order]
            else:
                # Cluster cell types based on expression array
                linkage_matrix = linkage(expr_array.T, method='average', metric='euclidean')
                col_order = leaves_list(linkage_matrix)
                data = [[row[i] for i in col_order] for row in data]
                cell_types = [cell_types[i] for i in col_order]
        
        # Calculate statistics (after all transformations)
        if plot_type == "heatmap":
            all_values = [v for row in data for v in row]
        else:
            all_values = [v['expression'] for row in data for v in row]
        
        stats = {
            "min_expression": float(np.min(all_values)),
            "max_expression": float(np.max(all_values)),
            "mean_expression": float(np.mean(all_values))
        }
        
        if plot_type == "dotplot":
            all_pcts = [v['percentage'] for row in data for v in row]
            stats["min_percentage"] = float(np.min(all_pcts))
            stats["max_percentage"] = float(np.max(all_pcts))
        
        return {
            "data": data,
            "genes": genes_found,
            "cell_types": cell_types,
            "genes_used": genes_found,
            "genes_not_found": genes_not_found,
            "plot_type": plot_type,
            "scale_expression": scale_expression,
            "cluster_rows": cluster_rows,
            "cluster_columns": cluster_columns,
            "stats": stats
        }
    
    def _detect_categorical_metadata(self, adata) -> Dict[str, List[str]]:
        """
        Detect categorical metadata columns and return their unique values
        
        Returns:
            Dict mapping column name to list of unique values (for categorical columns only)
        """
        categorical_metadata = {}
        
        for col in adata.obs.columns:
            # Get unique values first to check count
            unique_vals = adata.obs[col].unique()
            
            # Only include if there are reasonable number of unique values (< 1000)
            # This filters out things like cell IDs or batch columns with too many values
            if len(unique_vals) >= 1000:
                continue
            
            # Check if column is categorical, string, object, boolean, or numeric with few unique values
            is_categorical = pd.api.types.is_categorical_dtype(adata.obs[col])
            is_string = pd.api.types.is_object_dtype(adata.obs[col]) or pd.api.types.is_string_dtype(adata.obs[col])
            is_boolean = pd.api.types.is_bool_dtype(adata.obs[col])
            is_numeric_categorical = (pd.api.types.is_numeric_dtype(adata.obs[col]) and len(unique_vals) < 100)
            
            if is_categorical or is_string or is_boolean or is_numeric_categorical:
                # Convert to list and handle NaN values
                unique_list = [str(v) if not pd.isna(v) else "NA" for v in unique_vals]
                categorical_metadata[col] = sorted(unique_list)
        
        return categorical_metadata
    
    def get_dataset_info(self, dataset_id: str) -> Dict:
        """
        Get basic information about the h5ad dataset
        """
        # Check if we have a precomputed info file
        info_file = self.get_precomputed_dir(dataset_id) / "info.json"
        if info_file.exists():
            info = self._load_json(info_file)
            # Check for spatial features if not already in info
            if "spatial_features" not in info:
                spatial_info = self.detect_spatial_features(dataset_id)
                info["spatial_features"] = spatial_info
                self._save_json(info_file, info)
            # Check for categorical metadata if not already in info
            if "categorical_metadata" not in info:
                adata = self._load_h5ad(dataset_id)
                info["categorical_metadata"] = self._detect_categorical_metadata(adata)
                self._save_json(info_file, info)
            return info
        
        # Otherwise, load from h5ad
        adata = self._load_h5ad(dataset_id)
        
        # Detect available embeddings (check both h5ad and cached files)
        embeddings = []
        precomputed_dir = self.get_precomputed_dir(dataset_id)
        
        # Check cached embeddings first
        if (precomputed_dir / "umap.json").exists() or 'X_umap' in adata.obsm:
            embeddings.append('umap')
        if (precomputed_dir / "tsne.json").exists() or 'X_tsne' in adata.obsm:
            embeddings.append('tsne')
        if (precomputed_dir / "pca.json").exists() or 'X_pca' in adata.obsm:
            embeddings.append('pca')
        
        # Detect spatial features
        spatial_info = self.detect_spatial_features(dataset_id)
        
        # Detect categorical metadata
        categorical_metadata = self._detect_categorical_metadata(adata)
        
        info = {
            "dataset_id": dataset_id,
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "embeddings": embeddings,
            "metadata_columns": list(adata.obs.columns),
            "categorical_metadata": categorical_metadata,
            "gene_names": adata.var_names.tolist()[:100],  # Return first 100 genes as sample
            "spatial_features": spatial_info
        }
        
        # Cache the info
        self._save_json(info_file, info)
        
        return info
    
    def detect_spatial_features(self, dataset_id: str) -> Dict:
        """
        Auto-detect available spatial features in the dataset
        Returns: {
            "has_spatial_coordinates": bool,
            "has_svg": bool,
            "has_precomputed_deg": bool,
            "has_deconvolution": bool,
            "has_ccc": bool,
            "dataset_type": "visium" | "xenium" | "unknown" | None
        }
        """
        cache_key = f"{dataset_id}:spatial_features"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Check precomputed file
        spatial_info_file = self.get_precomputed_dir(dataset_id) / "spatial_info.json"
        if spatial_info_file.exists():
            info = self._load_json(spatial_info_file)
            self.cache[cache_key] = info
            return info
        
        # Load h5ad to detect features
        try:
            adata = self._load_h5ad(dataset_id)
            
            # Check for spatial coordinates
            has_spatial = any(key in adata.obsm for key in ['spatial', 'spatial_xy', 'X_spatial'])
            
            # Check for SVG data
            svg_cols = ['gft_score', 'svg_rank', 'cutoff_gft_score', 'pvalue', 'fdr']
            has_svg = any(col in adata.var.columns for col in svg_cols)
            
            # Check for precomputed DEG
            has_deg = "rank_genes_groups" in adata.uns
            
            # Check for deconvolution (Tangram)
            has_deconv = "tangram_ct_pred" in adata.obsm
            
            # Check for CCC (COMMOT)
            has_ccc = "commot_user_database" in adata.uns or "commot" in adata.uns
            
            # Check for spatial image
            has_image = False
            try:
                if "spatial" in adata.uns:
                    spatial_data = adata.uns["spatial"]
                    if isinstance(spatial_data, dict):
                        for sample_data in spatial_data.values():
                            if isinstance(sample_data, dict) and "images" in sample_data:
                                has_image = True
                                break
            except:
                pass
            
            # Detect dataset type
            dataset_type = None
            if has_deconv:
                dataset_type = "visium"
            elif has_ccc:
                dataset_type = "xenium"
            elif has_spatial:
                dataset_type = "unknown"  # Has spatial but can't determine type
            
            info = {
                "has_spatial_coordinates": has_spatial,
                "has_svg": has_svg,
                "has_precomputed_deg": has_deg,
                "has_deconvolution": has_deconv,
                "has_ccc": has_ccc,
                "has_spatial_image": has_image,
                "dataset_type": dataset_type
            }
            
            # Cache
            self._save_json(spatial_info_file, info)
            self.cache[cache_key] = info
            return info
            
        except Exception as e:
            logger.warning(f"Error detecting spatial features for {dataset_id}: {e}")
            return {
                "has_spatial_coordinates": False,
                "has_svg": False,
                "has_precomputed_deg": False,
                "has_deconvolution": False,
                "has_ccc": False,
                "dataset_type": None
            }
    
    def get_spatial_coordinates(self, dataset_id: str) -> Dict:
        """
        Get spatial coordinates (x, y) for spots/cells
        Checks multiple possible keys: spatial, spatial_xy, X_spatial
        Returns: {coordinates: [[x, y], ...], cell_ids: [...]}
        """
        cache_key = f"{dataset_id}:spatial_coords"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Check precomputed file
        spatial_file = self.get_precomputed_dir(dataset_id) / "spatial_coordinates.json"
        if spatial_file.exists():
            data = self._load_json(spatial_file)
            self.cache[cache_key] = data
            return data
        
        # Extract from h5ad
        logger.info(f"Extracting spatial coordinates from h5ad for {dataset_id}")
        adata = self._load_h5ad(dataset_id)
        
        # Try multiple possible keys
        spatial_key = None
        for key in ['spatial', 'spatial_xy', 'X_spatial']:
            if key in adata.obsm:
                spatial_key = key
                break
        
        if not spatial_key:
            raise ValueError(f"Spatial coordinates not found in dataset {dataset_id}. Checked: spatial, spatial_xy, X_spatial")
        
        spatial_coords = adata.obsm[spatial_key]
        
        # Handle 2D or 3D coordinates (take first 2 dimensions)
        if spatial_coords.shape[1] >= 2:
            coordinates = spatial_coords[:, :2].tolist()
        else:
            raise ValueError(f"Spatial coordinates must have at least 2 dimensions, got {spatial_coords.shape[1]}")
        
        cell_ids = adata.obs_names.tolist()
        
        data = {
            "coordinates": coordinates,
            "cell_ids": cell_ids,
            "spatial_key": spatial_key
        }
        
        # Cache
        self._save_json(spatial_file, data)
        self.cache[cache_key] = data
        
        return data
    
    def get_spatial_plot_data(
        self,
        dataset_id: str,
        genes: Optional[List[str]] = None,
        metadata_cols: Optional[List[str]] = None
    ) -> Dict:
        """
        Get combined spatial plot data (coordinates + gene expression + metadata)
        Returns everything needed for spatial visualization in one call
        """
        result = {}
        
        # Get spatial coordinates
        spatial_data = self.get_spatial_coordinates(dataset_id)
        result["coordinates"] = spatial_data["coordinates"]
        result["cell_ids"] = spatial_data["cell_ids"]
        
        # Get gene expression if requested
        if genes:
            result["genes"] = {}
            for gene in genes:
                try:
                    gene_data = self.get_gene_expression(dataset_id, gene)
                    result["genes"][gene] = gene_data["expression"]
                except ValueError as e:
                    logger.warning(f"Gene {gene} not found: {e}")
                    # Return empty list to satisfy response schema and signal missing gene
                    result["genes"][gene] = []
        
        # Get metadata if requested
        if metadata_cols:
            metadata = self.get_metadata(dataset_id, metadata_cols)
            result["metadata"] = metadata["data"]
        
        return result
    
    def get_svg_data(
        self,
        dataset_id: str,
        top_n: Optional[int] = None,
        min_score: Optional[float] = None
    ) -> Dict:
        """
        Get Spatially Variable Genes (SVG) data
        Extracts from ad.var columns: gft_score, svg_rank, cutoff_gft_score, pvalue, fdr
        
        Args:
            dataset_id: Dataset identifier
            top_n: Return only top N genes by rank (optional)
            min_score: Minimum gft_score threshold (optional)
        
        Returns: {
            "genes": [{gene, gft_score, svg_rank, cutoff_gft_score, pvalue, fdr}, ...],
            "total_genes": int,
            "filtered_count": int
        }
        """
        cache_key = f"{dataset_id}:svg"
        if cache_key in self.cache and top_n is None and min_score is None:
            return self.cache[cache_key]
        
        # Check precomputed file
        svg_file = self.get_precomputed_dir(dataset_id) / "svg_data.json"
        if svg_file.exists():
            svg_data = self._load_json(svg_file)
        else:
            # Extract from h5ad
            logger.info(f"Extracting SVG data from h5ad for {dataset_id}")
            adata = self._load_h5ad(dataset_id)
            
            # Check for SVG columns
            svg_cols = ['gft_score', 'svg_rank', 'cutoff_gft_score', 'pvalue', 'fdr']
            available_cols = [col for col in svg_cols if col in adata.var.columns]
            
            if not available_cols:
                raise ValueError(f"SVG data not found in dataset {dataset_id}. Expected columns: {svg_cols}")
            
            # Extract SVG data for all genes
            genes_list = []
            for gene in adata.var_names:
                gene_data = {"gene": gene}
                for col in svg_cols:
                    if col in adata.var.columns:
                        value = adata.var.loc[gene, col]
                        # Handle NaN values
                        if pd.isna(value):
                            gene_data[col] = None
                        elif isinstance(value, bool):
                            # Convert boolean to float (True -> 1.0, False -> 0.0)
                            gene_data[col] = float(value)
                        elif isinstance(value, (int, float, np.number)):
                            gene_data[col] = float(value)
                        elif isinstance(value, str):
                            # Try to convert string to float, otherwise keep as None
                            try:
                                gene_data[col] = float(value)
                            except (ValueError, TypeError):
                                gene_data[col] = None
                        else:
                            gene_data[col] = None
                    else:
                        gene_data[col] = None
                
                genes_list.append(gene_data)
            
            # Sort by svg_rank if available, otherwise by gft_score
            if 'svg_rank' in available_cols:
                genes_list.sort(key=lambda x: x['svg_rank'] if x['svg_rank'] is not None else float('inf'))
            elif 'gft_score' in available_cols:
                genes_list.sort(key=lambda x: x['gft_score'] if x['gft_score'] is not None else float('-inf'), reverse=True)
            
            svg_data = {
                "genes": genes_list,
                "total_genes": len(genes_list),
                "available_columns": available_cols
            }
            
            # Cache
            self._save_json(svg_file, svg_data)
        
        # Apply filters
        filtered_genes = svg_data["genes"].copy()
        
        if min_score is not None:
            filtered_genes = [g for g in filtered_genes if g.get('gft_score') is not None and g['gft_score'] >= min_score]
        
        if top_n is not None:
            filtered_genes = filtered_genes[:top_n]
        
        result = {
            "genes": filtered_genes,
            "total_genes": svg_data["total_genes"],
            "filtered_count": len(filtered_genes),
            "available_columns": svg_data.get("available_columns", [])
        }
        
        # Cache full result if no filters
        if top_n is None and min_score is None:
            self.cache[cache_key] = result
        
        return result
    
    def get_precomputed_deg(self, dataset_id: str, group: Optional[str] = None) -> Dict:
        """
        Get precomputed DEG results from ad.uns["rank_genes_groups"]
        
        Args:
            dataset_id: Dataset identifier
            group: Optional group name to filter (e.g., '0', '1', etc.)
        
        Returns: {
            "groups": [...],
            "params": {...},
            "genes": {group: [{gene, scores, logfoldchanges, pvals, pvals_adj}, ...]},
            "method": str
        }
        """
        cache_key = f"{dataset_id}:deg_precomputed"
        if cache_key in self.cache and group is None:
            return self.cache[cache_key]
        
        # Check precomputed file
        deg_file = self.get_precomputed_dir(dataset_id) / "deg_precomputed.json"
        if deg_file.exists():
            deg_data = self._load_json(deg_file)
        else:
            # Extract from h5ad
            logger.info(f"Extracting precomputed DEG from h5ad for {dataset_id}")
            adata = self._load_h5ad(dataset_id)
            
            if "rank_genes_groups" not in adata.uns:
                raise ValueError(f"Precomputed DEG not found in dataset {dataset_id}. Missing 'rank_genes_groups' in ad.uns")
            
            rank_genes = adata.uns["rank_genes_groups"]
            
            # Extract groups
            groups = rank_genes['names'].dtype.names if hasattr(rank_genes['names'], 'dtype') else list(rank_genes['names'].dtype.names) if hasattr(rank_genes['names'], 'dtype') else []
            if not groups:
                # Try alternative structure
                if isinstance(rank_genes['names'], np.ndarray):
                    groups = list(range(rank_genes['names'].shape[1]))
                else:
                    groups = list(rank_genes.get('names', {}).keys()) if isinstance(rank_genes.get('names'), dict) else []
            
            # Extract parameters
            params = rank_genes.get('params', {})
            
            # Extract gene data for each group
            genes_by_group = {}
            for group_name in groups:
                try:
                    # Get gene names, scores, logfoldchanges, pvals, pvals_adj
                    gene_names = rank_genes['names'][group_name] if isinstance(rank_genes['names'], np.ndarray) else rank_genes['names'].get(group_name, [])
                    scores = rank_genes.get('scores', {})
                    logfoldchanges = rank_genes.get('logfoldchanges', {})
                    pvals = rank_genes.get('pvals', {})
                    pvals_adj = rank_genes.get('pvals_adj', {})
                    
                    # Convert to list format
                    genes_list = []
                    for i, gene in enumerate(gene_names):
                        if pd.isna(gene) or gene == '':
                            continue
                        
                        gene_data = {"gene": str(gene)}
                        
                        if isinstance(scores, np.ndarray) and scores.shape[0] > i:
                            gene_data["score"] = float(scores[i]) if not pd.isna(scores[i]) else None
                        elif isinstance(scores, dict):
                            gene_data["score"] = float(scores.get(gene, None)) if scores.get(gene) is not None else None
                        
                        if isinstance(logfoldchanges, np.ndarray) and logfoldchanges.shape[0] > i:
                            gene_data["logfoldchange"] = float(logfoldchanges[i]) if not pd.isna(logfoldchanges[i]) else None
                        elif isinstance(logfoldchanges, dict):
                            gene_data["logfoldchange"] = float(logfoldchanges.get(gene, None)) if logfoldchanges.get(gene) is not None else None
                        
                        if isinstance(pvals, np.ndarray) and pvals.shape[0] > i:
                            gene_data["pval"] = float(pvals[i]) if not pd.isna(pvals[i]) else None
                        elif isinstance(pvals, dict):
                            gene_data["pval"] = float(pvals.get(gene, None)) if pvals.get(gene) is not None else None
                        
                        if isinstance(pvals_adj, np.ndarray) and pvals_adj.shape[0] > i:
                            gene_data["pval_adj"] = float(pvals_adj[i]) if not pd.isna(pvals_adj[i]) else None
                        elif isinstance(pvals_adj, dict):
                            gene_data["pval_adj"] = float(pvals_adj.get(gene, None)) if pvals_adj.get(gene) is not None else None
                        
                        genes_list.append(gene_data)
                    
                    genes_by_group[str(group_name)] = genes_list
                except Exception as e:
                    logger.warning(f"Error extracting DEG for group {group_name}: {e}")
                    genes_by_group[str(group_name)] = []
            
            deg_data = {
                "groups": [str(g) for g in groups],
                "params": {k: str(v) if not isinstance(v, (int, float, str, bool, type(None))) else v for k, v in params.items()},
                "genes": genes_by_group,
                "method": rank_genes.get('params', {}).get('method', 'unknown')
            }
            
            # Cache
            self._save_json(deg_file, deg_data)
        
        # Filter by group if requested
        if group is not None:
            if group not in deg_data["groups"]:
                raise ValueError(f"Group '{group}' not found. Available groups: {deg_data['groups']}")
            
            return {
                "groups": deg_data["groups"],
                "params": deg_data["params"],
                "genes": {group: deg_data["genes"].get(group, [])},
                "method": deg_data["method"]
            }
        
        # Cache full result
        self.cache[cache_key] = deg_data
        return deg_data
    
    def get_deconvolution_predictions(self, dataset_id: str) -> Dict:
        """
        Get deconvolution (Tangram) cell type predictions
        Extracts from ad.obsm["tangram_ct_pred"]
        
        Returns: {
            "predictions": {cell_type: [probabilities], ...},
            "cell_ids": [...],
            "cell_types": [...]
        }
        """
        cache_key = f"{dataset_id}:deconvolution"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Check precomputed file
        deconv_file = self.get_precomputed_dir(dataset_id) / "deconvolution.json"
        if deconv_file.exists():
            data = self._load_json(deconv_file)
            self.cache[cache_key] = data
            return data
        
        # Extract from h5ad
        logger.info(f"Extracting deconvolution predictions from h5ad for {dataset_id}")
        adata = self._load_h5ad(dataset_id)
        
        if "tangram_ct_pred" not in adata.obsm:
            raise ValueError(f"Deconvolution predictions not found in dataset {dataset_id}. Missing 'tangram_ct_pred' in ad.obsm")
        
        pred_matrix = adata.obsm["tangram_ct_pred"]
        
        # Convert to DataFrame for easier handling
        if isinstance(pred_matrix, np.ndarray):
            # Get cell type names from columns if available
            if hasattr(adata, 'uns') and 'tangram_ct_pred' in adata.uns:
                cell_types = adata.uns.get('tangram_ct_pred', {}).get('columns', [])
            else:
                # Try to infer from shape or use indices
                cell_types = [f"cell_type_{i}" for i in range(pred_matrix.shape[1])]
            
            pred_df = pd.DataFrame(pred_matrix, columns=cell_types)
        else:
            pred_df = pred_matrix
        
        # Convert to dict format
        predictions = {}
        for col in pred_df.columns:
            predictions[col] = pred_df[col].tolist()
        
        cell_ids = adata.obs_names.tolist()
        
        data = {
            "predictions": predictions,
            "cell_ids": cell_ids,
            "cell_types": list(pred_df.columns)
        }
        
        # Cache
        self._save_json(deconv_file, data)
        self.cache[cache_key] = data
        
        return data
    
    def get_ccc_interactions(self, dataset_id: str) -> Dict:
        """
        Get Cell-Cell Communication (CCC) interaction data
        Extracts from ad.uns['commot_user_database'] or ad.uns['commot']
        
        Returns: {
            "interactions": [...],
            "ligands": [...],
            "receptors": [...],
            "cell_types": [...],
            "data": {...}
        }
        """
        cache_key = f"{dataset_id}:ccc"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Check precomputed file
        ccc_file = self.get_precomputed_dir(dataset_id) / "ccc_interactions.json"
        if ccc_file.exists():
            data = self._load_json(ccc_file)
            # If interactions are empty but df_ligrec is available, reconstruct summary
            if (
                isinstance(data, dict)
                and not data.get("interactions")
                and isinstance(data.get("data"), dict)
            ):
                df_ligrec = data["data"].get("info", {}).get("df_ligrec", [])
                if isinstance(df_ligrec, list) and df_ligrec:
                    data["interactions"] = df_ligrec
                    data["ligands"] = sorted({r.get("ligand") for r in df_ligrec if r.get("ligand")})
                    data["receptors"] = sorted({r.get("receptor") for r in df_ligrec if r.get("receptor")})
                    # cell_types may be absent in this structure; leave empty to avoid guesswork
                    self._save_json(ccc_file, data)
            self.cache[cache_key] = data
            return data
        
        # Extract from h5ad
        logger.info(f"Extracting CCC interactions from h5ad for {dataset_id}")
        adata = self._load_h5ad(dataset_id)
        
        # Try multiple possible keys
        ccc_data = None
        if "commot_user_database" in adata.uns:
            ccc_data = adata.uns["commot_user_database"]
        elif "commot" in adata.uns:
            ccc_data = adata.uns["commot"]
        else:
            raise ValueError(f"CCC data not found in dataset {dataset_id}. Checked: commot_user_database, commot")
        
        # CCC data structure can vary, try to extract common fields
        result = {
            "data": ccc_data,
            "interactions": [],
            "ligands": [],
            "receptors": [],
            "cell_types": []
        }
        
        # Try to extract structured data if it's a dict or DataFrame
        if isinstance(ccc_data, dict):
            # Common COMMOT structure
            if "interactions" in ccc_data:
                result["interactions"] = ccc_data["interactions"]
            if "ligands" in ccc_data:
                result["ligands"] = ccc_data["ligands"]
            if "receptors" in ccc_data:
                result["receptors"] = ccc_data["receptors"]
            if "cell_types" in ccc_data:
                result["cell_types"] = ccc_data["cell_types"]
            
            # Fallback: df_ligrec inside info
            if not result["interactions"] and isinstance(ccc_data.get("info"), dict):
                df_ligrec = ccc_data["info"].get("df_ligrec", [])
                if isinstance(df_ligrec, list) and df_ligrec:
                    result["interactions"] = df_ligrec
                    result["ligands"] = sorted({r.get("ligand") for r in df_ligrec if r.get("ligand")})
                    result["receptors"] = sorted({r.get("receptor") for r in df_ligrec if r.get("receptor")})
        elif isinstance(ccc_data, pd.DataFrame):
            # If it's a DataFrame, convert to dict
            result["data"] = ccc_data.to_dict(orient='records')
            # Try to infer columns
            if "ligand" in ccc_data.columns:
                result["ligands"] = ccc_data["ligand"].unique().tolist()
            if "receptor" in ccc_data.columns:
                result["receptors"] = ccc_data["receptor"].unique().tolist()
        
        # Cache
        self._save_json(ccc_file, result)
        self.cache[cache_key] = result
        
        return result
    
    def get_spatial_image_info(self, dataset_id: str, sample_key: Optional[str] = None) -> Dict:
        """
        Get spatial image metadata and available samples
        
        Returns information about available spatial images including:
        - Available sample keys (slices/regions)
        - Image dimensions
        - Scale factors for coordinate transformation
        - Image format
        
        Args:
            dataset_id: Dataset identifier
            sample_key: Optional specific sample key (if None, returns info for all samples)
        
        Returns: {
            "samples": {
                "sample_key": {
                    "image_keys": ["hires", ...],
                    "image_shapes": {"hires": [height, width, channels]},
                    "scalefactors": {...},
                    "image_paths": {"hires": "path/to/image.png"}
                }
            },
            "default_sample": "sample_key"
        }
        """
        cache_key = f"{dataset_id}:spatial_image_info"
        if cache_key in self.cache and sample_key is None:
            return self.cache[cache_key]
        
        # Check precomputed file
        image_info_file = self.get_precomputed_dir(dataset_id) / "spatial_image_info.json"
        if image_info_file.exists():
            info = self._load_json(image_info_file)
            if sample_key:
                if sample_key in info.get("samples", {}):
                    return {
                        "samples": {sample_key: info["samples"][sample_key]},
                        "default_sample": sample_key
                    }
                else:
                    raise ValueError(f"Sample key '{sample_key}' not found. Available: {list(info.get('samples', {}).keys())}")
            self.cache[cache_key] = info
            return info
        
        # Extract from h5ad
        logger.info(f"Extracting spatial image info from h5ad for {dataset_id}")
        adata = self._load_h5ad(dataset_id)
        
        if "spatial" not in adata.uns:
            raise ValueError(f"Spatial image data not found in dataset {dataset_id}")
        
        spatial_data = adata.uns["spatial"]
        if not isinstance(spatial_data, dict):
            raise ValueError(f"Unexpected spatial data format in dataset {dataset_id}")
        
        samples_info = {}
        default_sample = None
        
        for key, sample_data in spatial_data.items():
            if not isinstance(sample_data, dict):
                continue
            
            sample_info = {
                "image_keys": [],
                "image_shapes": {},
                "scalefactors": {},
                "image_paths": {}
            }
            
            # Extract images
            if "images" in sample_data:
                images = sample_data["images"]
                if isinstance(images, dict):
                    for img_key, img_array in images.items():
                        if isinstance(img_array, np.ndarray):
                            sample_info["image_keys"].append(img_key)
                            sample_info["image_shapes"][img_key] = list(img_array.shape)
                            
                            # Save image to cache
                            image_path = self._save_spatial_image(dataset_id, key, img_key, img_array)
                            # Store relative path from precomputed directory
                            rel_path = image_path.relative_to(self.get_precomputed_dir(dataset_id))
                            sample_info["image_paths"][img_key] = str(rel_path)
            
            # Extract scale factors
            if "scalefactors" in sample_data:
                scalefactors = sample_data["scalefactors"]
                if isinstance(scalefactors, dict):
                    # Convert to JSON-serializable
                    sample_info["scalefactors"] = {
                        k: float(v) if isinstance(v, (int, float, np.number)) else v
                        for k, v in scalefactors.items()
                    }
            
            if sample_info["image_keys"]:
                samples_info[key] = sample_info
                if default_sample is None:
                    default_sample = key
        
        if not samples_info:
            raise ValueError(f"No spatial images found in dataset {dataset_id}")
        
        result = {
            "samples": samples_info,
            "default_sample": default_sample
        }
        
        # Cache
        self._save_json(image_info_file, result)
        self.cache[cache_key] = result
        
        # Filter by sample_key if requested
        if sample_key:
            if sample_key not in samples_info:
                raise ValueError(f"Sample key '{sample_key}' not found. Available: {list(samples_info.keys())}")
            return {
                "samples": {sample_key: samples_info[sample_key]},
                "default_sample": sample_key
            }
        
        return result
    
    def _save_spatial_image(self, dataset_id: str, sample_key: str, image_key: str, image_array: np.ndarray) -> Path:
        """
        Save spatial image array to PNG file
        
        Args:
            dataset_id: Dataset identifier
            sample_key: Sample/slice key
            image_key: Image key (e.g., 'hires')
            image_array: Image array (numpy array)
        
        Returns: Path to saved image file
        """
        image_dir = self.get_precomputed_dir(dataset_id) / "spatial_images" / sample_key
        image_dir.mkdir(parents=True, exist_ok=True)
        
        image_path = image_dir / f"{image_key}.png"
        
        # Skip if already exists
        if image_path.exists():
            return image_path
        
        # Normalize image array
        if image_array.dtype != np.uint8:
            # Normalize to 0-255 range
            if image_array.max() > 1.0:
                # Assume already in 0-255 range, just convert type
                image_array = image_array.astype(np.uint8)
            else:
                # Normalize from 0-1 to 0-255
                image_array = (image_array * 255).astype(np.uint8)
        else:
            image_array = image_array.copy()
        
        # Handle different image shapes
        if len(image_array.shape) == 2:
            # Grayscale
            img = Image.fromarray(image_array, mode='L')
        elif len(image_array.shape) == 3:
            if image_array.shape[2] == 3:
                # RGB
                img = Image.fromarray(image_array, mode='RGB')
            elif image_array.shape[2] == 4:
                # RGBA
                img = Image.fromarray(image_array, mode='RGBA')
            else:
                # Take first 3 channels
                img = Image.fromarray(image_array[:, :, :3], mode='RGB')
        else:
            raise ValueError(f"Unsupported image shape: {image_array.shape}")
        
        # Save as PNG
        img.save(image_path, "PNG")
        logger.info(f"Saved spatial image to {image_path}")
        
        return image_path
    
    def get_spatial_coordinates_transformed(
        self,
        dataset_id: str,
        sample_key: Optional[str] = None,
        image_key: str = "hires"
    ) -> Dict:
        """
        Get spatial coordinates transformed to match image space
        
        This applies the scale factor transformation so coordinates align with the image.
        Formula: image_coord = fullres_coord * tissue_hires_scalef
        
        Args:
            dataset_id: Dataset identifier
            sample_key: Optional sample key (uses default if None)
            image_key: Image key to transform for (default: "hires")
        
        Returns: {
            "coordinates": [[x, y], ...],  # Transformed to image space
            "cell_ids": [...],
            "sample_key": "sample_key",
            "image_key": "hires",
            "scale_factor": float,
            "original_coordinates": [[x, y], ...]  # Original full-res coords
        }
        """
        # Get spatial coordinates (full resolution)
        spatial_data = self.get_spatial_coordinates(dataset_id)
        original_coords = np.array(spatial_data["coordinates"])
        
        # Get image info to find scale factor
        image_info = self.get_spatial_image_info(dataset_id, sample_key)
        
        if sample_key is None:
            sample_key = image_info["default_sample"]
        
        if sample_key not in image_info["samples"]:
            raise ValueError(f"Sample key '{sample_key}' not found")
        
        sample_info = image_info["samples"][sample_key]
        scalefactors = sample_info.get("scalefactors", {})
        
        # Get scale factor for the image key
        # Common keys: tissue_hires_scalef, tissue_lowres_scalef
        scale_factor = 1.0
        if image_key == "hires":
            scale_factor = scalefactors.get("tissue_hires_scalef", 1.0)
        elif image_key == "lowres":
            scale_factor = scalefactors.get("tissue_lowres_scalef", 1.0)
        else:
            # Try to find a matching scale factor
            scale_key = f"tissue_{image_key}_scalef"
            scale_factor = scalefactors.get(scale_key, scalefactors.get("tissue_hires_scalef", 1.0))
        
        # Transform coordinates
        transformed_coords = (original_coords * scale_factor).tolist()
        
        return {
            "coordinates": transformed_coords,
            "cell_ids": spatial_data["cell_ids"],
            "sample_key": sample_key,
            "image_key": image_key,
            "scale_factor": float(scale_factor),
            "original_coordinates": original_coords.tolist(),
            "image_shape": sample_info.get("image_shapes", {}).get(image_key, None)
        }

    def get_ccc_network_from_csv(self, dataset_id: str, pathway: Optional[str] = None,
                                   min_prob: float = 0, source_cluster: Optional[str] = None,
                                   target_cluster: Optional[str] = None) -> Dict:
        """
        Get CCC network data from precomputed CSV file (CellChat format)

        The CSV should have columns: source, target, ligand, receptor, prob, pval,
        interaction_name, interaction_name_2, pathway_name, annotation, evidence

        Args:
            dataset_id: Dataset identifier
            pathway: Filter by pathway name (optional)
            min_prob: Minimum probability threshold
            source_cluster: Filter by source cluster (optional)
            target_cluster: Filter by target cluster (optional)

        Returns: {
            "nodes": [{"id": "C0", "label": "C0", "group": "C0"}, ...],
            "edges": [{"from": "C0", "to": "C1", "ligand": "...", "receptor": "...",
                       "pathway": "...", "prob": 0.1, ...}, ...]
        }
        """
        cache_key = f"{dataset_id}:ccc_csv:{pathway}:{min_prob}:{source_cluster}:{target_cluster}"
        if cache_key in self.cache:
            return self.cache[cache_key]

        # Check for CSV file
        ccc_csv_file = self.get_precomputed_dir(dataset_id) / "ccc" / "communications.csv"

        if not ccc_csv_file.exists():
            raise ValueError(f"CCC CSV file not found at {ccc_csv_file}")

        # Load CSV
        logger.info(f"Loading CCC data from CSV for {dataset_id}")
        df = pd.read_csv(ccc_csv_file)

        # Apply filters
        if pathway and pathway != 'all':
            df = df[df['pathway_name'] == pathway]

        if min_prob > 0:
            df = df[df['prob'] >= min_prob]

        if source_cluster:
            df = df[df['source'] == source_cluster]

        if target_cluster:
            df = df[df['target'] == target_cluster]

        # Build nodes and edges
        nodes_set = set()
        edges = []

        for _, row in df.iterrows():
            source = row['source']
            target = row['target']
            nodes_set.add(source)
            nodes_set.add(target)

            edges.append({
                "from": source,
                "to": target,
                "ligand": row.get('ligand', ''),
                "receptor": row.get('receptor', ''),
                "pathway": row.get('pathway_name', ''),
                "prob": float(row.get('prob', 0)) if pd.notna(row.get('prob')) else 0,
                "pval": float(row.get('pval', 1)) if pd.notna(row.get('pval')) else 1,
                "interaction_name": row.get('interaction_name', ''),
                "annotation": row.get('annotation', ''),
                "evidence": row.get('evidence', '')
            })

        nodes = [{"id": n, "label": n, "group": n} for n in sorted(nodes_set)]

        # Get unique pathways for filtering
        pathways = sorted(df['pathway_name'].dropna().unique().tolist())

        result = {
            "nodes": nodes,
            "edges": edges,
            "pathways": pathways,
            "total_interactions": len(edges)
        }

        self.cache[cache_key] = result
        return result

    def get_regulon_network_from_csv(self, dataset_id: str, cluster: str = '0',
                                      tf: Optional[str] = None, max_targets: int = 100) -> Dict:
        """
        Get Regulon TF-target network data from precomputed CSV files

        CSV files should be named: Cluster_cluster_{n}_regulon_network.csv
        Each CSV has columns: TF, Target

        Args:
            dataset_id: Dataset identifier
            cluster: Cluster number (default '0')
            tf: Filter by specific transcription factor (optional)
            max_targets: Maximum number of targets per TF

        Returns: {
            "nodes": [{"id": "TF1", "label": "TF1", "type": "tf"},
                      {"id": "GENE1", "label": "GENE1", "type": "target"}, ...],
            "edges": [{"tf": "TF1", "target": "GENE1"}, ...],
            "tfs": ["TF1", "TF2", ...],
            "available_clusters": ["0", "1", ...]
        }
        """
        cache_key = f"{dataset_id}:regulon:{cluster}:{tf}:{max_targets}"
        if cache_key in self.cache:
            return self.cache[cache_key]

        regulon_dir = self.get_precomputed_dir(dataset_id) / "regulon"

        if not regulon_dir.exists():
            raise ValueError(f"Regulon directory not found at {regulon_dir}")

        # Find available cluster files
        cluster_files = list(regulon_dir.glob("Cluster_cluster_*_regulon_network.csv"))
        available_clusters = []
        for f in cluster_files:
            # Extract cluster number from filename
            match = f.stem.replace("Cluster_cluster_", "").replace("_regulon_network", "")
            available_clusters.append(match)

        available_clusters = sorted(available_clusters, key=lambda x: int(x) if x.isdigit() else x)

        # Load specific cluster file
        cluster_file = regulon_dir / f"Cluster_cluster_{cluster}_regulon_network.csv"

        if not cluster_file.exists():
            # Try without the "Cluster_" prefix
            cluster_file = regulon_dir / f"cluster_{cluster}_regulon_network.csv"

        if not cluster_file.exists():
            raise ValueError(f"Regulon file not found for cluster {cluster}. Available: {available_clusters}")

        logger.info(f"Loading regulon data from CSV for {dataset_id}, cluster {cluster}")
        df = pd.read_csv(cluster_file)

        # Filter by TF if specified
        if tf and tf != 'all':
            df = df[df['TF'] == tf]

        # Build nodes and edges
        tf_set = set()
        target_set = set()
        edges = []

        # Count targets per TF for limiting
        tf_target_counts = {}

        for _, row in df.iterrows():
            tf_name = row['TF']
            target_name = row['Target']

            # Limit targets per TF
            tf_target_counts[tf_name] = tf_target_counts.get(tf_name, 0) + 1
            if tf_target_counts[tf_name] > max_targets:
                continue

            tf_set.add(tf_name)
            target_set.add(target_name)
            edges.append({
                "tf": tf_name,
                "target": target_name
            })

        # Build nodes list
        nodes = []
        for tf_name in sorted(tf_set):
            nodes.append({"id": tf_name, "label": tf_name, "type": "tf"})
        for target_name in sorted(target_set):
            if target_name not in tf_set:  # Avoid duplicates if a TF is also a target
                nodes.append({"id": target_name, "label": target_name, "type": "target"})

        # Get TF statistics (total target count before limiting)
        tf_stats = df.groupby('TF').size().to_dict()
        tfs = sorted(tf_stats.keys(), key=lambda x: tf_stats[x], reverse=True)

        result = {
            "nodes": nodes,
            "edges": edges,
            "tfs": tfs,
            "tf_stats": tf_stats,
            "available_clusters": available_clusters,
            "cluster": cluster,
            "total_edges": len(edges)
        }

        self.cache[cache_key] = result
        return result

