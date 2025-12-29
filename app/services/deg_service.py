"""
Differential Expression Gene (DEG) Analysis Service
Performs DEG analysis between two conditions/groups
"""
import logging
from typing import Optional, List, Dict, Any
from pathlib import Path
import json

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class DEGService:
    """Service for differential expression analysis"""
    
    def __init__(self, h5ad_base_path: str = "h5ad", cache_base_path: str = "h5ad/deg_cache"):
        self.h5ad_base_path = Path(h5ad_base_path)
        self.cache_base_path = Path(cache_base_path)
        self.cache_base_path.mkdir(parents=True, exist_ok=True)
        self._metadata_cache = {}  # Cache for metadata-only DataFrames
    
    def _generate_cache_key(self, params: Dict) -> str:
        """Generate unique cache key for DEG analysis"""
        import hashlib
        key_string = json.dumps(params, sort_keys=True)
        return hashlib.md5(key_string.encode()).hexdigest()
    
    def _load_h5ad(self, dataset_id: str):
        """Load h5ad file"""
        import anndata
        h5ad_path = self.h5ad_base_path / "raw" / f"{dataset_id}.h5ad"
        if not h5ad_path.exists():
            raise FileNotFoundError(f"H5AD file not found: {h5ad_path}")
        return anndata.read_h5ad(h5ad_path)
    
    def _load_metadata_only(self, dataset_id: str) -> pd.DataFrame:
        """
        Load only metadata (obs) from h5ad file or precomputed cache.
        Uses caching to avoid repeated file reads.
        """
        # Check cache first
        if dataset_id in self._metadata_cache:
            logger.debug(f"Using cached metadata for {dataset_id}")
            return self._metadata_cache[dataset_id]
        
        logger.info(f"Loading metadata for {dataset_id}")
        
        # Try to load from precomputed metadata.json first (fastest)
        precomputed_dir = self.h5ad_base_path / "precomputed" / dataset_id
        metadata_file = precomputed_dir / "metadata.json"
        
        if metadata_file.exists():
            logger.info(f"Loading metadata from precomputed cache: {metadata_file}")
            with open(metadata_file, 'r') as f:
                metadata_json = json.load(f)
            metadata_df = pd.DataFrame(metadata_json['data'])
            metadata_df.index = metadata_json['cell_ids']
        else:
            # Fall back to loading from h5ad (slower but works)
            logger.info(f"Loading metadata from h5ad file (no precomputed cache found)")
            adata = self._load_h5ad(dataset_id)
            metadata_df = adata.obs.copy()
        
        # Cache it
        self._metadata_cache[dataset_id] = metadata_df
        logger.info(f"Cached metadata for {dataset_id} ({len(metadata_df)} cells, {len(metadata_df.columns)} columns)")
        
        return metadata_df
    
    def _filter_cells(self, adata, cell_type: Optional[str] = None,
                     cell_types: Optional[List[str]] = None,
                     metadata_filters: Optional[Dict[str, str]] = None):
        """Filter cells based on criteria"""
        mask = np.ones(adata.n_obs, dtype=bool)
        
        applicable_cell_types: Optional[List[str]] = None
        if cell_types:
            cleaned = [ct for ct in cell_types if ct]
            applicable_cell_types = cleaned or None
        if applicable_cell_types is None and cell_type:
            applicable_cell_types = [cell_type]
        
        if applicable_cell_types and 'cell_type' in adata.obs.columns:
            mask &= adata.obs['cell_type'].isin(applicable_cell_types)
        
        if metadata_filters:
            for key, value in metadata_filters.items():
                if key in adata.obs.columns:
                    mask &= (adata.obs[key] == value)
        
        return adata[mask].copy()
    
    def _normalize_cell_types(self, *values) -> Optional[List[str]]:
        """Normalize cell type inputs (lists or comma-separated strings)"""
        normalized: List[str] = []
        seen = set()
        
        for value in values:
            if not value:
                continue
            if isinstance(value, str):
                candidates = [v.strip() for v in value.split(",") if v.strip()]
            else:
                candidates = [v.strip() for v in value if isinstance(v, str) and v.strip()]
            
            for candidate in candidates:
                if candidate and candidate not in seen:
                    seen.add(candidate)
                    normalized.append(candidate)
        
        return normalized or None
    
    def _should_log_transform(self, data: np.ndarray) -> bool:
        """Determine if log1p transformation should be applied"""
        if data.size == 0:
            return False
        max_val = float(np.max(data))
        # Heuristic: if maximum value is large, assume raw counts and log-transform
        return max_val > 20
    
    def _perform_ttest(self, group1_data: np.ndarray, group2_data: np.ndarray) -> Dict:
        """Perform t-test between two groups"""
        from scipy import stats
        
        # Handle sparse matrices
        if hasattr(group1_data, 'toarray'):
            group1_data = group1_data.toarray()
        if hasattr(group2_data, 'toarray'):
            group2_data = group2_data.toarray()
        
        # Calculate statistics
        mean1 = np.mean(group1_data, axis=0)
        mean2 = np.mean(group2_data, axis=0)
        
        log2fc = np.log2((mean1 + 1e-9) / (mean2 + 1e-9))
        
        # Perform t-test
        t_stat, p_values = stats.ttest_ind(group1_data, group2_data, axis=0)
        
        # Replace NaN with 1.0 (no difference)
        p_values = np.nan_to_num(p_values, nan=1.0)
        
        return {
            'mean_group1': mean1,
            'mean_group2': mean2,
            'log2_fold_change': log2fc,
            't_statistic': t_stat,
            'p_value': p_values
        }
    
    def _adjust_pvalues(self, p_values: np.ndarray, method: str = 'fdr_bh') -> np.ndarray:
        """Adjust p-values for multiple testing"""
        from statsmodels.stats.multitest import multipletests
        
        # Handle edge cases
        if len(p_values) == 0:
            return np.array([])
        
        _, p_adjusted, _, _ = multipletests(p_values, method=method)
        return p_adjusted
    
    def analyze_deg_between_datasets(
        self,
        dataset_id1: str,
        dataset_id2: str,
        cell_type: Optional[str] = None,
        cell_type2: Optional[str] = None,
        cell_types_group1: Optional[List[str]] = None,
        cell_types_group2: Optional[List[str]] = None,
        metadata_filters1: Optional[Dict[str, str]] = None,
        metadata_filters2: Optional[Dict[str, str]] = None,
        min_pct: float = 0.1,
        logfc_threshold: float = 0.25,
        p_value_threshold: float = 0.05,
        top_n: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Perform DEG analysis between two datasets
        
        Args:
            dataset_id1: First dataset ID (condition 1)
            dataset_id2: Second dataset ID (condition 2)
            cell_type: Optional cell type to filter (e.g., "Neuron")
            metadata_filters1: Additional filters for dataset 1
            metadata_filters2: Additional filters for dataset 2
            min_pct: Minimum fraction of cells expressing the gene in either group
            logfc_threshold: Minimum log2 fold change threshold
            p_value_threshold: P-value threshold for significance
            top_n: Return only top N genes (by p-value)
        
        Returns:
            Dictionary with DEG results
        """
        logger.info(f"Starting DEG analysis: {dataset_id1} vs {dataset_id2}")
        
        # Check cache
        cache_params = {
            'dataset_id1': dataset_id1,
            'dataset_id2': dataset_id2,
            'cell_type': cell_type,
            'cell_type2': cell_type2,
            'cell_types_group1': cell_types_group1,
            'cell_types_group2': cell_types_group2,
            'metadata_filters1': metadata_filters1,
            'metadata_filters2': metadata_filters2,
            'min_pct': min_pct,
            'logfc_threshold': logfc_threshold,
            'p_value_threshold': p_value_threshold
        }
        cache_key = self._generate_cache_key(cache_params)
        cache_file = self.cache_base_path / f"{cache_key}.json"
        
        if cache_file.exists():
            logger.info(f"Loading from cache: {cache_key}")
            with open(cache_file, 'r') as f:
                result = json.load(f)
                if top_n:
                    result['genes'] = result['genes'][:top_n]
                return result
        
        # Load datasets
        logger.info(f"Loading dataset: {dataset_id1}")
        adata1 = self._load_h5ad(dataset_id1)
        logger.info(f"Loading dataset: {dataset_id2}")
        adata2 = self._load_h5ad(dataset_id2)
        
        # Determine cell types per group
        group1_cell_types = self._normalize_cell_types(cell_types_group1, cell_type)
        group2_cell_types = self._normalize_cell_types(cell_types_group2, cell_type2)
        if not group2_cell_types:
            group2_cell_types = self._normalize_cell_types(cell_type2, cell_type)
        
        # Filter cells
        if group1_cell_types or cell_type or metadata_filters1:
            logger.info(f"Filtering cells in {dataset_id1}")
            adata1 = self._filter_cells(
                adata1,
                cell_type=cell_type if not group1_cell_types else None,
                cell_types=group1_cell_types,
                metadata_filters=metadata_filters1
            )
        if group2_cell_types or cell_type or cell_type2 or metadata_filters2:
            logger.info(f"Filtering cells in {dataset_id2}")
            adata2 = self._filter_cells(
                adata2,
                cell_type=(cell_type2 or cell_type) if not group2_cell_types else None,
                cell_types=group2_cell_types,
                metadata_filters=metadata_filters2
            )
        
        # Check if we have cells left
        if adata1.n_obs == 0 or adata2.n_obs == 0:
            raise ValueError("No cells remaining after filtering")
        
        logger.info(f"Cells after filtering: {adata1.n_obs} vs {adata2.n_obs}")
        
        # Ensure same genes
        common_genes = adata1.var_names.intersection(adata2.var_names)
        adata1 = adata1[:, common_genes]
        adata2 = adata2[:, common_genes]
        
        logger.info(f"Common genes: {len(common_genes)}")
        
        # Get expression data
        X1 = adata1.X.toarray() if hasattr(adata1.X, 'toarray') else np.array(adata1.X)
        X2 = adata2.X.toarray() if hasattr(adata2.X, 'toarray') else np.array(adata2.X)
        
        # Calculate percentage of cells expressing each gene (using raw counts)
        pct1 = np.sum(X1 > 0, axis=0) / X1.shape[0]
        pct2 = np.sum(X2 > 0, axis=0) / X2.shape[0]
        
        # Apply log1p transform if needed
        if self._should_log_transform(X1) or self._should_log_transform(X2):
            logger.info("Applying log1p transformation to expression matrices")
            X1 = np.log1p(np.clip(X1, a_min=0, a_max=None))
            X2 = np.log1p(np.clip(X2, a_min=0, a_max=None))
        
        # Filter by minimum percentage
        expressed_mask = (pct1 >= min_pct) | (pct2 >= min_pct)
        
        # Perform statistical test
        logger.info("Performing statistical test...")
        stats_results = self._perform_ttest(X1[:, expressed_mask], X2[:, expressed_mask])
        
        # Adjust p-values
        logger.info("Adjusting p-values...")
        p_adjusted = self._adjust_pvalues(stats_results['p_value'])
        
        # Compile results
        expressed_genes = common_genes[expressed_mask]
        results_df = pd.DataFrame({
            'gene': expressed_genes,
            'mean_group1': stats_results['mean_group1'],
            'mean_group2': stats_results['mean_group2'],
            'log2_fold_change': stats_results['log2_fold_change'],
            'pct_group1': pct1[expressed_mask],
            'pct_group2': pct2[expressed_mask],
            'p_value': stats_results['p_value'],
            'p_value_adj': p_adjusted
        })
        
        # Filter by significance
        significant = (
            (np.abs(results_df['log2_fold_change']) >= logfc_threshold) &
            (results_df['p_value_adj'] < p_value_threshold)
        )
        
        results_df['significant'] = significant
        
        # Sort by absolute log2 fold change and p-value
        results_df['abs_log2fc'] = np.abs(results_df['log2_fold_change'])
        results_df = results_df.sort_values(['significant', 'p_value_adj', 'abs_log2fc'], 
                                            ascending=[False, True, False])
        
        # Prepare output
        genes_list = []
        for _, row in results_df.iterrows():
            genes_list.append({
                'gene': row['gene'],
                'log2_fold_change': float(row['log2_fold_change']),
                'mean_expr_group1': float(row['mean_group1']),
                'mean_expr_group2': float(row['mean_group2']),
                'pct_group1': float(row['pct_group1']),
                'pct_group2': float(row['pct_group2']),
                'p_value': float(row['p_value']),
                'p_value_adj': float(row['p_value_adj']),
                'significant': bool(row['significant'])
            })
        
        result = {
            'comparison': {
                'group1': {
                    'dataset_id': dataset_id1,
                    'n_cells': int(adata1.n_obs),
                    'cell_type': cell_type,
                    'cell_types': group1_cell_types,
                    'filters': metadata_filters1
                },
                'group2': {
                    'dataset_id': dataset_id2,
                    'n_cells': int(adata2.n_obs),
                    'cell_type': cell_type2 or cell_type,
                    'cell_types': group2_cell_types,
                    'filters': metadata_filters2
                }
            },
            'parameters': {
                'min_pct': min_pct,
                'logfc_threshold': logfc_threshold,
                'p_value_threshold': p_value_threshold
            },
            'summary': {
                'total_genes_tested': int(np.sum(expressed_mask)),
                'significant_genes': int(np.sum(significant)),
                'upregulated_in_group1': int(np.sum((results_df['log2_fold_change'] > logfc_threshold) & significant)),
                'upregulated_in_group2': int(np.sum((results_df['log2_fold_change'] < -logfc_threshold) & significant))
            },
            'genes': genes_list
        }
        
        # Cache results
        logger.info(f"Caching results: {cache_key}")
        with open(cache_file, 'w') as f:
            json.dump(result, f)
        
        # Apply top_n filter if requested
        if top_n:
            result['genes'] = result['genes'][:top_n]
        
        logger.info("DEG analysis completed")
        return result
    
    def analyze_deg_within_dataset(
        self,
        dataset_id: str,
        group1_filters: Dict[str, str],
        group2_filters: Dict[str, str],
        min_pct: float = 0.1,
        logfc_threshold: float = 0.25,
        p_value_threshold: float = 0.05,
        top_n: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Perform DEG analysis within a single dataset between two groups
        
        Args:
            dataset_id: Dataset ID
            group1_filters: Filters for group 1 (e.g., {"cell_type": "Neuron", "condition": "Control"})
            group2_filters: Filters for group 2 (e.g., {"cell_type": "Neuron", "condition": "Disease"})
            min_pct: Minimum fraction of cells expressing the gene
            logfc_threshold: Minimum log2 fold change threshold
            p_value_threshold: P-value threshold
            top_n: Return only top N genes
        
        Returns:
            Dictionary with DEG results
        """
        logger.info(f"Starting within-dataset DEG analysis: {dataset_id}")
        
        # Load dataset
        adata = self._load_h5ad(dataset_id)
        
        # Filter for each group
        logger.info(f"Filtering group 1: {group1_filters}")
        adata1 = self._filter_cells(adata, metadata_filters=group1_filters)
        logger.info(f"Filtering group 2: {group2_filters}")
        adata2 = self._filter_cells(adata, metadata_filters=group2_filters)
        
        if adata1.n_obs == 0 or adata2.n_obs == 0:
            raise ValueError("No cells remaining after filtering")
        
        logger.info(f"Cells: {adata1.n_obs} vs {adata2.n_obs}")
        
        # Perform analysis (similar to between-dataset analysis)
        X1 = adata1.X.toarray() if hasattr(adata1.X, 'toarray') else np.array(adata1.X)
        X2 = adata2.X.toarray() if hasattr(adata2.X, 'toarray') else np.array(adata2.X)
        
        pct1 = np.sum(X1 > 0, axis=0) / X1.shape[0]
        pct2 = np.sum(X2 > 0, axis=0) / X2.shape[0]
        
        if self._should_log_transform(X1) or self._should_log_transform(X2):
            logger.info("Applying log1p transformation to expression matrices (within-dataset)")
            X1 = np.log1p(np.clip(X1, a_min=0, a_max=None))
            X2 = np.log1p(np.clip(X2, a_min=0, a_max=None))
        
        expressed_mask = (pct1 >= min_pct) | (pct2 >= min_pct)
        
        logger.info("Performing statistical test...")
        stats_results = self._perform_ttest(X1[:, expressed_mask], X2[:, expressed_mask])
        
        logger.info("Adjusting p-values...")
        p_adjusted = self._adjust_pvalues(stats_results['p_value'])
        
        # Compile results
        expressed_genes = adata.var_names[expressed_mask]
        results_df = pd.DataFrame({
            'gene': expressed_genes,
            'mean_group1': stats_results['mean_group1'],
            'mean_group2': stats_results['mean_group2'],
            'log2_fold_change': stats_results['log2_fold_change'],
            'pct_group1': pct1[expressed_mask],
            'pct_group2': pct2[expressed_mask],
            'p_value': stats_results['p_value'],
            'p_value_adj': p_adjusted
        })
        
        significant = (
            (np.abs(results_df['log2_fold_change']) >= logfc_threshold) &
            (results_df['p_value_adj'] < p_value_threshold)
        )
        
        results_df['significant'] = significant
        results_df['abs_log2fc'] = np.abs(results_df['log2_fold_change'])
        results_df = results_df.sort_values(['significant', 'p_value_adj', 'abs_log2fc'], 
                                            ascending=[False, True, False])
        
        genes_list = []
        for _, row in results_df.iterrows():
            genes_list.append({
                'gene': row['gene'],
                'log2_fold_change': float(row['log2_fold_change']),
                'mean_expr_group1': float(row['mean_group1']),
                'mean_expr_group2': float(row['mean_group2']),
                'pct_group1': float(row['pct_group1']),
                'pct_group2': float(row['pct_group2']),
                'p_value': float(row['p_value']),
                'p_value_adj': float(row['p_value_adj']),
                'significant': bool(row['significant'])
            })
        
        result = {
            'comparison': {
                'dataset_id': dataset_id,
                'group1': {
                    'filters': group1_filters,
                    'n_cells': int(adata1.n_obs)
                },
                'group2': {
                    'filters': group2_filters,
                    'n_cells': int(adata2.n_obs)
                }
            },
            'parameters': {
                'min_pct': min_pct,
                'logfc_threshold': logfc_threshold,
                'p_value_threshold': p_value_threshold
            },
            'summary': {
                'total_genes_tested': int(np.sum(expressed_mask)),
                'significant_genes': int(np.sum(significant)),
                'upregulated_in_group1': int(np.sum((results_df['log2_fold_change'] > logfc_threshold) & significant)),
                'upregulated_in_group2': int(np.sum((results_df['log2_fold_change'] < -logfc_threshold) & significant))
            },
            'genes': genes_list
        }
        
        if top_n:
            result['genes'] = result['genes'][:top_n]
        
        logger.info("DEG analysis completed")
        return result
    
    def _select_cells_multi_filter(self, adata, filters: List[Dict[str, Any]]) -> np.ndarray:
        """
        Select cells matching all filter conditions using multi-filter logic.
        
        Args:
            adata: AnnData object
            filters: List of {column, values} dicts
            
        Returns:
            Boolean mask of matching cells
            
        Logic:
            - Multiple filters use AND logic
            - Multiple values within a filter use OR logic
            
        Example:
            filters = [
                {"column": "cell_type", "values": ["Astrocyte", "Oligodendrocyte"]},
                {"column": "sex", "values": ["Male"]}
            ]
            Result: (cell_type == "Astrocyte" OR cell_type == "Oligodendrocyte") AND (sex == "Male")
        """
        mask = np.ones(adata.n_obs, dtype=bool)
        
        for f in filters:
            column = f.get('column')
            values = f.get('values', [])
            
            if not column or not values:
                continue
                
            if column not in adata.obs.columns:
                raise ValueError(f"Column '{column}' not found in dataset. "
                               f"Available columns: {list(adata.obs.columns)}")
            
            # OR logic within values, AND logic between filters
            filter_mask = adata.obs[column].isin(values)
            mask = mask & filter_mask
        
        return mask
    
    def _select_cells_multi_filter_from_metadata(self, metadata_df: pd.DataFrame, filters: List[Dict[str, Any]]) -> np.ndarray:
        """
        Select cells from metadata DataFrame using multi-filter logic.
        Optimized version that works with metadata-only DataFrame.
        """
        mask = np.ones(len(metadata_df), dtype=bool)
        
        for f in filters:
            column = f.get('column')
            values = f.get('values', [])
            
            if not column or not values:
                continue
                
            if column not in metadata_df.columns:
                raise ValueError(f"Column '{column}' not found in dataset. "
                               f"Available columns: {list(metadata_df.columns)}")
            
            # OR logic within values, AND logic between filters
            filter_mask = metadata_df[column].isin(values)
            mask = mask & filter_mask
        
        return mask
    
    def preview_cell_counts(
        self,
        dataset_id: str,
        group1_filters: List[Dict[str, Any]],
        group2_filters: List[Dict[str, Any]]
    ) -> Dict:
        """
        Preview cell counts for filter criteria without running full DEG analysis.
        
        OPTIMIZED: Only loads metadata (not expression data) for fast response.
        
        Args:
            dataset_id: Dataset identifier
            group1_filters: List of filters for group 1
            group2_filters: List of filters for group 2
            
        Returns:
            {
                "group1_count": int,
                "group2_count": int,
                "overlap_count": int,
                "warnings": [...]
            }
        """
        logger.info(f"Previewing cell counts for {dataset_id} (metadata-only, fast)")
        
        # Load only metadata (much faster than loading full h5ad)
        metadata_df = self._load_metadata_only(dataset_id)
        
        # Get cell masks using metadata-only method
        mask1 = self._select_cells_multi_filter_from_metadata(metadata_df, group1_filters)
        mask2 = self._select_cells_multi_filter_from_metadata(metadata_df, group2_filters)
        
        # Count cells
        group1_count = int(np.sum(mask1))
        group2_count = int(np.sum(mask2))
        
        # Check for overlap
        overlap = mask1 & mask2
        overlap_count = int(np.sum(overlap))
        
        # Generate warnings
        warnings = []
        if group1_count == 0:
            warnings.append("No cells match the filters for Group 1")
        if group2_count == 0:
            warnings.append("No cells match the filters for Group 2")
        if overlap_count > 0:
            warnings.append(f"{overlap_count} cells are in both groups (overlapping)")
        if group1_count < 10:
            warnings.append(f"Group 1 has very few cells ({group1_count}). Results may be unreliable.")
        if group2_count < 10:
            warnings.append(f"Group 2 has very few cells ({group2_count}). Results may be unreliable.")
        
        return {
            "group1_count": group1_count,
            "group2_count": group2_count,
            "overlap_count": overlap_count,
            "warnings": warnings
        }
    
    def _generate_group_description(self, filters: List[Dict[str, Any]]) -> str:
        """Generate human-readable description of filter group"""
        if not filters:
            return "All cells"
        
        descriptions = []
        for f in filters:
            column = f.get('column', '')
            values = f.get('values', [])
            if len(values) == 1:
                descriptions.append(f"{column}={values[0]}")
            elif len(values) > 1:
                descriptions.append(f"{column} in [{', '.join(values[:3])}{'...' if len(values) > 3 else ''}]")
        
        return " AND ".join(descriptions) if descriptions else "All cells"
    
    def atlas_compare(
        self,
        dataset_id: str,
        group1_filters: List[Dict[str, Any]],
        group2_filters: List[Dict[str, Any]],
        logfc_threshold: float = 0.25,
        p_value_threshold: float = 0.05,
        min_pct: float = 0.1,
        top_n: Optional[int] = None
    ) -> Dict:
        """
        Perform advanced DEG analysis between two groups defined by multi-filter criteria.
        
        Args:
            dataset_id: Dataset identifier
            group1_filters: List of {column, values} filters for group 1
            group2_filters: List of {column, values} filters for group 2
            logfc_threshold: Minimum log2 fold change
            p_value_threshold: Adjusted p-value threshold
            min_pct: Minimum fraction of cells expressing gene
            top_n: Optional limit on number of genes returned
            
        Returns:
            {
                "summary": {...},
                "group1_description": {...},
                "group2_description": {...},
                "genes": [...]
            }
        """
        logger.info(f"Starting atlas DEG comparison for {dataset_id}")
        logger.info(f"Group 1 filters: {group1_filters}")
        logger.info(f"Group 2 filters: {group2_filters}")
        
        # Load h5ad
        adata = self._load_h5ad(dataset_id)
        
        # Get cell masks
        mask1 = self._select_cells_multi_filter(adata, group1_filters)
        mask2 = self._select_cells_multi_filter(adata, group2_filters)
        
        # Count cells
        group1_count = int(np.sum(mask1))
        group2_count = int(np.sum(mask2))
        
        logger.info(f"Group 1: {group1_count} cells, Group 2: {group2_count} cells")
        
        # Validate groups
        if group1_count == 0:
            raise ValueError("No cells match the specified filters for Group 1")
        if group2_count == 0:
            raise ValueError("No cells match the specified filters for Group 2")
        
        # Check for overlap
        overlap = mask1 & mask2
        overlap_count = int(np.sum(overlap))
        if overlap_count > 0:
            raise ValueError(f"Groups have overlapping cells. {overlap_count} cells are in both groups.")
        
        # Subset data
        group1_cells = adata[mask1].copy()
        group2_cells = adata[mask2].copy()
        
        # Run DEG analysis (reuse existing logic)
        logger.info("Running DEG analysis...")
        
        # Get expression matrices
        X1 = group1_cells.X
        X2 = group2_cells.X
        
        # Handle sparse matrices
        if hasattr(X1, 'toarray'):
            X1 = X1.toarray()
        if hasattr(X2, 'toarray'):
            X2 = X2.toarray()
        
        # Apply log transformation if needed
        if self._should_log_transform(X1):
            X1 = np.log1p(X1)
        if self._should_log_transform(X2):
            X2 = np.log1p(X2)
        
        # Calculate statistics for each gene
        mean1 = np.mean(X1, axis=0)
        mean2 = np.mean(X2, axis=0)
        
        # Calculate percentage of cells expressing each gene
        pct1 = np.sum(X1 > 0, axis=0) / X1.shape[0]
        pct2 = np.sum(X2 > 0, axis=0) / X2.shape[0]
        
        # Filter genes by minimum percentage
        expressed_mask = (pct1 >= min_pct) | (pct2 >= min_pct)
        
        # Calculate log2 fold change
        log2fc = mean2 - mean1  # Already in log space
        
        # Handle NaN and infinity in log2fc
        log2fc = np.nan_to_num(log2fc, nan=0.0, posinf=0.0, neginf=0.0)
        
        # Perform t-test
        from scipy import stats
        t_stats, p_values = stats.ttest_ind(X2, X1, axis=0, equal_var=False)
        
        # Handle NaN p-values (can occur with zero variance or insufficient data)
        p_values = np.nan_to_num(p_values, nan=1.0, posinf=1.0, neginf=1.0)
        
        # Adjust p-values (Benjamini-Hochberg)
        from statsmodels.stats.multitest import multipletests
        _, p_adj, _, _ = multipletests(p_values, method='fdr_bh')
        
        # Handle NaN in adjusted p-values
        p_adj = np.nan_to_num(p_adj, nan=1.0, posinf=1.0, neginf=1.0)
        
        # Create results dataframe
        results_df = pd.DataFrame({
            'gene': adata.var_names,
            'log2_fold_change': log2fc,
            'p_value': p_values,
            'p_value_adj': p_adj,
            'mean_expr_group1': mean1,
            'mean_expr_group2': mean2,
            'pct_group1': pct1,
            'pct_group2': pct2,
            'expressed': expressed_mask
        })
        
        # Filter and mark significant genes
        results_df['significant'] = (
            (np.abs(results_df['log2_fold_change']) >= logfc_threshold) &
            (results_df['p_value_adj'] <= p_value_threshold) &
            results_df['expressed']
        )
        
        # Sort by significance
        results_df = results_df.sort_values('p_value_adj')
        
        # Convert to list of dicts, handling NaN and infinity values
        genes_list = []
        for _, row in results_df.iterrows():
            if row['expressed']:  # Only include expressed genes
                # Helper function to safely convert floats (handle NaN and inf)
                def safe_float(val):
                    if pd.isna(val) or np.isinf(val):
                        return 0.0
                    return float(val)
                
                genes_list.append({
                    'gene': row['gene'],
                    'log2_fold_change': safe_float(row['log2_fold_change']),
                    'p_value': safe_float(row['p_value']) if not pd.isna(row['p_value']) else 1.0,
                    'p_value_adj': safe_float(row['p_value_adj']) if not pd.isna(row['p_value_adj']) else 1.0,
                    'mean_expr_group1': safe_float(row['mean_expr_group1']),
                    'mean_expr_group2': safe_float(row['mean_expr_group2']),
                    'pct_group1': safe_float(row['pct_group1']),
                    'pct_group2': safe_float(row['pct_group2']),
                    'significant': bool(row['significant'])
                })
        
        # Generate descriptions
        group1_desc = self._generate_group_description(group1_filters)
        group2_desc = self._generate_group_description(group2_filters)
        comparison_desc = f"{group1_desc} vs {group2_desc}"
        
        # Build result
        result = {
            'summary': {
                'group1_cells': group1_count,
                'group2_cells': group2_count,
                'total_genes_tested': int(np.sum(expressed_mask)),
                'significant_genes': int(np.sum(results_df['significant'])),
                'comparison_description': comparison_desc
            },
            'group1_description': {
                'filters_applied': group1_filters,
                'cell_count': group1_count
            },
            'group2_description': {
                'filters_applied': group2_filters,
                'cell_count': group2_count
            },
            'genes': genes_list
        }
        
        if top_n:
            result['genes'] = result['genes'][:top_n]
        
        logger.info(f"Atlas DEG analysis completed. Found {result['summary']['significant_genes']} significant genes.")
        return result

