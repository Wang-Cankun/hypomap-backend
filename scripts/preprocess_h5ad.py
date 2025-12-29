#!/usr/bin/env python3
"""
Script to preprocess h5ad files for fast web serving
Extracts and caches commonly accessed data (embeddings, metadata, gene stats)
"""
import argparse
import json
import logging
import os
import sqlite3
from pathlib import Path
from typing import Dict, Any

import numpy as np
import pandas as pd
from PIL import Image

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def convert_to_json_serializable(obj):
    """Recursively convert objects to JSON-serializable format"""
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient='records')
    elif isinstance(obj, pd.Series):
        return obj.tolist()
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.floating)):
        return float(obj)
    elif isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_to_json_serializable(item) for item in obj]
    elif isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    else:
        return str(obj)

def save_json(file_path: Path, data: Dict) -> None:
    """Save data as JSON"""
    file_path.parent.mkdir(parents=True, exist_ok=True)
    # Convert to JSON-serializable format
    serializable_data = convert_to_json_serializable(data)
    with open(file_path, 'w') as f:
        json.dump(serializable_data, f)
    logger.info(f"Saved to {file_path}")


def preprocess_h5ad(h5ad_path: str, dataset_id: str, output_base: str = "h5ad/precomputed"):
    """
    Preprocess h5ad file and extract commonly accessed data
    
    Args:
        h5ad_path: Path to the h5ad file
        dataset_id: Dataset identifier
        output_base: Base directory for output files
    """
    try:
        import anndata
    except ImportError:
        logger.error("anndata package not installed. Run: pip install anndata")
        return False
    
    logger.info(f"Loading h5ad file: {h5ad_path}")
    
    try:
        adata = anndata.read_h5ad(h5ad_path)
    except Exception as e:
        logger.error(f"Error loading h5ad file: {e}")
        return False
    
    logger.info(f"Dataset shape: {adata.shape} (n_obs × n_vars)")
    logger.info(f"Number of cells: {adata.n_obs}")
    logger.info(f"Number of genes: {adata.n_vars}")
    
    # Create output directory
    output_dir = Path(output_base) / dataset_id
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # 1. Extract embeddings
    logger.info("Extracting embeddings...")
    embeddings = []
    
    embedding_map = {
        'X_umap': 'umap',
        'X_tsne': 'tsne',
        'X_pca': 'pca'
    }
    
    for obsm_key, embedding_name in embedding_map.items():
        if obsm_key in adata.obsm:
            logger.info(f"  Found {embedding_name} embedding")
            embeddings.append(embedding_name)
            
            coordinates = adata.obsm[obsm_key].tolist()
            cell_ids = adata.obs_names.tolist()
            
            embedding_data = {
                "embedding_type": embedding_name,
                "coordinates": coordinates,
                "cell_ids": cell_ids
            }
            
            save_json(output_dir / f"{embedding_name}.json", embedding_data)
    
    if not embeddings:
        logger.warning("No embeddings found in the h5ad file")
    
    # 2. Extract metadata
    logger.info("Extracting metadata...")
    metadata_dict = adata.obs.to_dict(orient='list')
    cell_ids = adata.obs_names.tolist()
    
    # Convert any non-serializable types
    for key, values in metadata_dict.items():
        metadata_dict[key] = [
            None if pd.isna(v) else (str(v) if not isinstance(v, (int, float, str, bool)) else v)
            for v in values
        ]
    
    metadata = {
        "columns": list(adata.obs.columns),
        "data": metadata_dict,
        "cell_ids": cell_ids
    }
    
    save_json(output_dir / "metadata.json", metadata)
    logger.info(f"  Extracted {len(metadata['columns'])} metadata columns")
    
    # 3. Compute gene statistics
    logger.info("Computing gene statistics...")
    
    # Handle sparse matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)
    
    gene_stats = {}
    for i, gene in enumerate(adata.var_names):
        if i % 1000 == 0:
            logger.info(f"  Processed {i}/{adata.n_vars} genes")
        
        gene_expr = X[:, i]
        gene_stats[gene] = {
            "mean_expression": float(np.mean(gene_expr)),
            "max_expression": float(np.max(gene_expr)),
            "expressed_cells": int(np.sum(gene_expr > 0)),
            "percent_expressed": float(np.sum(gene_expr > 0) / len(gene_expr) * 100)
        }
    
    # save_json(output_dir / "gene_stats.json", gene_stats)
    logger.info(f"  Computed statistics for {len(gene_stats)} genes")
    
    # 4. Extract spatial coordinates
    logger.info("Checking for spatial coordinates...")
    spatial_key = None
    for key in ['spatial', 'spatial_xy', 'X_spatial']:
        if key in adata.obsm:
            spatial_key = key
            logger.info(f"  Found spatial coordinates in {key}")
            break
    
    if spatial_key:
        spatial_coords = adata.obsm[spatial_key]
        if spatial_coords.shape[1] >= 2:
            coordinates = spatial_coords[:, :2].tolist()
            cell_ids = adata.obs_names.tolist()
            
            spatial_data = {
                "coordinates": coordinates,
                "cell_ids": cell_ids,
                "spatial_key": spatial_key
            }
            
            save_json(output_dir / "spatial_coordinates.json", spatial_data)
            logger.info(f"  Extracted spatial coordinates for {len(cell_ids)} cells")
    
    # 5. Extract SVG data
    logger.info("Checking for Spatially Variable Genes (SVG)...")
    svg_cols = ['gft_score', 'svg_rank', 'cutoff_gft_score', 'pvalue', 'fdr']
    available_svg_cols = [col for col in svg_cols if col in adata.var.columns]
    
    if available_svg_cols:
        logger.info(f"  Found SVG columns: {', '.join(available_svg_cols)}")
        genes_list = []
        for gene in adata.var_names:
            gene_data = {"gene": gene}
            for col in svg_cols:
                if col in adata.var.columns:
                    value = adata.var.loc[gene, col]
                    if pd.isna(value):
                        gene_data[col] = None
                    elif isinstance(value, (int, float, np.number)):
                        gene_data[col] = float(value)
                    elif isinstance(value, bool):
                        gene_data[col] = bool(value)
                    elif isinstance(value, str):
                        gene_data[col] = str(value)
                    else:
                        # Convert other types to string
                        gene_data[col] = str(value)
                else:
                    gene_data[col] = None
            genes_list.append(gene_data)
        
        # Sort by svg_rank if available, otherwise by gft_score
        if 'svg_rank' in available_svg_cols:
            genes_list.sort(key=lambda x: x['svg_rank'] if x['svg_rank'] is not None else float('inf'))
        elif 'gft_score' in available_svg_cols:
            genes_list.sort(key=lambda x: x['gft_score'] if x['gft_score'] is not None else float('-inf'), reverse=True)
        
        svg_data = {
            "genes": genes_list,
            "total_genes": len(genes_list),
            "available_columns": available_svg_cols
        }
        
        save_json(output_dir / "svg_data.json", svg_data)
        logger.info(f"  Extracted SVG data for {len(genes_list)} genes")
    
    # 6. Extract precomputed DEG
    logger.info("Checking for precomputed DEG...")
    if "rank_genes_groups" in adata.uns:
        logger.info("  Found precomputed DEG in rank_genes_groups")
        rank_genes = adata.uns["rank_genes_groups"]
        
        # Extract groups
        groups = rank_genes['names'].dtype.names if hasattr(rank_genes['names'], 'dtype') else []
        if not groups:
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
                gene_names = rank_genes['names'][group_name] if isinstance(rank_genes['names'], np.ndarray) else rank_genes['names'].get(group_name, [])
                scores = rank_genes.get('scores', {})
                logfoldchanges = rank_genes.get('logfoldchanges', {})
                pvals = rank_genes.get('pvals', {})
                pvals_adj = rank_genes.get('pvals_adj', {})
                
                genes_list = []
                for i, gene in enumerate(gene_names):
                    if pd.isna(gene) or gene == '':
                        continue
                    
                    gene_data = {"gene": str(gene)}
                    
                    try:
                        if isinstance(scores, np.ndarray) and scores.shape[0] > i:
                            val = scores[i]
                            if not pd.isna(val) and not isinstance(val, (tuple, list)):
                                gene_data["score"] = float(val)
                            else:
                                gene_data["score"] = None
                        elif isinstance(scores, dict):
                            val = scores.get(gene, None)
                            if val is not None and not isinstance(val, (tuple, list)):
                                gene_data["score"] = float(val)
                            else:
                                gene_data["score"] = None
                        else:
                            gene_data["score"] = None
                    except (ValueError, TypeError):
                        gene_data["score"] = None
                    
                    try:
                        if isinstance(logfoldchanges, np.ndarray) and logfoldchanges.shape[0] > i:
                            val = logfoldchanges[i]
                            if not pd.isna(val) and not isinstance(val, (tuple, list)):
                                gene_data["logfoldchange"] = float(val)
                            else:
                                gene_data["logfoldchange"] = None
                        elif isinstance(logfoldchanges, dict):
                            val = logfoldchanges.get(gene, None)
                            if val is not None and not isinstance(val, (tuple, list)):
                                gene_data["logfoldchange"] = float(val)
                            else:
                                gene_data["logfoldchange"] = None
                        else:
                            gene_data["logfoldchange"] = None
                    except (ValueError, TypeError):
                        gene_data["logfoldchange"] = None
                    
                    try:
                        if isinstance(pvals, np.ndarray) and pvals.shape[0] > i:
                            val = pvals[i]
                            if not pd.isna(val) and not isinstance(val, (tuple, list)):
                                gene_data["pval"] = float(val)
                            else:
                                gene_data["pval"] = None
                        elif isinstance(pvals, dict):
                            val = pvals.get(gene, None)
                            if val is not None and not isinstance(val, (tuple, list)):
                                gene_data["pval"] = float(val)
                            else:
                                gene_data["pval"] = None
                        else:
                            gene_data["pval"] = None
                    except (ValueError, TypeError):
                        gene_data["pval"] = None
                    
                    try:
                        if isinstance(pvals_adj, np.ndarray) and pvals_adj.shape[0] > i:
                            val = pvals_adj[i]
                            if not pd.isna(val) and not isinstance(val, (tuple, list)):
                                gene_data["pval_adj"] = float(val)
                            else:
                                gene_data["pval_adj"] = None
                        elif isinstance(pvals_adj, dict):
                            val = pvals_adj.get(gene, None)
                            if val is not None and not isinstance(val, (tuple, list)):
                                gene_data["pval_adj"] = float(val)
                            else:
                                gene_data["pval_adj"] = None
                        else:
                            gene_data["pval_adj"] = None
                    except (ValueError, TypeError):
                        gene_data["pval_adj"] = None
                    
                    genes_list.append(gene_data)
                
                genes_by_group[str(group_name)] = genes_list
            except Exception as e:
                logger.warning(f"  Error extracting DEG for group {group_name}: {e}")
                genes_by_group[str(group_name)] = []
        
        deg_data = {
            "groups": [str(g) for g in groups],
            "params": {k: str(v) if not isinstance(v, (int, float, str, bool, type(None))) else v for k, v in params.items()},
            "genes": genes_by_group,
            "method": rank_genes.get('params', {}).get('method', 'unknown')
        }
        
        save_json(output_dir / "deg_precomputed.json", deg_data)
        logger.info(f"  Extracted DEG for {len(groups)} groups")
    
    # 7. Extract deconvolution predictions
    logger.info("Checking for deconvolution predictions...")
    if "tangram_ct_pred" in adata.obsm:
        logger.info("  Found Tangram deconvolution predictions")
        pred_matrix = adata.obsm["tangram_ct_pred"]
        
        if isinstance(pred_matrix, np.ndarray):
            cell_types = [f"cell_type_{i}" for i in range(pred_matrix.shape[1])]
            pred_df = pd.DataFrame(pred_matrix, columns=cell_types)
        else:
            pred_df = pred_matrix
            cell_types = list(pred_df.columns)
        
        predictions = {}
        for col in pred_df.columns:
            predictions[col] = pred_df[col].tolist()
        
        cell_ids = adata.obs_names.tolist()
        
        deconv_data = {
            "predictions": predictions,
            "cell_ids": cell_ids,
            "cell_types": cell_types
        }
        
        save_json(output_dir / "deconvolution.json", deconv_data)
        logger.info(f"  Extracted deconvolution predictions for {len(cell_types)} cell types")
    
    # 8. Extract CCC interactions
    logger.info("Checking for Cell-Cell Communication (CCC) data...")
    ccc_data = None
    if "commot_user_database" in adata.uns:
        logger.info("  Found CCC data in commot_user_database")
        ccc_data = adata.uns["commot_user_database"]
    elif "commot" in adata.uns:
        logger.info("  Found CCC data in commot")
        ccc_data = adata.uns["commot"]
    
    if ccc_data is not None:
        result = {
            "data": None,
            "interactions": [],
            "ligands": [],
            "receptors": [],
            "cell_types": []
        }
        
        if isinstance(ccc_data, dict):
            # Convert dict values that might be DataFrames
            serializable_data = {}
            for k, v in ccc_data.items():
                if isinstance(v, pd.DataFrame):
                    serializable_data[k] = v.to_dict(orient='records')
                elif isinstance(v, (list, dict, str, int, float, bool, type(None))):
                    serializable_data[k] = v
                else:
                    serializable_data[k] = str(v)
            result["data"] = serializable_data
            
            if "interactions" in ccc_data:
                val = ccc_data["interactions"]
                if isinstance(val, pd.DataFrame):
                    result["interactions"] = val.to_dict(orient='records')
                elif isinstance(val, list):
                    result["interactions"] = val
            if "ligands" in ccc_data:
                val = ccc_data["ligands"]
                if isinstance(val, pd.Series):
                    result["ligands"] = val.tolist()
                elif isinstance(val, list):
                    result["ligands"] = val
            if "receptors" in ccc_data:
                val = ccc_data["receptors"]
                if isinstance(val, pd.Series):
                    result["receptors"] = val.tolist()
                elif isinstance(val, list):
                    result["receptors"] = val
            if "cell_types" in ccc_data:
                val = ccc_data["cell_types"]
                if isinstance(val, pd.Series):
                    result["cell_types"] = val.tolist()
                elif isinstance(val, list):
                    result["cell_types"] = val
        elif isinstance(ccc_data, pd.DataFrame):
            result["data"] = ccc_data.to_dict(orient='records')
            if "ligand" in ccc_data.columns:
                result["ligands"] = ccc_data["ligand"].unique().tolist()
            if "receptor" in ccc_data.columns:
                result["receptors"] = ccc_data["receptor"].unique().tolist()
        else:
            # For other types, convert to string representation
            result["data"] = str(ccc_data)
        
        save_json(output_dir / "ccc_interactions.json", result)
        logger.info(f"  Extracted CCC interaction data")
    
    # 9. Extract spatial images
    logger.info("Checking for spatial images...")
    has_image = False
    if "spatial" in adata.uns:
        spatial_data = adata.uns["spatial"]
        if isinstance(spatial_data, dict):
            samples_info = {}
            default_sample = None
            
            for sample_key, sample_data in spatial_data.items():
                if not isinstance(sample_data, dict) or "images" not in sample_data:
                    continue
                
                has_image = True
                images = sample_data["images"]
                if isinstance(images, dict):
                    sample_info = {
                        "image_keys": [],
                        "image_shapes": {},
                        "scalefactors": {},
                        "image_paths": {}
                    }
                    
                    for img_key, img_array in images.items():
                        if isinstance(img_array, np.ndarray):
                            sample_info["image_keys"].append(img_key)
                            sample_info["image_shapes"][img_key] = list(img_array.shape)
                            
                            # Save image
                            image_dir = output_dir / "spatial_images" / sample_key
                            image_dir.mkdir(parents=True, exist_ok=True)
                            image_path = image_dir / f"{img_key}.png"
                            
                            # Normalize and save
                            if img_array.dtype != np.uint8:
                                if img_array.max() > 1.0:
                                    img_array = img_array.astype(np.uint8)
                                else:
                                    img_array = (img_array * 255).astype(np.uint8)
                            
                            if len(img_array.shape) == 2:
                                img = Image.fromarray(img_array, mode='L')
                            elif len(img_array.shape) == 3:
                                if img_array.shape[2] == 3:
                                    img = Image.fromarray(img_array, mode='RGB')
                                elif img_array.shape[2] == 4:
                                    img = Image.fromarray(img_array, mode='RGBA')
                                else:
                                    img = Image.fromarray(img_array[:, :, :3], mode='RGB')
                            else:
                                continue
                            
                            img.save(image_path, "PNG")
                            sample_info["image_paths"][img_key] = str(image_path.relative_to(output_dir))
                    
                    if "scalefactors" in sample_data:
                        scalefactors = sample_data["scalefactors"]
                        if isinstance(scalefactors, dict):
                            sample_info["scalefactors"] = {
                                k: float(v) if isinstance(v, (int, float, np.number)) else v
                                for k, v in scalefactors.items()
                            }
                    
                    if sample_info["image_keys"]:
                        samples_info[sample_key] = sample_info
                        if default_sample is None:
                            default_sample = sample_key
            
            if samples_info:
                image_info = {
                    "samples": samples_info,
                    "default_sample": default_sample
                }
                save_json(output_dir / "spatial_image_info.json", image_info)
                logger.info(f"  Extracted {len(samples_info)} spatial image samples")
    
    # 10. Detect spatial features and create spatial info
    logger.info("Detecting spatial features...")
    has_spatial = spatial_key is not None
    has_svg = len(available_svg_cols) > 0
    has_deg = "rank_genes_groups" in adata.uns
    has_deconv = "tangram_ct_pred" in adata.obsm
    has_ccc = ccc_data is not None
    
    dataset_type = None
    if has_deconv:
        dataset_type = "visium"
    elif has_ccc:
        dataset_type = "xenium"
    elif has_spatial:
        dataset_type = "unknown"
    
    spatial_info = {
        "has_spatial_coordinates": has_spatial,
        "has_svg": has_svg,
        "has_precomputed_deg": has_deg,
        "has_deconvolution": has_deconv,
        "has_ccc": has_ccc,
        "has_spatial_image": has_image,
        "dataset_type": dataset_type
    }
    
    save_json(output_dir / "spatial_info.json", spatial_info)

    # 11. Create SQLite database for fast search
    logger.info("Creating SQLite database for fast search...")
    try:
        db_path = output_dir / "data.db"
        # Remove existing db if it exists to start fresh
        if db_path.exists():
            db_path.unlink()
            
        conn = sqlite3.connect(str(db_path))
        
        # Create gene_stats table with virtual columns
        conn.execute("""
            CREATE TABLE gene_stats (
                gene TEXT PRIMARY KEY,
                stats_json JSON,
                
                -- Virtual computed columns for fast indexing
                mean_expression REAL GENERATED ALWAYS AS (json_extract(stats_json, '$.mean_expression')) VIRTUAL,
                max_expression REAL GENERATED ALWAYS AS (json_extract(stats_json, '$.max_expression')) VIRTUAL,
                expressed_cells INTEGER GENERATED ALWAYS AS (json_extract(stats_json, '$.expressed_cells')) VIRTUAL,
                percent_expressed REAL GENERATED ALWAYS AS (json_extract(stats_json, '$.percent_expressed')) VIRTUAL,
                
                -- Case-insensitive search column
                gene_lower TEXT GENERATED ALWAYS AS (LOWER(gene)) VIRTUAL
            )
        """)
        
        # Create indexes
        conn.execute("CREATE INDEX idx_gene_search ON gene_stats(gene)")
        conn.execute("CREATE INDEX idx_gene_lower ON gene_stats(gene_lower)")
        conn.execute("CREATE INDEX idx_mean_expr ON gene_stats(mean_expression DESC)")
        conn.execute("CREATE INDEX idx_expressed_cells ON gene_stats(expressed_cells DESC)")
        
        # Insert data
        # gene_stats dict keys are gene names, values are stats dicts
        data_to_insert = [(gene, json.dumps(stats)) for gene, stats in gene_stats.items()]
        
        conn.executemany(
            "INSERT INTO gene_stats (gene, stats_json) VALUES (?, ?)",
            data_to_insert
        )
        
        conn.commit()
        conn.close()
        logger.info(f"  Created SQLite database with {len(data_to_insert)} gene stats")
        
    except Exception as e:
        logger.error(f"  Error creating SQLite database: {e}")
        # Don't fail the whole process, just log error

    
    # 10. Create dataset info file
    logger.info("Creating dataset info...")
    dataset_info = {
        "dataset_id": dataset_id,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "embeddings": embeddings,
        "metadata_columns": list(adata.obs.columns),
        "spatial_features": spatial_info,
        "precomputed": True
    }
    
    save_json(output_dir / "info.json", dataset_info)
    
    logger.info("=" * 60)
    logger.info("Preprocessing completed successfully!")
    logger.info(f"Dataset: {dataset_id}")
    logger.info(f"Cells: {adata.n_obs:,}")
    logger.info(f"Genes: {adata.n_vars:,}")
    logger.info(f"Embeddings: {', '.join(embeddings) if embeddings else 'None'}")
    logger.info(f"Metadata columns: {len(metadata['columns'])}")
    logger.info(f"Spatial features: {', '.join([k for k, v in spatial_info.items() if v and k != 'dataset_type']) if any(v for k, v in spatial_info.items() if k != 'dataset_type') else 'None'}")
    if spatial_info.get("dataset_type"):
        logger.info(f"Dataset type: {spatial_info['dataset_type']}")
    logger.info("=" * 60)
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess h5ad files for fast web serving"
    )
    parser.add_argument(
        "--h5ad",
        required=True,
        help="Path to the h5ad file"
    )
    parser.add_argument(
        "--dataset-id",
        help="Dataset identifier (default: filename without extension)"
    )
    parser.add_argument(
        "--output",
        default="h5ad/precomputed",
        help="Output base directory (default: h5ad/precomputed)"
    )
    
    args = parser.parse_args()
    
    # Determine dataset ID
    if args.dataset_id:
        dataset_id = args.dataset_id
    else:
        dataset_id = Path(args.h5ad).stem
    
    # Check if file exists
    if not Path(args.h5ad).exists():
        logger.error(f"H5AD file not found: {args.h5ad}")
        return 1
    
    # Preprocess
    success = preprocess_h5ad(args.h5ad, dataset_id, args.output)
    
    if success:
        logger.info(f"\nYou can now access this dataset via the API:")
        logger.info(f"  GET /api/v1/h5ad/{dataset_id}/info")
        logger.info(f"  GET /api/v1/h5ad/{dataset_id}/embedding/umap")
        logger.info(f"  GET /api/v1/h5ad/{dataset_id}/expression/{{gene}}")
        return 0
    else:
        return 1


if __name__ == "__main__":
    exit(main())

