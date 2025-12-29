#!/usr/bin/env python3
"""
Script to compute UMAP/tSNE embeddings for h5ad files that don't have them
"""
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def compute_embeddings(h5ad_path: str, output_path: str = None):
    """
    Compute PCA, neighbors, UMAP, and tSNE embeddings
    """
    try:
        import anndata
        import scanpy as sc
    except ImportError:
        logger.error("Required packages not installed. Run: pip install anndata scanpy")
        return False
    
    logger.info(f"Loading h5ad file: {h5ad_path}")
    adata = anndata.read_h5ad(h5ad_path)
    
    logger.info(f"Dataset shape: {adata.shape}")
    
    # Check if embeddings already exist
    has_umap = 'X_umap' in adata.obsm
    has_tsne = 'X_tsne' in adata.obsm
    has_pca = 'X_pca' in adata.obsm
    
    if has_umap and has_tsne and has_pca:
        logger.info("All embeddings already exist!")
        return True
    
    logger.info("Computing embeddings...")
    
    # Compute PCA if not present
    if not has_pca:
        logger.info("  Computing PCA...")
        sc.pp.pca(adata, n_comps=50)
    else:
        logger.info("  PCA already exists")
    
    # Compute neighbors
    logger.info("  Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    
    # Compute UMAP if not present
    if not has_umap:
        logger.info("  Computing UMAP...")
        sc.tl.umap(adata)
    else:
        logger.info("  UMAP already exists")
    
    # Compute tSNE if not present
    if not has_tsne:
        logger.info("  Computing tSNE...")
        sc.tl.tsne(adata, n_pcs=30)
    else:
        logger.info("  tSNE already exists")
    
    # Save
    if output_path is None:
        output_path = h5ad_path
    
    logger.info(f"Saving to: {output_path}")
    adata.write(output_path)
    
    logger.info("✓ Embeddings computed successfully!")
    logger.info(f"  PCA: shape {adata.obsm['X_pca'].shape}")
    logger.info(f"  UMAP: shape {adata.obsm['X_umap'].shape}")
    logger.info(f"  tSNE: shape {adata.obsm['X_tsne'].shape}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Compute embeddings (PCA, UMAP, tSNE) for h5ad files"
    )
    parser.add_argument(
        "--h5ad",
        required=True,
        help="Path to the h5ad file"
    )
    parser.add_argument(
        "--output",
        help="Output path (default: overwrite input file)"
    )
    
    args = parser.parse_args()
    
    success = compute_embeddings(args.h5ad, args.output)
    
    if success:
        logger.info("\nNow run the preprocessing script to cache the embeddings:")
        logger.info(f"  python scripts/preprocess_h5ad.py --h5ad {args.output or args.h5ad}")
        return 0
    else:
        return 1


if __name__ == "__main__":
    exit(main())

