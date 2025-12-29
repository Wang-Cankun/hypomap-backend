# Scripts

This directory contains utility scripts for data import, preprocessing, and maintenance.

## Import Scripts

### `import_csv.py`
Import scRNA-seq datasets from CSV files into the SQLite database.

**Usage:**
```bash
python scripts/import_csv.py --csv "data/ssKIND - scRNAseq_data_for_ssKIND.csv"
python scripts/import_csv.py --replace --csv "data/ssKIND - scRNAseq_data_for_ssKIND.csv"
```

### `import_papers.py`
Import scRNA-seq research papers from CSV files.

**Usage:**
```bash
python scripts/import_papers.py --csv "data/ssKIND - scRNAseq_paper_for_ssKIND.csv"
```

### `import_spatial_data.py`
Import spatial transcriptomics datasets from CSV files.

**Usage:**
```bash
python scripts/import_spatial_data.py --csv "data/ssKIND - Spatial data for ssKIND.csv"
```

### `import_spatial_papers.py`
Import spatial transcriptomics research papers from CSV files.

**Usage:**
```bash
python scripts/import_spatial_papers.py --csv "data/ssKIND - Spatial paper for ssKIND.csv"
```

## Preprocessing Scripts

### `preprocess_h5ad.py`
Preprocess H5AD files and extract commonly accessed data for fast web serving.

**What it does:**
- Extracts embeddings (UMAP, tSNE, PCA) as JSON
- Caches metadata
- Computes gene statistics
- Creates SQLite database for fast gene search
- Extracts spatial coordinates, SVG data, DEG, deconvolution, CCC interactions
- Saves spatial images

**Usage:**
```bash
python scripts/preprocess_h5ad.py --h5ad h5ad/raw/dataset_id.h5ad --dataset-id dataset_id
```

**Output:**
Creates precomputed data in `h5ad/precomputed/{dataset_id}/`:
- `umap.json`, `tsne.json`, `pca.json` - Embeddings
- `metadata.json` - Cell metadata
- `data.db` - SQLite database for gene search
- `gene_expression/` - Individual gene expression files
- `spatial_coordinates.json` - Spatial coordinates (if available)
- `svg_data.json` - Spatially Variable Genes (if available)
- And more...

### `compute_embeddings.py`
Compute PCA, UMAP, and tSNE embeddings for H5AD files that don't have them.

**Usage:**
```bash
python scripts/compute_embeddings.py --h5ad h5ad/raw/dataset_id.h5ad
python scripts/compute_embeddings.py --h5ad h5ad/raw/dataset_id.h5ad --output h5ad/raw/dataset_id_with_embeddings.h5ad
```

**Note:** After computing embeddings, run `preprocess_h5ad.py` to cache them.

## Combined Import

For convenience, use the root-level `import_all_data.py` to import all data types at once:

```bash
python import_all_data.py --replace
```

This script calls all the individual import scripts in the correct order.

