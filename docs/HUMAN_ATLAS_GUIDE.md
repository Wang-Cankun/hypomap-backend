# Human Subset Atlas Integration Guide

## Overview

`human_subset.h5ad` is a large demo atlas (200k cells × 5k genes) that showcases how to interact with a multi-cell-type dataset through the `ssKIND-backend` APIs. It includes UMAP + tSNE embeddings and the `cell_type` annotation column for coloring plots.

This guide explains how to:

- Access atlas metadata and embeddings
- Fetch per-gene expression
- Generate heatmap / dotplot summaries with optional scaling + clustering
- Visualize tSNE colored by `cell_type`

## Dataset Prep

The file lives at `h5ad/raw/human_subset.h5ad`. Run:

```bash
/Users/wang.13246/.local/share/mamba/envs/sskind/bin/python scripts/preprocess_h5ad.py \
  --h5ad h5ad/raw/human_subset.h5ad \
  --dataset-id human_subset
```

This creates `h5ad/precomputed/human_subset/` with embeddings, metadata, and gene stats.

## Key Facts

- Dataset ID: `human_subset`
- Cells: 200,000
- Genes: 5,000
- Embeddings: `umap`, `tsne`
- Metadata column to use for coloring: `cell_type`

## API Calls

### Info

```bash
curl -s http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/info | jq
```

The info endpoint now includes a `categorical_metadata` field that lists all categorical metadata columns along with their unique values. This is useful for the frontend to provide selection options for which metadata to display:

```bash
# Get available categorical metadata columns
curl -s http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/info | jq '.categorical_metadata | keys'

# Get unique values for a specific categorical column
curl -s http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/info | jq '.categorical_metadata.cell_type'
```

### Embeddings

```bash
# tSNE coordinates (200k points)
curl -s http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/embedding/tsne \
  | jq '{embedding_type, coords: ( .coordinates | length ), cells: ( .cell_ids | length )}'
```

Use the `cell_ids` order to align with metadata or heatmap output when plotting.

### Metadata

You can query any categorical metadata column, not just `cell_type`:

```bash
# Query cell_type
curl -s "http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/metadata?columns=cell_type" \
  | jq '.data.cell_type[:10]'

# Query multiple categorical columns
curl -s "http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/metadata?columns=cell_type,supercluster_name,cluster_name" \
  | jq '.data | keys'

# Query any available categorical column (see info endpoint for full list)
curl -s "http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/metadata?columns=_scvi_labels" \
  | jq '.data._scvi_labels[:10]'
```

### Per-gene expression

```bash
curl -s http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/expression/APOE | jq '.stats'
```

### Heatmap / Dotplot

#### Heatmap with scaling + clustering

```bash
curl -s -X POST http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/heatmap \
  -H "Content-Type: application/json" \
  -d '{
        "gene_list": ["APOE", "TREM2", "MBP"],
        "plot_type": "heatmap",
        "scale_expression": true,
        "cluster_rows": true,
        "cluster_columns": true
      }' | jq '{genes, cell_types: .cell_types[:5], stats}'
```

#### Dotplot (expression + percent expressing)

```bash
curl -s -X POST http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/heatmap \
  -H "Content-Type: application/json" \
  -d '{
        "gene_list": ["APOE", "TREM2", "MBP"],
        "plot_type": "dotplot"
      }' | jq '.data[0][0]'
```

### Combined Plot Data (embedding + metadata + expression)

You can query with any categorical metadata column(s):

```bash
# Query with cell_type
curl -s "http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/plot-data?embedding=tsne&genes=APOE&metadata=cell_type" \
  | jq '{embedding_type: .embedding.embedding_type, points: (.embedding.coordinates | length)}'

# Query with multiple metadata columns
curl -s "http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/plot-data?embedding=umap&metadata=cell_type,supercluster_name,cluster_name" \
  | jq '.metadata | keys'
```

## Frontend Integration Example

### Step 1: Get available categorical metadata options

```javascript
// Fetch dataset info to get available categorical metadata
const info = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/info"
).then((r) => r.json());

// Display available options to user
const metadataOptions = Object.keys(info.categorical_metadata);
// ["Pubmed_id", "_scvi_labels", "cell_type", "cluster_name", "supercluster_name", ...]

// Get unique values for a selected metadata column (e.g., for dropdown)
const cellTypeOptions = info.categorical_metadata.cell_type;
// ["Amygdala excitatory", "Astrocyte", "Bergmann glia", ...]
```

### Step 2: Query plot data with selected metadata

```javascript
// User selects metadata columns to display
const selectedMetadata = ["cell_type", "supercluster_name", "cluster_name"];

// Fetch plot data
const plotData = await fetch(
  `http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/plot-data?` +
    `embedding=umap&metadata=${selectedMetadata.join(",")}`
).then((r) => r.json());

// plotData.coordinates: [[x1, y1], [x2, y2], ...]
// plotData.metadata.cell_type: ["Neuron", "Astrocyte", ...]
// plotData.metadata.supercluster_name: ["Excitatory", "Inhibitory", ...]
```

## Plotting tSNE in Scanpy

```python
import scanpy as sc
adata = sc.read_h5ad("h5ad/raw/human_subset.h5ad")
sc.pl.tsne(adata, color="cell_type", frameon=False, show=False)
```

Use the API's `embedding/tsne` + `metadata` endpoints if you need the same visualization client-side.

## Notes

- The API preserves gene + cell type order as they appear in the dataset unless clustering is requested.
- The `categorical_metadata` field in the info endpoint provides all available categorical columns with their unique values, making it easy for the frontend to provide selection options.
- While `cell_type` is commonly used for coloring, you can use any categorical metadata column (e.g., `supercluster_name`, `cluster_name`, `_scvi_labels`, etc.).
- Categorical metadata detection includes:
  - String/object columns with < 1000 unique values
  - Boolean columns
  - Numeric columns with < 100 unique values (e.g., cluster labels)
- Remember to document new datasets (add them to README/API summaries as needed).
