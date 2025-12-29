# H5AD Visualization API Guide

## Overview

The ssKIND backend now supports efficient visualization of single-cell h5ad data through a RESTful API with lazy loading and caching strategies.

## Architecture Summary

### Design Principles
1. **Read-Only Access**: h5ad files are never modified, only read
2. **Precomputation**: Frequently accessed data (UMAP, metadata, gene stats) are cached as JSON
3. **Lazy Loading**: Gene expression loaded on-demand and cached after first access
4. **In-Memory Cache**: Hot data kept in memory for fast repeated access

### File Structure
```
h5ad/
├── raw/                           # Original h5ad files (read-only)
│   └── AD093044.h5ad
└── precomputed/                   # Cached JSON data
    └── AD093044/
        ├── info.json             # Dataset metadata
        ├── umap.json             # UMAP coordinates (cached)
        ├── tsne.json             # tSNE coordinates (cached)
        ├── pca.json              # PCA coordinates (cached)
        ├── metadata.json         # Cell metadata (cached)
        ├── gene_stats.json       # Gene statistics (cached)
        └── gene_expression/      # On-demand gene expression cache
            ├── APOE.json
            ├── APP.json
            └── ...
```

## API Endpoints

### 1. List Available Datasets
```bash
GET /api/v1/h5ad/datasets
```

**Response:**
```json
[
  {
    "dataset_id": "AD093044",
    "n_cells": 7952,
    "n_genes": 33538,
    "embeddings": ["umap", "tsne", "pca"],
    "metadata_columns": ["cell_type", "cluster_name", ...],
    "gene_names": ["APOE", "APP", ...]
  }
]
```

### 2. Get Dataset Info
```bash
GET /api/v1/h5ad/{dataset_id}/info
```

**Example:**
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/info"
```

### 3. Get Embedding Coordinates
```bash
GET /api/v1/h5ad/{dataset_id}/embedding/{type}
```

**Types:** `umap`, `tsne`, `pca`

**Example:**
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/embedding/umap"
```

**Response:**
```json
{
  "embedding_type": "umap",
  "coordinates": [
    [11.93, -0.69],
    [21.42, 9.13],
    ...
  ],
  "cell_ids": ["cell_1", "cell_2", ...]
}
```

### 4. Get Gene Expression
```bash
GET /api/v1/h5ad/{dataset_id}/expression/{gene}
```

**Example:**
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/expression/APOE"
```

**Response:**
```json
{
  "gene": "APOE",
  "expression": [0.0, 0.0, 1.5, 2.3, ...],
  "cell_ids": ["cell_1", "cell_2", ...],
  "stats": {
    "mean": 0.109,
    "max": 10.0,
    "min": 0.0,
    "expressed_cells": 568,
    "mean_expressed": 1.53
  }
}
```

### 5. Get Cell Metadata
```bash
GET /api/v1/h5ad/{dataset_id}/metadata?columns=cell_type,cluster_name
```

**Example:**
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/metadata?columns=cell_type"
```

**Response:**
```json
{
  "columns": ["cell_type"],
  "data": {
    "cell_type": ["Neuron", "Astrocyte", "Microglia", ...]
  },
  "cell_ids": ["cell_1", "cell_2", ...]
}
```

### 6. Get Combined Plot Data (Optimized)
```bash
GET /api/v1/h5ad/{dataset_id}/plot-data?embedding=umap&genes=APOE,APP&metadata=cell_type
```

**Example:**
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE,APP&metadata=cell_type"
```

**Response:**
```json
{
  "coordinates": [[11.93, -0.69], ...],
  "genes": {
    "APOE": [0.0, 0.0, 1.5, ...],
    "APP": [1.2, 0.5, 2.1, ...]
  },
  "metadata": {
    "cell_type": ["Neuron", "Astrocyte", ...]
  },
  "cell_ids": ["cell_1", "cell_2", ...]
}
```

This endpoint returns everything needed for visualization in a single request!

### 7. Search Genes
```bash
GET /api/v1/h5ad/{dataset_id}/genes/search?q=APO&limit=10
```

**Example:**
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/genes/search?q=APO&limit=10"
```

**Response:**
```json
[
  {
    "symbol": "APOE",
    "mean_expression": 0.109,
    "expressed_cells": 568
  },
  {
    "symbol": "APOD",
    "mean_expression": 0.299,
    "expressed_cells": 1674
  }
]
```

## Frontend Integration Examples

### Example 1: Basic UMAP Visualization
```javascript
// Step 1: Get UMAP coordinates
const response = await fetch(
  'http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/embedding/umap'
);
const data = await response.json();

// Step 2: Plot using your favorite library (Plotly, D3, etc.)
const x = data.coordinates.map(c => c[0]);
const y = data.coordinates.map(c => c[1]);

Plotly.newPlot('myDiv', [{
  x: x,
  y: y,
  mode: 'markers',
  type: 'scatter'
}]);
```

### Example 2: Gene Expression Overlay
```javascript
// Get UMAP + gene expression in one call
const response = await fetch(
  'http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE'
);
const data = await response.json();

// Plot with color by gene expression
const x = data.coordinates.map(c => c[0]);
const y = data.coordinates.map(c => c[1]);
const color = data.genes.APOE;

Plotly.newPlot('myDiv', [{
  x: x,
  y: y,
  mode: 'markers',
  type: 'scatter',
  marker: {
    color: color,
    colorscale: 'Viridis',
    showscale: true
  }
}]);
```

### Example 3: Interactive Gene Search
```javascript
// Gene search with autocomplete
async function searchGenes(query) {
  const response = await fetch(
    `http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/genes/search?q=${query}&limit=20`
  );
  return await response.json();
}

// Usage in autocomplete component
const suggestions = await searchGenes('APO');
// Returns: [{symbol: "APOE", ...}, {symbol: "APOD", ...}, ...]
```

### Example 4: Cell Type Coloring
```javascript
// Get UMAP + cell type metadata
const response = await fetch(
  'http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&metadata=cell_type'
);
const data = await response.json();

// Create color map for cell types
const cellTypes = [...new Set(data.metadata.cell_type)];
const colorMap = {};
cellTypes.forEach((type, i) => {
  colorMap[type] = i;
});

// Plot with color by cell type
const x = data.coordinates.map(c => c[0]);
const y = data.coordinates.map(c => c[1]);
const color = data.metadata.cell_type.map(t => colorMap[t]);

Plotly.newPlot('myDiv', [{
  x: x,
  y: y,
  mode: 'markers',
  type: 'scatter',
  marker: {
    color: color,
    colorscale: 'Portland'
  }
}]);
```

## Performance Characteristics

### First Access (Cold Cache)
- **Dataset info**: ~50-100ms
- **UMAP coordinates**: ~100-200ms (7,952 cells)
- **Gene expression**: ~200-500ms per gene

### Subsequent Access (Warm Cache)
- **Cached embeddings**: ~10-20ms
- **Cached gene expression**: ~10-20ms
- **In-memory cache**: <5ms

### Memory Usage
- **Precomputed files**: ~10-50MB per dataset (depends on size)
- **In-memory cache**: Automatically managed, keeps hot data

## Adding New H5AD Files

### Step 1: Place h5ad file
```bash
cp your_dataset.h5ad h5ad/raw/
```

### Step 2: Preprocess (cache embeddings and metadata)
```bash
python scripts/preprocess_h5ad.py --h5ad h5ad/raw/your_dataset.h5ad
```

This will:
- Extract UMAP/tSNE/PCA coordinates → JSON
- Extract cell metadata → JSON
- Compute gene statistics → JSON
- Create dataset info → JSON

### Step 3: Use the API
```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/datasets"
# Your dataset will appear in the list
```

## Best Practices

### For Frontend Development

1. **Use the combined endpoint** (`/plot-data`) for initial visualization to minimize requests
2. **Cache on the frontend** when possible (e.g., UMAP coordinates don't change)
3. **Lazy load gene expression** - only fetch when user selects a gene
4. **Use gene search** for autocomplete instead of fetching all gene names
5. **Consider subsampling** for very large datasets (>100k cells) for initial rendering

### For Backend Performance

1. **Precompute datasets** before deploying to production
2. **Monitor cache size** - gene expression cache can grow large
3. **Consider Redis** for distributed deployment (optional)
4. **Use ETag/Cache-Control** headers for browser caching (future enhancement)

## Troubleshooting

### Dataset not appearing in list
- Check file is in `h5ad/raw/` directory
- Run preprocessing script
- Check PM2 logs: `pm2 logs sskind-backend`

### Slow gene expression queries
- First access is always slower (loads from h5ad)
- Subsequent accesses use cache (much faster)
- Consider precomputing popular genes

### Out of memory
- Reduce cache size in service configuration
- Precompute more data to disk instead of memory
- Consider pagination for very large datasets

## API Documentation

Full interactive API documentation available at:
```
http://localhost:9117/sskind-backend/docs
```

Look for the "H5AD Data" section for interactive testing.

