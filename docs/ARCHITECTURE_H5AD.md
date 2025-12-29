# H5AD Data Architecture Design

## Overview
Architecture for efficiently serving single-cell h5ad data to frontend for visualization and analysis.

## Key Design Principles

### 1. **Lazy Loading & Caching**
- Don't load entire h5ad files into memory
- Cache frequently accessed data (UMAP coordinates, metadata)
- Use Redis or in-memory cache for hot data

### 2. **Precomputed Data Strategy**
- Pre-extract and store commonly accessed data:
  - UMAP/tSNE coordinates → JSON/Parquet
  - Cell metadata → JSON/Parquet
  - Gene statistics (mean, variance) → JSON
  
### 3. **On-Demand Gene Expression**
- Load gene expression only when requested
- Support batch gene queries
- Cache recent gene expression queries

## Architecture Components

### A. Database Layer (SQLite)
Store metadata and precomputed summaries:

```sql
h5ad_files:
  - id
  - dataset_id (links to scrna_datasets)
  - file_path
  - n_cells
  - n_genes
  - has_umap
  - has_tsne
  - has_pca
  - available_embeddings (JSON)
  - cell_type_column
  - available_metadata_columns (JSON)
  - created_at

h5ad_precomputed:
  - id
  - h5ad_file_id
  - data_type (umap_coords, metadata, gene_stats)
  - file_path (path to cached JSON/Parquet)
  - created_at
```

### B. File System Structure
```
h5ad/
  ├── raw/                    # Original h5ad files
  │   └── AD093044.h5ad
  │
  └── precomputed/           # Cached data for fast access
      └── AD093044/
          ├── metadata.json          # Cell metadata
          ├── umap.json              # UMAP coordinates
          ├── tsne.json              # tSNE if available
          ├── pca.json               # PCA if available
          ├── gene_stats.json        # Gene statistics
          └── gene_expression/       # Individual gene expression
              ├── APOE.json
              ├── APP.json
              └── ...
```

### C. API Endpoints

#### 1. List Available Datasets
```
GET /api/v1/h5ad/datasets
Response: [
  {
    "id": "AD093044",
    "dataset_id": "AD093044",
    "n_cells": 50000,
    "n_genes": 20000,
    "embeddings": ["umap", "tsne"],
    "metadata_columns": ["cell_type", "cluster", "sample"]
  }
]
```

#### 2. Get Dataset Metadata
```
GET /api/v1/h5ad/{dataset_id}/metadata
Optional params: ?columns=cell_type,cluster
Response: {
  "cells": ["cell_1", "cell_2", ...],
  "metadata": {
    "cell_type": ["Neuron", "Astrocyte", ...],
    "cluster": ["1", "2", ...]
  }
}
```

#### 3. Get Embedding Coordinates
```
GET /api/v1/h5ad/{dataset_id}/embedding/{type}
where type = umap|tsne|pca
Optional params: ?limit=10000&offset=0
Response: {
  "embedding_type": "umap",
  "coordinates": [
    [x1, y1],
    [x2, y2],
    ...
  ],
  "cell_ids": ["cell_1", "cell_2", ...]
}
```

#### 4. Get Gene Expression
```
GET /api/v1/h5ad/{dataset_id}/expression/{gene_symbol}
Response: {
  "gene": "APOE",
  "expression": [0.5, 1.2, 0.0, ...],
  "cell_ids": ["cell_1", "cell_2", ...],
  "stats": {
    "mean": 0.45,
    "max": 5.2,
    "min": 0.0,
    "expressed_cells": 15234
  }
}
```

#### 5. Get Multiple Genes (Batch)
```
POST /api/v1/h5ad/{dataset_id}/expression/batch
Body: {
  "genes": ["APOE", "APP", "MAPT"]
}
Response: {
  "APOE": [...],
  "APP": [...],
  "MAPT": [...]
}
```

#### 6. Search Genes
```
GET /api/v1/h5ad/{dataset_id}/genes?search=APO
Response: [
  {"symbol": "APOE", "mean_expression": 0.45},
  {"symbol": "APOA1", "mean_expression": 0.12}
]
```

#### 7. Get Combined Data (Optimized for Plotting)
```
GET /api/v1/h5ad/{dataset_id}/plot-data
Params: 
  - embedding: umap|tsne
  - genes: APOE,APP (optional)
  - metadata: cell_type,cluster (optional)
  - subsample: 10000 (optional, for large datasets)

Response: {
  "coordinates": [[x, y], ...],
  "genes": {
    "APOE": [...],
    "APP": [...]
  },
  "metadata": {
    "cell_type": [...],
    "cluster": [...]
  },
  "cell_ids": [...]
}
```

## Implementation Strategy

### Phase 1: Precomputation Script
Create a script to preprocess h5ad files:

```python
# preprocess_h5ad.py
def preprocess_h5ad(h5ad_path, dataset_id):
    """
    Extract and cache commonly accessed data
    """
    adata = anndata.read_h5ad(h5ad_path)
    
    output_dir = f"h5ad/precomputed/{dataset_id}"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Extract UMAP/embeddings
    if 'X_umap' in adata.obsm:
        save_json(f"{output_dir}/umap.json", {
            "coordinates": adata.obsm['X_umap'].tolist(),
            "cell_ids": adata.obs_names.tolist()
        })
    
    # 2. Extract metadata
    save_json(f"{output_dir}/metadata.json", {
        "columns": adata.obs.columns.tolist(),
        "data": adata.obs.to_dict(orient='list'),
        "cell_ids": adata.obs_names.tolist()
    })
    
    # 3. Compute gene statistics
    gene_stats = compute_gene_stats(adata)
    save_json(f"{output_dir}/gene_stats.json", gene_stats)
    
    # 4. Update database
    register_h5ad_in_db(dataset_id, adata, output_dir)
```

### Phase 2: H5AD Service Layer
```python
# app/services/h5ad_service.py
class H5ADService:
    def __init__(self):
        self.cache = {}  # In-memory cache
    
    def get_embedding(self, dataset_id: str, embedding_type: str):
        # Check cache first
        cache_key = f"{dataset_id}:{embedding_type}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Load from precomputed file
        file_path = f"h5ad/precomputed/{dataset_id}/{embedding_type}.json"
        if os.path.exists(file_path):
            data = load_json(file_path)
            self.cache[cache_key] = data
            return data
        
        # Fallback: load from h5ad
        return self._extract_from_h5ad(dataset_id, embedding_type)
    
    def get_gene_expression(self, dataset_id: str, gene: str):
        # Check precomputed cache
        cache_file = f"h5ad/precomputed/{dataset_id}/gene_expression/{gene}.json"
        if os.path.exists(cache_file):
            return load_json(cache_file)
        
        # Load from h5ad on-demand
        adata = self._load_h5ad(dataset_id)
        if gene not in adata.var_names:
            raise ValueError(f"Gene {gene} not found")
        
        expression = adata[:, gene].X.toarray().flatten().tolist()
        
        # Cache for future use
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        save_json(cache_file, {
            "gene": gene,
            "expression": expression,
            "cell_ids": adata.obs_names.tolist()
        })
        
        return expression
```

### Phase 3: Optimization Strategies

#### For Large Datasets (>100k cells):
1. **Subsampling**: Return subset of points initially, load more on zoom
2. **Decimation**: Use algorithms to preserve visual structure while reducing points
3. **Spatial indexing**: Use QuadTree/KD-tree for region-based queries
4. **Parquet format**: Use Parquet instead of JSON for better compression

#### For Gene Expression:
1. **Sparse matrix handling**: Most gene expression is sparse (zeros)
2. **Send only non-zero values**: `{indices: [1,5,10], values: [1.2, 0.5, 2.1]}`
3. **Quantization**: Reduce precision for visualization (2-3 decimals)

#### Caching Strategy:
1. **LRU Cache**: Keep most recently used gene expressions in memory
2. **Redis**: For distributed/production deployment
3. **ETag/If-None-Match**: Browser caching for static data

## Technology Stack

### Required Python Packages:
```
anndata>=0.9.0
scanpy>=1.9.0
h5py>=3.8.0
numpy>=1.24.0
pandas>=2.0.0
pyarrow>=12.0.0  # For Parquet support
```

### Optional (for production):
```
redis>=4.5.0  # Caching
fastapi-cache2>=0.2.0  # FastAPI caching
```

## Performance Considerations

### Initial Load Time:
- Metadata + UMAP: ~50-100ms (10k cells)
- Gene expression: ~20-50ms per gene (cached)

### Memory Usage:
- Don't keep h5ad in memory
- Open, extract, close immediately
- Cache only JSON/serialized data

### Scalability:
- Each h5ad file is independent
- Can serve multiple datasets simultaneously
- Consider async/background jobs for precomputation

## Frontend Integration

### Initial Plot:
```javascript
// 1. Fetch UMAP coordinates
const umap = await fetch(`/api/v1/h5ad/${datasetId}/embedding/umap`)

// 2. Plot points
plotScatter(umap.coordinates)

// 3. User selects gene
const expression = await fetch(`/api/v1/h5ad/${datasetId}/expression/APOE`)

// 4. Color points by expression
updateScatterColor(expression.expression)
```

### Optimized (Single Request):
```javascript
const plotData = await fetch(
  `/api/v1/h5ad/${datasetId}/plot-data?embedding=umap&genes=APOE&metadata=cell_type`
)
// Contains everything needed for initial visualization
```

## Security Considerations

1. **File path validation**: Prevent directory traversal
2. **Rate limiting**: Limit gene expression queries
3. **Dataset access control**: Link to user permissions if needed
4. **File size limits**: Reject extremely large query responses

## Next Steps

1. Install required packages (anndata, scanpy, h5py)
2. Create database models for h5ad_files
3. Create precomputation script
4. Implement H5ADService class
5. Add API endpoints
6. Test with AD093044.h5ad
7. Add caching layer
8. Optimize for production

