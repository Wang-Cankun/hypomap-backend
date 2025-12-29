# ssKIND Backend API Summary

## Quick Overview

The ssKIND backend provides a comprehensive REST API for accessing and analyzing single-cell and spatial transcriptomics data. All endpoints return JSON and support CORS for frontend integration.

**Base URL**: `http://localhost:9117/sskind-backend/api/v1`

**Interactive Docs**: `http://localhost:9117/sskind-backend/docs`

---

## API Endpoints Summary

### 📊 Datasets API

#### scRNA-seq Datasets
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/datasets/` | GET | List all datasets |
| `/datasets/{id}` | GET | Get specific dataset |
| `/datasets/stats/` | GET | Get statistics |
| `/datasets/search/disease/{disease}` | GET | Search by disease |
| `/datasets/search/species/{species}` | GET | Search by species |
| `/datasets/search/brain-region/{region}` | GET | Search by brain region |

#### scRNA-seq Papers
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/papers/` | GET | List all papers |
| `/papers/{id}` | GET | Get specific paper |
| `/papers/stats/` | GET | Get statistics |
| `/papers/search/disease/{disease}` | GET | Search by disease |
| `/papers/search/species/{species}` | GET | Search by species |
| `/papers/search/brain-region/{region}` | GET | Search by brain region |

#### Spatial Datasets
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/spatial-datasets/` | GET | List all spatial datasets |
| `/spatial-datasets/{id}` | GET | Get specific dataset |
| `/spatial-datasets/stats/` | GET | Get statistics |
| `/spatial-datasets/search/disease/{disease}` | GET | Search by disease |
| `/spatial-datasets/search/species/{species}` | GET | Search by species |
| `/spatial-datasets/search/brain-region/{region}` | GET | Search by brain region |
| `/spatial-datasets/search/methodology/{method}` | GET | Search by methodology |
| `/spatial-datasets/search/study/{study_id}` | GET | Get all datasets from study |

#### Spatial Papers
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/spatial-papers/` | GET | List all spatial papers |
| `/spatial-papers/{id}` | GET | Get specific paper |
| `/spatial-papers/stats/` | GET | Get statistics |
| `/spatial-papers/search/disease/{disease}` | GET | Search by disease |
| `/spatial-papers/search/species/{species}` | GET | Search by species |
| `/spatial-papers/search/brain-region/{region}` | GET | Search by brain region |
| `/spatial-papers/search/methodology/{method}` | GET | Search by methodology |

---

### 🧬 H5AD Data API

#### Dataset Information
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/h5ad/{id}/info` | GET | Get dataset metadata |
| `/h5ad/{id}/genes/search` | GET | Search genes (autocomplete) |

#### Embeddings
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/h5ad/{id}/embedding/{type}` | GET | Get UMAP/tSNE/PCA coordinates |

**Parameters**: `type` = `umap`, `tsne`, or `pca`

#### Gene Expression
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/h5ad/{id}/expression/{gene}` | GET | Get gene expression (normalized) |
| `/h5ad/{id}/module-score` | POST | Calculate module score for gene list |

**Module Score Body**:
```json
{
  "gene_list": ["APOE", "APP", "MAPT"],
  "method": "mean",
  "use_raw": false
}
```

#### Metadata
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/h5ad/{id}/metadata` | GET | Get cell metadata |

**Query Parameters**: `columns` (comma-separated list)

#### Combined Data
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/h5ad/{id}/plot-data` | GET | Get embedding + expression + metadata in one call |

**Query Parameters**:
- `embedding`: `umap`, `tsne`, or `pca`
- `genes`: comma-separated gene list
- `metadata`: comma-separated metadata columns

---

### 📈 DEG Analysis API

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/deg/between-datasets` | POST | Compare two datasets |
| `/deg/within-dataset` | POST | Compare groups within dataset |
| `/deg/cell-types/{dataset_id}` | GET | Get available cell types |
| `/deg/metadata-columns/{dataset_id}` | GET | Get metadata columns |

#### Between Datasets Body
```json
{
  "dataset_id1": "AD093044",
  "dataset_id2": "AD093044",
  "cell_types_group1": ["Astrocyte", "CGE interneuron"],
  "cell_types_group2": ["Microglia", "Oligodendrocyte"],
  "min_pct": 0.1,
  "logfc_threshold": 0.25,
  "p_value_threshold": 0.05,
  "top_n": 50
}
```

#### Within Dataset Body
```json
{
  "dataset_id": "AD093044",
  "group1_filters": {"cell_type": "Microglia", "disease": "AD"},
  "group2_filters": {"cell_type": "Microglia", "disease": "Control"},
  "min_pct": 0.1,
  "logfc_threshold": 0.25,
  "p_value_threshold": 0.05
}
```

---

## Data Statistics (Current)

### Datasets Available
- **scRNA-seq Datasets**: 3,620 datasets
- **Spatial Datasets**: 1,380 datasets
- **scRNA-seq Papers**: 232 papers
- **Spatial Papers**: 170 papers

### H5AD Files Available
- **AD093044**: 7,952 cells (UMAP, tSNE, PCA)
- **AD093045**: 7,182 cells (UMAP, tSNE, PCA)
- **AD093046**: 6,536 cells (UMAP, tSNE, PCA)

**Total cells with h5ad**: 21,670 cells

---

## Response Formats

### Standard List Response
```json
[
  {
    "id": 1,
    "dataset_id": "AD093044",
    "disease": "AD",
    "species": "Human",
    "n_cells": 7952,
    "...": "..."
  }
]
```

### Embedding Response
```json
{
  "dataset_id": "AD093044",
  "embedding_type": "umap",
  "coordinates": [[1.2, 3.4], [5.6, 7.8], ...],
  "cell_ids": ["cell_1", "cell_2", ...]
}
```

### Gene Expression Response
```json
{
  "gene": "APOE",
  "expression": [0.0, 3.2, 1.5, ...],
  "cell_ids": ["cell_1", "cell_2", ...],
  "data_type": "normalized",
  "stats": {
    "mean": 2.03,
    "median": 1.20,
    "std": 2.37,
    "min": 0.0,
    "max": 20.2,
    "expressed_cells": 7436,
    "mean_expressed": 2.18
  }
}
```

### Module Score Response
```json
{
  "module_score": [0.6, 5.6, 1.4, ...],
  "cell_ids": ["cell_1", "cell_2", ...],
  "genes_used": ["APOE", "APP", "MAPT"],
  "genes_not_found": [],
  "method": "mean",
  "n_genes_used": 3,
  "data_type": "normalized",
  "stats": {
    "mean": 2.03,
    "median": 1.20,
    "std": 2.37,
    "min": 0.0,
    "max": 20.2
  }
}
```

### DEG Analysis Response
```json
{
  "comparison": {
    "group1": {"dataset_id": "AD093044", "n_cells": 7952},
    "group2": {"dataset_id": "AD093045", "n_cells": 7182}
  },
  "parameters": {
    "min_pct": 0.1,
    "logfc_threshold": 0.25,
    "p_value_threshold": 0.05
  },
  "summary": {
    "total_genes_tested": 15234,
    "significant_genes": 823,
    "upregulated_in_group1": 412,
    "upregulated_in_group2": 411
  },
  "genes": [
    {
      "gene": "APOE",
      "log2_fold_change": 2.34,
      "mean_expr_group1": 5.6,
      "mean_expr_group2": 2.1,
      "pct_group1": 0.95,
      "pct_group2": 0.82,
      "p_value": 1.2e-45,
      "p_value_adj": 3.4e-42,
      "significant": true
    }
  ]
}
```

---

## Performance

### Response Times (with caching)
- **Embeddings**: <100ms (cached), 5-10s (first access)
- **Gene Expression**: <100ms (cached), 2-5s (first access)
- **Module Score**: 1-3s (computed on-demand)
- **DEG Analysis**: 10-60s (first), <100ms (cached)
- **Dataset List**: <50ms
- **Search**: <100ms

### Caching Strategy
- **Embeddings**: Precomputed to JSON
- **Gene Expression**: Cached after first access
- **Module Score**: Not cached (fast computation)
- **DEG Results**: Cached by parameters hash

---

## Error Responses

### 404 Not Found
```json
{
  "detail": "Dataset AD999999 not found"
}
```

### 400 Bad Request
```json
{
  "detail": "Gene 'INVALID' not found in dataset"
}
```

### 500 Internal Server Error
```json
{
  "detail": "Internal server error: [error message]"
}
```

---

## Common Use Cases

### 1. Browse Datasets
```bash
# Get all Alzheimer's datasets
curl "http://localhost:9117/sskind-backend/api/v1/datasets/search/disease/AD"
```

### 2. Visualize UMAP
```bash
# Get UMAP coordinates
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/embedding/umap"
```

### 3. Plot Gene Expression on UMAP
```bash
# Get UMAP + gene expression in one call
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE"
```

### 4. Calculate Gene Signature
```bash
# Module score for AD risk genes
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/module-score" \
  -H "Content-Type: application/json" \
  -d '{"gene_list": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"]}'
```

### 5. Find DEGs
```bash
# Compare two datasets
curl -X POST "http://localhost:9117/sskind-backend/api/v1/deg/between-datasets" \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id1": "AD093044",
    "dataset_id2": "AD093045",
    "cell_type": "Microglia",
    "top_n": 50
  }'
```

---

## Frontend Integration

### React Example

```javascript
import { useState, useEffect } from 'react';

function DatasetBrowser() {
  const [datasets, setDatasets] = useState([]);
  
  useEffect(() => {
    fetch('http://localhost:9117/sskind-backend/api/v1/datasets/')
      .then(res => res.json())
      .then(data => setDatasets(data));
  }, []);
  
  return (
    <div>
      {datasets.map(d => (
        <div key={d.dataset_id}>
          {d.dataset_id}: {d.disease}, {d.n_cells} cells
        </div>
      ))}
    </div>
  );
}
```

### Vue Example

```javascript
<template>
  <div>
    <div v-for="dataset in datasets" :key="dataset.dataset_id">
      {{ dataset.dataset_id }}: {{ dataset.disease }}
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return { datasets: [] };
  },
  async mounted() {
    const res = await fetch(
      'http://localhost:9117/sskind-backend/api/v1/datasets/'
    );
    this.datasets = await res.json();
  }
}
</script>
```

### Plotly Visualization

```javascript
// Get UMAP + gene expression
const response = await fetch(
  'http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE&metadata=cell_type'
);
const data = await response.json();

// Plot with Plotly
Plotly.newPlot('plot', [{
  x: data.embedding.coordinates.map(c => c[0]),
  y: data.embedding.coordinates.map(c => c[1]),
  mode: 'markers',
  type: 'scattergl',
  marker: {
    color: data.gene_expression.APOE.expression,
    colorscale: 'Viridis',
    showscale: true
  }
}], {
  title: 'APOE Expression on UMAP',
  xaxis: { title: 'UMAP 1' },
  yaxis: { title: 'UMAP 2' }
});
```

---

## Best Practices

### 1. Use Combined Endpoints
✅ Good: `/h5ad/{id}/plot-data?embedding=umap&genes=APOE,APP`

❌ Bad: 3 separate requests for embedding, gene1, gene2

### 2. Implement Caching on Frontend
```javascript
const cache = new Map();

async function fetchWithCache(url) {
  if (cache.has(url)) return cache.get(url);
  const response = await fetch(url);
  const data = await response.json();
  cache.set(url, data);
  return data;
}
```

### 3. Use Pagination
✅ Good: `/datasets/?limit=100&skip=0`

❌ Bad: `/datasets/` (loads all)

### 4. Handle Errors
```javascript
async function safeRequest(url) {
  try {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }
    return await response.json();
  } catch (error) {
    console.error('Request failed:', error);
    return null;
  }
}
```

---

## Rate Limits & Quotas

- **No rate limits** currently enforced
- **Recommended**: Max 100 requests/minute
- **Large datasets**: Use pagination and caching

---

## Support & Documentation

- **Full Documentation**: [`docs/`](./README.md)
- **Interactive API Docs**: http://localhost:9117/sskind-backend/docs
- **Quick Start**: [H5AD Quick Start](./H5AD_QUICK_START.md)

---

**Last Updated**: November 2024

