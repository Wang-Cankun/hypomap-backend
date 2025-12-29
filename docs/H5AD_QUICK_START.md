# H5AD API Quick Start

## TL;DR - Get Started in 30 Seconds

```bash
# 1. Your h5ad file should be in h5ad/raw/
# 2. Preprocess it (one-time, extracts embeddings & metadata)
python scripts/preprocess_h5ad.py --h5ad h5ad/raw/AD093044.h5ad

# 3. Start server (if not already running)
pm2 start ecosystem.config.js

# 4. Test it
curl http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/embedding/umap
```

## Frontend: Show UMAP Plot (1 API Call)

```javascript
// Get everything you need in ONE request
const response = await fetch(
  'http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE&metadata=cell_type'
);
const data = await response.json();

// data contains:
// - coordinates: [[x, y], ...] for UMAP
// - genes: { APOE: [expr values] }
// - metadata: { cell_type: [...] }
// - cell_ids: [...]

// Plot with Plotly
Plotly.newPlot('myDiv', [{
  x: data.coordinates.map(c => c[0]),
  y: data.coordinates.map(c => c[1]),
  mode: 'markers',
  type: 'scattergl',
  marker: {
    color: data.genes.APOE,  // Color by gene expression
    colorscale: 'Viridis'
  }
}]);
```

## All API Endpoints

| What | Endpoint | Example |
|------|----------|---------|
| List datasets | `GET /h5ad/datasets` | List all available |
| Get UMAP | `GET /h5ad/{id}/embedding/umap` | Coordinates |
| Get tSNE | `GET /h5ad/{id}/embedding/tsne` | Coordinates |
| Gene expression | `GET /h5ad/{id}/expression/{gene}` | APOE values |
| Cell metadata | `GET /h5ad/{id}/metadata?columns=cell_type` | Metadata |
| Search genes | `GET /h5ad/{id}/genes/search?q=APO` | Find genes |
| **Combined (best)** | `GET /h5ad/{id}/plot-data?...` | Everything |

Base URL: `http://localhost:9117/sskind-backend/api/v1`

## Use Cases

### 1. Basic UMAP Scatter
```javascript
const umap = await fetch('/api/v1/h5ad/AD093044/embedding/umap').then(r => r.json());
// Plot umap.coordinates
```

### 2. Color by Gene Expression
```javascript
const data = await fetch(
  '/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE'
).then(r => r.json());
// Color by data.genes.APOE
```

### 3. Color by Cell Type
```javascript
const data = await fetch(
  '/api/v1/h5ad/AD093044/plot-data?embedding=umap&metadata=cell_type'
).then(r => r.json());
// Color by data.metadata.cell_type
```

### 4. Gene Autocomplete
```javascript
const genes = await fetch(
  '/api/v1/h5ad/AD093044/genes/search?q=APO&limit=20'
).then(r => r.json());
// genes = [{symbol: "APOE", mean_expression: 0.109, ...}, ...]
```

### 5. Multiple Genes
```javascript
const data = await fetch(
  '/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE,APP,MAPT'
).then(r => r.json());
// data.genes = { APOE: [...], APP: [...], MAPT: [...] }
```

## Performance

- **First access**: 100-500ms (loads from h5ad)
- **Cached**: 10-50ms (loads from JSON cache)
- **Dataset size**: Works great for <100K cells

## Architecture (Simple Version)

```
Frontend Request
    ↓
FastAPI Endpoint
    ↓
Check Cache (JSON files)
    ↓
If not cached: Load from h5ad → Cache → Return
If cached: Return immediately
```

## Files

```
h5ad/
├── raw/
│   └── AD093044.h5ad              ← Your original file (read-only)
└── precomputed/
    └── AD093044/
        ├── umap.json              ← Cached coordinates
        ├── metadata.json          ← Cached metadata
        ├── gene_stats.json        ← Gene info
        └── gene_expression/
            ├── APOE.json          ← Cached on first request
            └── APP.json
```

## Common Tasks

### Add new dataset
```bash
cp new_data.h5ad h5ad/raw/
python preprocess_h5ad.py --h5ad h5ad/raw/new_data.h5ad
pm2 restart sskind-backend
```

### Clear cache (force reprocessing)
```bash
rm -rf h5ad/precomputed/AD093044/
```

### Check what's available
```bash
curl http://localhost:9117/sskind-backend/api/v1/h5ad/datasets | python3 -m json.tool
```

## Full Documentation

- **Architecture**: See `ARCHITECTURE_H5AD.md`
- **API Guide**: See `H5AD_API_GUIDE.md`
- **Summary**: See `H5AD_ARCHITECTURE_SUMMARY.md`
- **Interactive Docs**: http://localhost:9117/sskind-backend/docs

## Need Help?

1. Check PM2 logs: `pm2 logs sskind-backend`
2. Test endpoint: `curl http://localhost:9117/sskind-backend/api/v1/h5ad/datasets`
3. Verify file exists: `ls -lh h5ad/raw/`
4. Check preprocessing: `ls -lh h5ad/precomputed/AD093044/`

