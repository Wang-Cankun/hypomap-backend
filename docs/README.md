# ssKIND Backend API Documentation

Welcome to the ssKIND backend API documentation! This comprehensive guide covers all features for serving and analyzing single-cell and spatial transcriptomics data.

## 📚 Documentation Index

### Getting Started

1. **[API Summary](./API_SUMMARY.md)** 📋
   - Complete endpoint reference
   - Response formats
   - Quick examples
   - Performance metrics

2. **[Quick Start Guide](./H5AD_QUICK_START.md)** ⚡
   - Get started in 30 seconds
   - Basic examples
   - Essential endpoints

### Core Features

3. **[Datasets API Guide](./DATASETS_API_GUIDE.md)** 📊
   - Browse and search datasets
   - Metadata queries
   - Paper information
   - scRNA-seq and spatial datasets

4. **[H5AD Data API Guide](./H5AD_API_GUIDE.md)** 🧬
   - UMAP/tSNE visualization
   - Gene expression queries
   - Metadata retrieval
   - Combined plot data
   - Gene search

### Advanced Analysis

5. **[DEG Analysis Guide](./DEG_ANALYSIS_GUIDE.md)** 📈
   - Differential expression analysis
   - Between-dataset comparisons
   - Within-dataset comparisons
   - Cell type markers
   - Statistical methods

6. **[Module Score Guide](./MODULE_SCORE_GUIDE.md)** 🎯
   - Gene signature scoring
   - Pathway activity
   - Multi-gene analysis
   - Custom scoring methods

### Performance & Optimization

12. **[SQLite Performance Report](./PERFORMANCE_REPORT.md)** ⚡
    - SQLite gene search performance benchmarks
    - Memory efficiency analysis
    - Query optimization details
    - Performance improvements

13. **[SQLite Migration Summary](./SQLITE_MIGRATION_SUMMARY.md)** 🚀
    - Migration implementation details
    - Performance verification results
    - Database schema and indexes
    - API compatibility

### Architecture & Implementation

7. **[H5AD Architecture Summary](./H5AD_ARCHITECTURE_SUMMARY.md)** 🏗️
   - System design
   - Data flow
   - Caching strategy
   - Performance metrics
   - Scalability

8. **[H5AD Architecture Details](./ARCHITECTURE_H5AD.md)** 🔧
   - Detailed technical design
   - Implementation strategy
   - Database schema
   - File organization

## 🚀 Quick Links

### Common Tasks

| Task | Endpoint | Documentation |
|------|----------|---------------|
| List all datasets | `GET /datasets/` | [Datasets API](./DATASETS_API_GUIDE.md#list-all-datasets) |
| Get UMAP coordinates | `GET /h5ad/{id}/embedding/umap` | [H5AD API](./H5AD_API_GUIDE.md#get-embedding-coordinates) |
| Get gene expression | `GET /h5ad/{id}/expression/{gene}` | [H5AD API](./H5AD_API_GUIDE.md#get-gene-expression) |
| Calculate DEG | `POST /deg/between-datasets` | [DEG Guide](./DEG_ANALYSIS_GUIDE.md#deg-between-two-datasets) |
| Module score | `POST /h5ad/{id}/module-score` | [Module Score](./MODULE_SCORE_GUIDE.md#endpoint) |
| Search genes | `GET /h5ad/{id}/genes/search` | [H5AD API](./H5AD_API_GUIDE.md#search-genes) |

### API Base URL

```
http://localhost:9117/sskind-backend/api/v1
```

### Interactive Documentation

```
http://localhost:9117/sskind-backend/docs
```

## 📖 Guide Selection Helper

### "I want to..."

- **Browse available datasets** → [Datasets API Guide](./DATASETS_API_GUIDE.md)
- **Visualize cells on UMAP** → [H5AD Quick Start](./H5AD_QUICK_START.md)
- **Find differentially expressed genes** → [DEG Analysis Guide](./DEG_ANALYSIS_GUIDE.md)
- **Calculate gene signatures** → [Module Score Guide](./MODULE_SCORE_GUIDE.md)
- **Understand the system architecture** → [Architecture Summary](./H5AD_ARCHITECTURE_SUMMARY.md)
- **Get started quickly** → [Quick Start Guide](./H5AD_QUICK_START.md)

## 🎯 Use Case Examples

### Research Scenarios

#### 1. Explore Dataset Collection
```javascript
// Find Alzheimer's disease datasets
const datasets = await fetch('/api/v1/datasets/search/disease/AD');
// → See: Datasets API Guide
```

#### 2. Visualize Single-Cell Data
```javascript
// Get UMAP + gene expression in one call
const data = await fetch(
  '/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE&metadata=cell_type'
);
// → See: H5AD API Guide
```

#### 3. Compare Conditions
```javascript
// Find genes differentially expressed between datasets
const deg = await fetch('/api/v1/deg/between-datasets', {
  method: 'POST',
  body: JSON.stringify({
    dataset_id1: 'AD093044',
    dataset_id2: 'AD093045',
    cell_type: 'Microglia'
  })
});
// → See: DEG Analysis Guide
```

#### 4. Score Gene Signatures
```javascript
// Calculate Alzheimer's risk score
const score = await fetch('/api/v1/h5ad/AD093044/module-score', {
  method: 'POST',
  body: JSON.stringify({
    gene_list: ['APOE', 'APP', 'MAPT', 'PSEN1', 'PSEN2']
  })
});
// → See: Module Score Guide
```

## 🛠️ Technology Stack

- **Framework**: FastAPI (Python)
- **Database**: SQLite with SQLAlchemy ORM
- **Data Format**: h5ad (AnnData)
- **Analysis**: scanpy, scipy, statsmodels
- **Caching**: JSON precomputation + in-memory cache

## 📊 Current Capabilities

### Datasets
- ✅ 3 scRNA-seq datasets (21,670 cells total)
- ✅ 1,380 spatial datasets
- ✅ 58 spatial papers
- ✅ 232 scRNA-seq papers

### Features
- ✅ UMAP/tSNE/PCA embeddings
- ✅ Gene expression queries (normalized counts)
- ✅ Metadata filtering
- ✅ Gene search with autocomplete
- ✅ DEG analysis (t-test with FDR correction)
- ✅ Module score calculation (3 methods)
- ✅ Automatic caching
- ✅ Fast response times (<100ms cached)

## 🎓 Learning Path

### For Frontend Developers
1. Start with [Quick Start Guide](./H5AD_QUICK_START.md)
2. Read [H5AD API Guide](./H5AD_API_GUIDE.md) for visualization
3. Check [Datasets API Guide](./DATASETS_API_GUIDE.md) for data browsing
4. Explore [DEG Guide](./DEG_ANALYSIS_GUIDE.md) for analysis

### For Bioinformaticians
1. Read [Architecture Summary](./H5AD_ARCHITECTURE_SUMMARY.md)
2. Deep dive into [DEG Analysis Guide](./DEG_ANALYSIS_GUIDE.md)
3. Explore [Module Score Guide](./MODULE_SCORE_GUIDE.md)
4. Check [Architecture Details](./ARCHITECTURE_H5AD.md) for implementation

### For Data Scientists
1. Start with [DEG Analysis Guide](./DEG_ANALYSIS_GUIDE.md)
2. Read [Module Score Guide](./MODULE_SCORE_GUIDE.md)
3. Check [H5AD API Guide](./H5AD_API_GUIDE.md) for data access
4. Review [Architecture Summary](./H5AD_ARCHITECTURE_SUMMARY.md)

## 🔍 Search Documentation

To find specific topics:

- **Endpoints**: Check the [Quick Links](#quick-links) table
- **Code examples**: All guides include frontend examples
- **Error handling**: See Troubleshooting sections in each guide
- **Performance**: See [Architecture Summary](./H5AD_ARCHITECTURE_SUMMARY.md)

## 📝 API Conventions

### Response Format
All endpoints return JSON with consistent structure:
- **Success**: HTTP 200 with data
- **Not Found**: HTTP 404 with error message
- **Bad Request**: HTTP 400 with error details
- **Server Error**: HTTP 500 with error message

### Naming Conventions
- `dataset_id`: Unique identifier for dataset
- `gene`: Gene symbol (e.g., "APOE")
- `cell_type`: Cell type annotation
- `embedding_type`: "umap", "tsne", or "pca"

### Data Types
- **Normalized counts**: Default for visualization (recommended)
- **Raw counts**: Available via `use_raw: true` parameter
- **Log-transformed**: Stored in h5ad `.X` matrix

## 🐛 Troubleshooting

### Common Issues

**"Dataset not found"**
- Check dataset ID spelling
- Use `/datasets/` to list available datasets
- See: [Datasets API Guide](./DATASETS_API_GUIDE.md)

**"Gene not found"**
- Check gene symbol (case-sensitive)
- Use `/h5ad/{id}/genes/search` for autocomplete
- See: [H5AD API Guide](./H5AD_API_GUIDE.md#search-genes)

**"Slow response"**
- First access computes and caches (5-30s)
- Subsequent access is fast (<100ms)
- See: [Architecture Summary](./H5AD_ARCHITECTURE_SUMMARY.md#performance)

## 🤝 Contributing

### Adding New Datasets
1. Place h5ad file in `h5ad/raw/`
2. Run preprocessing: `python scripts/preprocess_h5ad.py --h5ad h5ad/raw/your_file.h5ad`
3. Restart server: `pm2 restart sskind-backend`

### Documentation Updates
All documentation is in Markdown format in this `docs/` directory.

## 📞 Support

- **Interactive Docs**: http://localhost:9117/sskind-backend/docs
- **GitHub Issues**: (your repository)
- **Email**: (your contact)

## 🔄 Version History

### Current Version
- H5AD data serving with caching
- DEG analysis (between & within datasets)
- Module score calculation
- Full dataset metadata API
- Comprehensive documentation

### Future Plans
- [ ] Batch gene expression queries
- [ ] Clustering API
- [ ] Cell-cell communication analysis
- [ ] Trajectory analysis
- [ ] Gene set enrichment analysis

## 📄 License

(Your license information)

---

**Last Updated**: November 2024

**API Version**: 1.0.0

**Documentation Version**: 1.0.0

