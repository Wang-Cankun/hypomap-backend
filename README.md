# ssKIND Backend

A modern FastAPI backend for serving and analyzing single-cell RNA-seq and spatial transcriptomics data, featuring comprehensive APIs for data visualization, gene expression analysis, and differential expression.

## 📚 Documentation

**Complete API documentation is available in the [`docs/`](./docs/) folder:**

- **[📖 Full Documentation Index](./docs/README.md)** - Start here!
- **[⚡ Quick Start Guide](./docs/H5AD_QUICK_START.md)** - Get started in 30 seconds
- **[📊 Datasets API](./docs/DATASETS_API_GUIDE.md)** - Browse and search datasets
- **[🧬 H5AD Data API](./docs/H5AD_API_GUIDE.md)** - Visualizations and gene expression
- **[📈 DEG Analysis](./docs/DEG_ANALYSIS_GUIDE.md)** - Differential expression analysis
- **[🎯 Module Score](./docs/MODULE_SCORE_GUIDE.md)** - Gene signature scoring
- **[🏗️ Architecture](./docs/H5AD_ARCHITECTURE_SUMMARY.md)** - System design and performance

## Features

- 🚀 **FastAPI** - Modern, fast web framework for building APIs
- 🗄️ **SQLite** - Lightweight, serverless database
- 🧬 **H5AD Support** - Single-cell data in AnnData format
- 📊 **Gene Expression** - Normalized and raw count access
- 🎨 **UMAP/tSNE/PCA** - Precomputed embeddings with caching
- 📈 **DEG Analysis** - Statistical differential expression (t-test, FDR)
- 🎯 **Module Scores** - Gene signature and pathway activity
- 🔍 **Gene Search** - Autocomplete and fuzzy matching
- ⚡ **Smart Caching** - Fast response times (<100ms)
- 📚 **Auto-generated API docs** - Interactive docs at `/docs`
- 🔄 **CORS support** - Frontend-friendly
- 🏗️ **Modular structure** - Clean, organized codebase

## Project Structure

```
├── app/                    # Application code
│   ├── api/                # API endpoints
│   │   ├── endpoints.py    # Main API routes
│   │   ├── h5ad_endpoints.py
│   │   ├── deg_endpoints.py
│   │   └── api_metadata.py
│   ├── services/          # Business logic
│   │   ├── h5ad_service.py # H5AD data handling
│   │   └── deg_service.py  # DEG analysis
│   ├── models.py           # SQLAlchemy models
│   ├── schemas.py         # Pydantic schemas
│   ├── crud.py            # Database operations
│   ├── database.py        # Database configuration
│   └── config.py          # App configuration
├── scripts/                # Import and utility scripts
│   ├── import_csv.py       # scRNA-seq data import
│   ├── import_papers.py    # Papers import
│   ├── import_spatial_data.py
│   ├── import_spatial_papers.py
│   ├── preprocess_h5ad.py  # H5AD preprocessing
│   └── compute_embeddings.py  # Embedding computation
├── tests/                  # Test scripts
│   ├── test_sqlite_performance.py
│   ├── test_api_performance.py
│   ├── test_performance_comparison.py
│   └── test_spatial_api.py
├── docs/                   # Documentation
│   ├── README.md           # Documentation index
│   ├── H5AD_QUICK_START.md
│   ├── API_SUMMARY.md
│   ├── PERFORMANCE_REPORT.md
│   └── ... (see docs/README.md)
├── h5ad/                   # H5AD data storage
│   ├── raw/                # Original .h5ad files
│   └── precomputed/        # Processed data (JSON, SQLite)
├── data/                   # CSV data files
├── main.py                 # Application entry point
├── scripts/                 # Utility scripts
│   ├── preprocess_h5ad.py  # H5AD preprocessing
│   ├── compute_embeddings.py  # Embedding computation
│   ├── import_csv.py       # Data import scripts
│   ├── import_papers.py
│   ├── import_spatial_data.py
│   └── import_spatial_papers.py
├── import_all_data.py      # Combined import script
├── requirements.txt        # Python dependencies
└── README.md              # This file
```

## CSV Import

The application includes functionality to import scRNA-seq datasets from CSV files. The latest schema matches the following headers:

```text
Dataset_id,Public_dataset_id,Normation,Pubmed_id,Disease,Status,Control,Species,Brain_region,Treatment,Sex,Stage,Age,Zhu_Age,N_cells,Protocol,Methodology,Model
```

### Import Scripts

All import scripts are located in the `scripts/` directory:

**scRNA-seq Data Import:**
```bash
python scripts/import_csv.py --csv "data/ssKIND - scRNAseq_data_for_ssKIND.csv"
```

**Replace All Existing Data:**
```bash
# Drop and recreate table before import
python scripts/import_csv.py --replace --csv "data/ssKIND - scRNAseq_data_for_ssKIND.csv"
```

**Options:**
- `--replace`: Drops and recreates the table, then imports fresh data
- `--batch-size N`: Set batch size for database inserts (default: 100)

## Running the Application

```bash
python main.py
# or
python run.py
```

Docs: `http://localhost:8000/sskind-backend/docs`

## API Endpoints (Datasets)

- Get all datasets:
  ```bash
  curl "http://localhost:8000/sskind-backend/api/v1/datasets/"
  ```
- Get dataset stats:
  ```bash
  curl "http://localhost:8000/sskind-backend/api/v1/datasets/stats/"
  ```
- Search by disease:
  ```bash
  curl "http://localhost:8000/sskind-backend/api/v1/datasets/search/disease/AD"
  ```

## Notes

- If you change CSV headers, update `app/models.py` and `app/schemas.py` to match.
- Use `--replace` to safely refresh the table when the schema changes.

## Data Import

### Individual Import Scripts

All import scripts are in the `scripts/` directory:

**scRNA-seq Papers:**
```bash
python scripts/import_papers.py --csv "data/ssKIND - scRNAseq_paper_for_ssKIND.csv"
python scripts/import_papers.py --replace --csv "data/ssKIND - scRNAseq_paper_for_ssKIND.csv"
```

**Spatial Data:**
```bash
python scripts/import_spatial_data.py --csv "data/ssKIND - Spatial data for ssKIND.csv"
```

**Spatial Papers:**
```bash
python scripts/import_spatial_papers.py --csv "data/ssKIND - Spatial paper for ssKIND.csv"
```

### Combined Import

Import all data types at once:

```bash
python import_all_data.py --replace
```

## Papers API Endpoints

- List papers:
  ```bash
  curl "http://localhost:8000/sskind-backend/api/v1/papers/"
  ```
- Get one paper by id:
  ```bash
  curl "http://localhost:8000/sskind-backend/api/v1/papers/AD001"
  ```
- Search:

  ```bash
  curl "http://localhost:8000/sskind-backend/api/v1/papers/search/disease/AD"
  curl "http://localhost:8000/sskind-backend/api/v1/papers/search/species/Human"
  curl "http://localhost:8000/sskind-backend/api/v1/papers/search/brain-region/Prefrontal%20cortex"
  ```

- **Spatial Datasets API**:

  ```bash
  # Get all spatial datasets
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/"

  # Get specific spatial dataset
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/ST024001"

  # Search spatial datasets
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/search/disease/AD"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/search/species/human"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/search/brain-region/frontal%20cortex"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/search/methodology/Visium"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/search/study/ST024"

  # Get spatial datasets statistics
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-datasets/stats/"
  ```

- **Spatial Papers API**:

  ```bash
  # Get all spatial papers
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/"

  # Get specific spatial paper
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/ST024"

  # Search spatial papers
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/search/disease/AD"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/search/species/human"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/search/brain-region/frontal%20cortex"
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/search/methodology/Visium"

  # Get spatial papers statistics
  curl "http://localhost:8000/sskind-backend/api/v1/spatial-papers/stats/"
  ```

### Database Tables

The import scripts create the following database tables:

- `scrna_datasets`: scRNA-seq dataset information
- `scrna_papers`: scRNA-seq research papers
- `spatial_datasets`: Spatial transcriptomics dataset information
- `spatial_papers`: Spatial transcriptomics research papers

### ID Structure for Spatial Data

- **`data_id`**: Study/project identifier (e.g., ST024) - multiple datasets can share the same data_id
- **`dataset_id`**: Unique sample identifier (e.g., ST024001, ST024002) - each row has a unique dataset_id
- Use `dataset_id` to get specific datasets, `data_id` to get all datasets from a study

## Production Deployment with PM2

The project includes PM2 configuration for production deployment:

### PM2 Setup

1. **Install PM2** (if not already installed):

   ```bash
   npm install -g pm2
   ```

2. **Start the application**:

   ```bash
   # Using the management script (recommended)
   ./pm2-manage.sh start

   # Or directly with PM2
   pm2 start ecosystem.config.js --env production
   ```

3. **Management commands**:
   ```bash
   ./pm2-manage.sh status    # Check application status
   ./pm2-manage.sh logs      # View logs
   ./pm2-manage.sh restart   # Restart application
   ./pm2-manage.sh stop      # Stop application
   ./pm2-manage.sh delete    # Remove from PM2
   ```

### PM2 Features Configured

- **Port**: 9117 (as requested)
- **Host**: 0.0.0.0 (accepts connections from any IP)
- **Python Interpreter**: `~/micromamba/envs/sskind/bin/python`
- **Auto-restart**: Enabled with smart restart policies
- **Memory limit**: 1GB (restarts if exceeded)
- **Logging**: Separate error, output, and combined log files
- **Zero-downtime reload**: Use `pm2 reload sskind-backend`

### Environment Modes

- **Production**: `pm2 start ecosystem.config.js --env production`
- **Development**: `pm2 start ecosystem.config.js --env development`

The application will be available at: `http://localhost:9117/sskind-backend/`

---

## 🚀 Quick API Reference

### Common Endpoints

```bash
# Interactive API documentation
http://localhost:9117/sskind-backend/docs

# List datasets
curl "http://localhost:9117/sskind-backend/api/v1/datasets/"

# Get UMAP coordinates
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/embedding/umap"

# Get gene expression
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/expression/APOE"

# Calculate module score
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/module-score" \
  -H "Content-Type: application/json" \
  -d '{"gene_list": ["APOE", "APP", "MAPT"]}'

# Differential expression analysis
curl -X POST "http://localhost:9117/sskind-backend/api/v1/deg/between-datasets" \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id1": "AD093044",
    "dataset_id2": "AD093045",
    "cell_type": "Microglia"
  }'
```

### Frontend Integration Example

```javascript
// Get UMAP + gene expression for visualization
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/plot-data?embedding=umap&genes=APOE,APP&metadata=cell_type"
);
const data = await response.json();

// data.embedding: UMAP coordinates
// data.gene_expression: Gene expression values
// data.metadata: Cell type annotations
```

**📖 For complete documentation, examples, and advanced usage, see the [`docs/`](./docs/) folder.**

## Testing

Test scripts are located in the `tests/` directory:

```bash
# Performance tests
python tests/test_sqlite_performance.py
python tests/test_api_performance.py
python tests/test_performance_comparison.py

# API tests
python tests/test_spatial_api.py
```

See [`tests/README.md`](./tests/README.md) for detailed test documentation.

## H5AD Data Preprocessing

To preprocess H5AD files and generate optimized data:

```bash
python scripts/preprocess_h5ad.py --h5ad h5ad/raw/dataset_id.h5ad --dataset-id dataset_id
```

To compute embeddings (PCA, UMAP, tSNE) for h5ad files that don't have them:

```bash
python scripts/compute_embeddings.py --h5ad h5ad/raw/dataset_id.h5ad
```

This creates:
- Precomputed embeddings (UMAP, tSNE, PCA) as JSON
- Metadata cache
- Gene statistics in SQLite database (`data.db`) for fast search
- Individual gene expression files

See [`docs/PERFORMANCE_REPORT.md`](./docs/PERFORMANCE_REPORT.md) for performance details.
