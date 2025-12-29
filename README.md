# HypoMap Backend

A FastAPI backend for single-cell RNA-seq analysis with AI-powered hypothesis generation.

## Features

- **H5AD Support** - Single-cell data in AnnData format
- **UMAP/tSNE/PCA** - Precomputed embeddings with caching
- **Gene Expression** - Normalized and raw count access
- **DEG Analysis** - Differential expression with t-test and FDR
- **CCC Network** - Cell-cell communication data from CellChat
- **Regulon Network** - Transcription factor regulatory networks
- **AI Assistant** - Ollama-powered LLM for biological Q&A
- **Smart Caching** - Fast response times (<100ms)

## Tech Stack

- FastAPI
- SQLite
- AnnData/Scanpy
- Ollama API

## Getting Started

### Prerequisites

- Python 3.11+ (via micromamba)
- Ollama server (for AI features)

### Development

```bash
# Create environment
micromamba create -n hypomap python=3.11 -c conda-forge -y
micromamba activate hypomap

# Install dependencies
pip install -r requirements.txt

# Run development server (port 9120)
python main.py
```

API docs: http://localhost:9120/hypomap-backend/docs

### Production Deployment

See [docs/DEPLOYMENT.md](./docs/DEPLOYMENT.md) for full deployment guide.

```bash
# Start with PM2
pm2 start ecosystem.config.js

# Management commands
pm2 status                   # Check status
pm2 logs hypomap-backend     # View logs
pm2 restart hypomap-backend  # Restart
pm2 stop hypomap-backend     # Stop
```

## Project Structure

```
app/
├── api/
│   ├── h5ad_endpoints.py    # Data visualization APIs
│   ├── ai_endpoints.py      # AI chat endpoints
│   └── deg_endpoints.py     # Differential expression
├── services/
│   ├── h5ad_service.py      # H5AD data handling
│   ├── ollama_service.py    # Ollama AI integration
│   └── deg_service.py       # DEG analysis
├── config.py                # App configuration
└── app.py                   # FastAPI application

h5ad/
├── raw/                     # Original .h5ad files
└── precomputed/             # Processed data (JSON, SQLite)
    └── {dataset_id}/
        ├── ccc/             # CellChat data
        └── regulon/         # Regulon network data

docs/
└── DEPLOYMENT.md            # Production deployment guide
```

## API Endpoints

### Data Visualization
```bash
# Get UMAP coordinates
GET /hypomap-backend/api/v1/h5ad/{dataset_id}/embedding/umap

# Get gene expression
GET /hypomap-backend/api/v1/h5ad/{dataset_id}/expression/{gene}

# Get CCC network
GET /hypomap-backend/api/v1/h5ad/{dataset_id}/ccc/network

# Get regulon network
GET /hypomap-backend/api/v1/h5ad/{dataset_id}/regulon/network?cluster=0
```

### AI Assistant
```bash
# Ask a question
POST /hypomap-backend/api/v1/ai/ask
{
  "dataset_id": "hypomap_demo",
  "cluster": "0",
  "question": "What is the function of key genes in this cluster?",
  "model": "qwen3:30b"
}

# Get available models
GET /hypomap-backend/api/v1/ai/models
```

## Configuration

Environment variables in `.env`:
```env
PORT=9120
HOST=0.0.0.0
DEBUG=False
GLOBAL_PREFIX=/hypomap-backend
API_PREFIX=/api/v1
```

## Data Preprocessing

```bash
# Preprocess H5AD file
python scripts/preprocess_h5ad.py h5ad/raw/dataset.h5ad
```

This generates:
- Precomputed embeddings (UMAP, tSNE, PCA) as JSON
- Metadata cache
- Gene statistics in SQLite

## Related

- [HypoMap Frontend](https://github.com/Wang-Cankun/hypomap-frontend) - Vue 3 frontend

## License

MIT
