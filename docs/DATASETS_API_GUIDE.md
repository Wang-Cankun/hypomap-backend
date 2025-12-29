# Datasets API Guide

## Overview

The ssKIND backend provides comprehensive APIs for managing and querying single-cell RNA-seq and spatial transcriptomics datasets. This guide covers all dataset-related endpoints including metadata, papers, and dataset information.

## Base URL

```
http://localhost:9117/sskind-backend/api/v1
```

## Table of Contents

1. [scRNA-seq Datasets](#scrna-seq-datasets)
2. [scRNA-seq Papers](#scrna-seq-papers)
3. [Spatial Datasets](#spatial-datasets)
4. [Spatial Papers](#spatial-papers)
5. [Search & Filtering](#search--filtering)
6. [Statistics](#statistics)

---

## scRNA-seq Datasets

### List All Datasets

**GET** `/datasets/`

**Query Parameters:**

- `skip` (int): Number of records to skip (pagination)
- `limit` (int): Maximum records to return

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/datasets/?limit=10&skip=0"
```

**Response:**

```json
[
  {
    "id": 1,
    "dataset_id": "AD093044",
    "public_dataset_id": "GSE123456",
    "disease": "AD",
    "species": "Human",
    "brain_region": "Prefrontal cortex",
    "n_cells": 7952,
    "methodology": "10x Genomics",
    "created_at": "2024-01-01T00:00:00Z"
  }
]
```

### Get Specific Dataset

**GET** `/datasets/{dataset_id}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/datasets/AD093044"
```

### Get Dataset Statistics

**GET** `/datasets/stats/`

Returns summary statistics across all datasets.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/datasets/stats/"
```

**Response:**

```json
{
  "total_datasets": 3620,
  "total_cells": 2500000,
  "diseases": ["AD", "PD", "ALS", "HD"],
  "species": ["Human", "Mouse"],
  "methodologies": ["10x Genomics", "Smart-seq2", "Drop-seq"],
  "by_disease": {
    "AD": 1500,
    "PD": 800,
    "ALS": 500
  }
}
```

### Search Datasets

#### By Disease

**GET** `/datasets/search/disease/{disease}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/datasets/search/disease/AD"
```

#### By Species

**GET** `/datasets/search/species/{species}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/datasets/search/species/Human"
```

#### By Brain Region

**GET** `/datasets/search/brain-region/{region}`

**Example:**

```bash
# Note: URL encode spaces
curl "http://localhost:9117/sskind-backend/api/v1/datasets/search/brain-region/Prefrontal%20cortex"
```

---

## scRNA-seq Papers

### List All Papers

**GET** `/papers/`

**Query Parameters:**

- `skip` (int): Pagination offset
- `limit` (int): Maximum records

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/papers/?limit=10"
```

**Response:**

```json
[
  {
    "id": 1,
    "paper_id": "AD001",
    "title": "Single-cell analysis reveals...",
    "author": "Smith J, Doe A",
    "journal": "Nature",
    "pubmed_id": "12345678",
    "disease": "AD",
    "species": "Human",
    "methodology": "10x Genomics",
    "citation": "Smith et al. Nature 2023",
    "doi": "10.1038/...",
    "abstract": "We performed..."
  }
]
```

### Get Specific Paper

**GET** `/papers/{paper_id}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/papers/AD001"
```

### Search Papers

#### By Disease

**GET** `/papers/search/disease/{disease}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/papers/search/disease/AD"
```

#### By Species

**GET** `/papers/search/species/{species}`

#### By Brain Region

**GET** `/papers/search/brain-region/{region}`

### Get Paper Statistics

**GET** `/papers/stats/`

---

## Spatial Datasets

### List All Spatial Datasets

**GET** `/spatial-datasets/`

**Query Parameters:**

- `skip` (int): Pagination offset
- `limit` (int): Maximum records

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/spatial-datasets/?limit=10"
```

**Response:**

```json
[
  {
    "id": 1,
    "data_id": "ST001",
    "dataset_id": "ST001001",
    "public_dataset_id": "GSE...",
    "disease": "AD",
    "species": "Human",
    "brain_region": "Hippocampus",
    "methodology": "Visium",
    "n_cells": null,
    "slide_id": "151507"
  }
]
```

### Get Specific Spatial Dataset

**GET** `/spatial-datasets/{dataset_id}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/spatial-datasets/ST001001"
```

### Search Spatial Datasets

#### By Disease

**GET** `/spatial-datasets/search/disease/{disease}`

#### By Species

**GET** `/spatial-datasets/search/species/{species}`

#### By Brain Region

**GET** `/spatial-datasets/search/brain-region/{region}`

#### By Methodology

**GET** `/spatial-datasets/search/methodology/{methodology}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/spatial-datasets/search/methodology/Visium"
```

#### By Study ID

**GET** `/spatial-datasets/search/study/{study_id}`

Get all datasets from a specific study (same `data_id`).

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/spatial-datasets/search/study/ST001"
```

### Get Spatial Dataset Statistics

**GET** `/spatial-datasets/stats/`

**Response:**

```json
{
  "total_datasets": 1380,
  "diseases": ["AD", "PD", "Control"],
  "species": ["Human", "Mouse"],
  "methodologies": ["Visium", "Slide-seq", "MERFISH"],
  "by_methodology": {
    "Visium": 800,
    "Slide-seq": 400,
    "MERFISH": 180
  }
}
```

---

## Spatial Papers

### List All Spatial Papers

**GET** `/spatial-papers/`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/spatial-papers/"
```

### Get Specific Spatial Paper

**GET** `/spatial-papers/{paper_id}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/spatial-papers/ST001"
```

### Search Spatial Papers

#### By Disease

**GET** `/spatial-papers/search/disease/{disease}`

#### By Species

**GET** `/spatial-papers/search/species/{species}`

#### By Brain Region

**GET** `/spatial-papers/search/brain-region/{region}`

#### By Methodology

**GET** `/spatial-papers/search/methodology/{methodology}`

### Get Spatial Paper Statistics

**GET** `/spatial-papers/stats/`

---

## Frontend Integration Examples

### 1. Load and Display Datasets

```javascript
async function loadDatasets() {
  const response = await fetch(
    "http://localhost:9117/sskind-backend/api/v1/datasets/?limit=100"
  );
  const datasets = await response.json();

  // Display in table
  datasets.forEach((dataset) => {
    console.log(
      `${dataset.dataset_id}: ${dataset.disease}, ${dataset.n_cells} cells`
    );
  });
}
```

### 2. Search with Filters

```javascript
async function searchDatasets(disease, species) {
  // Search by disease
  const diseaseResponse = await fetch(
    `http://localhost:9117/sskind-backend/api/v1/datasets/search/disease/${disease}`
  );
  const diseaseDatasets = await diseaseResponse.json();

  // Filter by species on frontend
  const filtered = diseaseDatasets.filter((d) => d.species === species);

  return filtered;
}

// Usage
const adHumanDatasets = await searchDatasets("AD", "Human");
```

### 3. Dataset Browser with Pagination

```javascript
class DatasetBrowser {
  constructor() {
    this.page = 0;
    this.pageSize = 20;
  }

  async loadPage() {
    const skip = this.page * this.pageSize;
    const response = await fetch(
      `http://localhost:9117/sskind-backend/api/v1/datasets/?skip=${skip}&limit=${this.pageSize}`
    );
    return await response.json();
  }

  async nextPage() {
    this.page++;
    return await this.loadPage();
  }

  async prevPage() {
    if (this.page > 0) {
      this.page--;
    }
    return await this.loadPage();
  }
}

// Usage
const browser = new DatasetBrowser();
const datasets = await browser.loadPage();
```

### 4. Get Dataset with Associated Paper

```javascript
async function getDatasetWithPaper(datasetId) {
  // Get dataset
  const datasetResponse = await fetch(
    `http://localhost:9117/sskind-backend/api/v1/datasets/${datasetId}`
  );
  const dataset = await datasetResponse.json();

  // Get associated paper using pubmed_id
  const paperResponse = await fetch(
    `http://localhost:9117/sskind-backend/api/v1/papers/`
  );
  const papers = await paperResponse.json();
  const paper = papers.find((p) => p.pubmed_id === dataset.pubmed_id);

  return {
    dataset,
    paper,
  };
}
```

### 5. Statistics Dashboard

```javascript
async function loadStatistics() {
  // Get scRNA-seq stats
  const scrnaStatsResponse = await fetch(
    "http://localhost:9117/sskind-backend/api/v1/datasets/stats/"
  );
  const scrnaStats = await scrnaStatsResponse.json();

  // Get spatial stats
  const spatialStatsResponse = await fetch(
    "http://localhost:9117/sskind-backend/api/v1/spatial-datasets/stats/"
  );
  const spatialStats = await spatialStatsResponse.json();

  return {
    scrna: scrnaStats,
    spatial: spatialStats,
    total: {
      datasets: scrnaStats.total_datasets + spatialStats.total_datasets,
      cells: scrnaStats.total_cells,
    },
  };
}

// Display
const stats = await loadStatistics();
console.log(`Total Datasets: ${stats.total.datasets}`);
console.log(`Total Cells: ${stats.total.cells.toLocaleString()}`);
```

### 6. Interactive Search UI

```javascript
class DatasetSearch {
  constructor() {
    this.filters = {
      disease: null,
      species: null,
      brainRegion: null,
    };
  }

  async search() {
    let results = [];

    // Start with disease filter if set
    if (this.filters.disease) {
      const response = await fetch(
        `http://localhost:9117/sskind-backend/api/v1/datasets/search/disease/${this.filters.disease}`
      );
      results = await response.json();
    } else {
      // Get all
      const response = await fetch(
        "http://localhost:9117/sskind-backend/api/v1/datasets/"
      );
      results = await response.json();
    }

    // Apply additional filters on frontend
    if (this.filters.species) {
      results = results.filter((d) => d.species === this.filters.species);
    }

    if (this.filters.brainRegion) {
      results = results.filter(
        (d) =>
          d.brain_region && d.brain_region.includes(this.filters.brainRegion)
      );
    }

    return results;
  }

  setFilter(type, value) {
    this.filters[type] = value;
  }

  clearFilters() {
    this.filters = { disease: null, species: null, brainRegion: null };
  }
}

// Usage
const search = new DatasetSearch();
search.setFilter("disease", "AD");
search.setFilter("species", "Human");
const results = await search.search();
```

### 7. Spatial Dataset by Study

```javascript
async function getStudyDatasets(studyId) {
  // Get all datasets from a study
  const response = await fetch(
    `http://localhost:9117/sskind-backend/api/v1/spatial-datasets/search/study/${studyId}`
  );
  const datasets = await response.json();

  console.log(`Study ${studyId} has ${datasets.length} datasets`);

  return datasets;
}

// Get all slides from ST001 study
const st001_datasets = await getStudyDatasets("ST001");
```

## Data Fields Reference

### scRNA-seq Dataset Fields

| Field               | Type     | Description                             |
| ------------------- | -------- | --------------------------------------- |
| `id`                | int      | Database primary key                    |
| `dataset_id`        | string   | Unique dataset identifier               |
| `public_dataset_id` | string   | Public repository ID (e.g., GSE number) |
| `normation`         | string   | Normalization method                    |
| `pubmed_id`         | string   | PubMed ID                               |
| `disease`           | string   | Disease type (AD, PD, etc.)             |
| `status`            | string   | Disease status                          |
| `control`           | string   | Control type                            |
| `species`           | string   | Species (Human, Mouse)                  |
| `brain_region`      | string   | Brain region studied                    |
| `treatment`         | string   | Treatment applied                       |
| `sex`               | string   | Sex of samples                          |
| `stage`             | string   | Disease stage                           |
| `age`               | string   | Age of samples                          |
| `n_cells`           | int      | Number of cells                         |
| `protocol`          | string   | Experimental protocol                   |
| `methodology`       | string   | Sequencing methodology                  |
| `model`             | string   | Disease model                           |
| `created_at`        | datetime | Creation timestamp                      |

### Spatial Dataset Fields

| Field               | Type   | Description                               |
| ------------------- | ------ | ----------------------------------------- |
| `id`                | int    | Database primary key                      |
| `data_id`           | string | Study/project identifier (e.g., ST001)    |
| `dataset_id`        | string | Unique sample identifier (e.g., ST001001) |
| `data_path`         | string | Path to data files                        |
| `public_dataset_id` | string | Public repository ID                      |
| `slide_id`          | string | Slide identifier                          |
| `disease`           | string | Disease type                              |
| `status`            | string | Disease status                            |
| `species`           | string | Species                                   |
| `brain_region`      | string | Brain region                              |
| `treatment`         | string | Treatment                                 |
| `sex`               | string | Sex                                       |
| `stage`             | string | Stage                                     |
| `age`               | string | Age                                       |
| `n_cells`           | int    | Number of spots/cells                     |
| `protocol`          | string | Protocol                                  |
| `methodology`       | string | Spatial method (Visium, Slide-seq, etc.)  |
| `model`             | string | Disease model                             |

### Paper Fields

| Field            | Type   | Description             |
| ---------------- | ------ | ----------------------- |
| `id`             | int    | Database primary key    |
| `paper_id`       | string | Unique paper identifier |
| `title`          | string | Paper title             |
| `author`         | string | Authors                 |
| `journal`        | string | Journal name            |
| `citation`       | string | Full citation           |
| `abstract`       | string | Paper abstract          |
| `pubmed_id`      | string | PubMed ID               |
| `doi`            | string | DOI                     |
| `disease`        | string | Disease studied         |
| `species`        | string | Species                 |
| `methodology`    | string | Methodology used        |
| `public_data_id` | string | Data repository link    |

## Best Practices

### 1. Pagination

Always use pagination for large datasets:

```javascript
// Good: Paginated
const response = await fetch("/api/v1/datasets/?limit=100&skip=0");

// Avoid: Loading all at once (slow for large datasets)
const response = await fetch("/api/v1/datasets/");
```

### 2. Caching

Cache frequently accessed data:

```javascript
class DatasetCache {
  constructor() {
    this.cache = new Map();
    this.ttl = 5 * 60 * 1000; // 5 minutes
  }

  async get(datasetId) {
    const cached = this.cache.get(datasetId);
    if (cached && Date.now() - cached.timestamp < this.ttl) {
      return cached.data;
    }

    const response = await fetch(`/api/v1/datasets/${datasetId}`);
    const data = await response.json();

    this.cache.set(datasetId, {
      data,
      timestamp: Date.now(),
    });

    return data;
  }
}
```

### 3. Error Handling

```javascript
async function getDatasetSafe(datasetId) {
  try {
    const response = await fetch(`/api/v1/datasets/${datasetId}`);

    if (!response.ok) {
      if (response.status === 404) {
        throw new Error(`Dataset ${datasetId} not found`);
      }
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }

    return await response.json();
  } catch (error) {
    console.error("Error fetching dataset:", error);
    throw error;
  }
}
```

### 4. Search Optimization

For complex searches, filter on backend then refine on frontend:

```javascript
// 1. Backend filter (fast, reduces data transfer)
const response = await fetch("/api/v1/datasets/search/disease/AD");
const adDatasets = await response.json();

// 2. Frontend filter (for fields without backend search)
const filtered = adDatasets.filter(
  (d) => d.n_cells > 1000 && d.methodology === "10x Genomics"
);
```

## API Rate Limits

- No rate limits currently implemented
- Recommended: Max 100 requests per minute per client
- Use pagination and caching to minimize requests

## Troubleshooting

### No results returned

- Check spelling (case-sensitive for some fields)
- Use `/stats/` endpoints to see available values
- Try broader search terms

### Slow response

- Use pagination (`limit` parameter)
- Cache results
- Search by specific criteria instead of loading all

### 404 Not Found

- Verify dataset/paper ID is correct
- Check if dataset exists using list endpoint
- Use search endpoint to find correct ID

## Related Documentation

- [H5AD Data API](./H5AD_API_GUIDE.md) - For analyzing specific datasets
- [DEG Analysis](./DEG_ANALYSIS_GUIDE.md) - Differential expression
- [Module Score](./MODULE_SCORE_GUIDE.md) - Gene signature scoring

## API Documentation

Interactive API documentation:

```
http://localhost:9117/sskind-backend/docs
```

Browse all endpoints with interactive testing capabilities.
