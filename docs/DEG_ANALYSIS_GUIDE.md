# DEG (Differential Expression Gene) Analysis API Guide

## Overview

The ssKIND backend now supports comprehensive DEG analysis for single-cell data with two modes:

1. **Between Datasets**: Compare two different datasets (e.g., control vs disease)
2. **Within Dataset**: Compare two groups within the same dataset (e.g., cell type A vs cell type B)

## Features

✅ T-test based statistical analysis
✅ Multiple testing correction (FDR/Benjamini-Hochberg)
✅ Cell type filtering (single or multiple cell types per group)
✅ Flexible metadata filtering
✅ Automatic caching of results
✅ Top N gene filtering
✅ Comprehensive statistics (fold change, p-values, expression percentages)
✅ Automatic log1p transformation when raw counts detected

## API Endpoints

### 1. DEG Between Two Datasets

**Endpoint:** `POST /api/v1/deg/between-datasets`

**Use Cases:**

- Compare control vs disease datasets
- Compare different experimental conditions
- Compare same cell types across datasets

**Example Request:**

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/deg/between-datasets" \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id1": "AD093044",
    "dataset_id2": "AD093044",
    "cell_types_group1": ["Astrocyte", "CGE interneuron"],
    "cell_types_group2": [
      "Cerebellar inhibitory",
      "Committed oligodendrocyte precursor",
      "Deep-layer corticothalamic and 6b",
      "Deep-layer intratelencephalic",
      "Deep-layer near-projecting",
      "Fibroblast",
      "LAMP5-LHX6 and Chandelier",
      "MGE interneuron",
      "Microglia",
      "Miscellaneous",
      "Oligodendrocyte",
      "Oligodendrocyte precursor",
      "Splatter",
      "Thalamic excitatory",
      "Unknown",
      "Upper rhombic lip",
      "Upper-layer intratelencephalic",
      "Vascular"
    ],
    "min_pct": 0.1,
    "logfc_threshold": 0.25,
    "p_value_threshold": 0.05,
    "top_n": 200
  }'
```

**Parameters:**

| Parameter           | Type   | Required | Default | Description                               |
| ------------------- | ------ | -------- | ------- | ----------------------------------------- |
| `dataset_id1`       | string | ✅       | -       | First dataset (group 1)                   |
| `dataset_id2`       | string | ✅       | -       | Second dataset (group 2)                  |
| `cell_type`         | string | ❌       | null    | Cell type(s) for group 1 (comma-separated) |
| `cell_type2`        | string | ❌       | null    | Cell type(s) for group 2 (comma-separated) |
| `cell_types_group1` | array  | ❌       | null    | Array of cell types for group 1           |
| `cell_types_group2` | array  | ❌       | null    | Array of cell types for group 2           |
| `metadata_filters1` | object | ❌       | null    | Additional filters for dataset 1          |
| `metadata_filters2` | object | ❌       | null    | Additional filters for dataset 2          |
| `min_pct`           | float  | ❌       | 0.1     | Minimum fraction of cells expressing gene |
| `logfc_threshold`   | float  | ❌       | 0.25    | Minimum log2 fold change for significance |
| `p_value_threshold` | float  | ❌       | 0.05    | Adjusted p-value threshold                |
| `top_n`             | int    | ❌       | null    | Return only top N genes                   |

> **Tip:** You can pass multiple cell types either as arrays (`"cell_types_group1": ["Astrocyte","CGE interneuron"]`) or as comma-separated strings (`"cell_type": "Astrocyte,CGE interneuron"`). The backend automatically normalizes the inputs.

**Response:**

```json
{
  "comparison": {
    "group1": {
      "dataset_id": "AD093044",
      "n_cells": 401,
      "cell_type": "Microglia",
      "filters": null
    },
    "group2": {
      "dataset_id": "AD093045",
      "n_cells": 454,
      "cell_type": "Microglia",
      "filters": null
    }
  },
  "parameters": {
    "min_pct": 0.1,
    "logfc_threshold": 0.5,
    "p_value_threshold": 0.05
  },
  "summary": {
    "total_genes_tested": 3949,
    "significant_genes": 788,
    "upregulated_in_group1": 183,
    "upregulated_in_group2": 605
  },
  "genes": [
    {
      "gene": "MALAT1",
      "log2_fold_change": -1.35,
      "mean_expr_group1": 18.58,
      "mean_expr_group2": 47.40,
      "pct_group1": 0.998,
      "pct_group2": 1.0,
      "p_value": 0.0,
      "p_value_adj": 0.0,
      "significant": true
    },
    ...
  ]
}
```

### 2. DEG Within Dataset

**Endpoint:** `POST /api/v1/deg/within-dataset`

**Use Cases:**

- Compare different cell types
- Compare different clusters
- Compare conditions within same dataset

**Example Request:**

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/deg/within-dataset" \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id": "AD093044",
    "group1_filters": {"cell_type": "Microglia"},
    "group2_filters": {"cell_type": "Astrocyte"},
    "min_pct": 0.1,
    "logfc_threshold": 1.0,
    "p_value_threshold": 0.01,
    "top_n": 50
  }'
```

**Parameters:**

| Parameter           | Type   | Required | Default | Description                                              |
| ------------------- | ------ | -------- | ------- | -------------------------------------------------------- |
| `dataset_id`        | string | ✅       | -       | Dataset to analyze                                       |
| `group1_filters`    | object | ✅       | -       | Filters for group 1 (e.g., `{"cell_type": "Neuron"}`)    |
| `group2_filters`    | object | ✅       | -       | Filters for group 2 (e.g., `{"cell_type": "Astrocyte"}`) |
| `min_pct`           | float  | ❌       | 0.1     | Minimum fraction of cells expressing gene                |
| `logfc_threshold`   | float  | ❌       | 0.25    | Minimum log2 fold change                                 |
| `p_value_threshold` | float  | ❌       | 0.05    | Adjusted p-value threshold                               |
| `top_n`             | int    | ❌       | null    | Return only top N genes                                  |

### 3. Get Available Cell Types

**Endpoint:** `GET /api/v1/deg/cell-types/{dataset_id}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/deg/cell-types/AD093044"
```

**Response:**

```json
{
  "dataset_id": "AD093044",
  "cell_types": [
    "Astrocyte",
    "Microglia",
    "Neuron",
    "Oligodendrocyte",
    ...
  ],
  "counts": {
    "Microglia": 401,
    "Astrocyte": 119,
    ...
  }
}
```

### 4. Get Metadata Columns

**Endpoint:** `GET /api/v1/deg/metadata-columns/{dataset_id}`

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/deg/metadata-columns/AD093044"
```

**Response:**

```json
{
  "dataset_id": "AD093044",
  "metadata_columns": {
    "cell_type": {
      "unique_values": ["Microglia", "Astrocyte", ...],
      "n_unique": 20,
      "counts": {...}
    },
    "cluster_name": {
      "unique_values": ["Cluster1", "Cluster2", ...],
      "n_unique": 15,
      "counts": {...}
    }
  }
}
```

## Understanding the Results

### Gene Result Fields

| Field              | Description                                                                             |
| ------------------ | --------------------------------------------------------------------------------------- |
| `gene`             | Gene symbol                                                                             |
| `log2_fold_change` | Log2(Group1/Group2). Positive = upregulated in group1, Negative = upregulated in group2 |
| `mean_expr_group1` | Mean expression level in group 1                                                        |
| `mean_expr_group2` | Mean expression level in group 2                                                        |
| `pct_group1`       | Percentage of cells expressing gene in group 1 (0-1)                                    |
| `pct_group2`       | Percentage of cells expressing gene in group 2 (0-1)                                    |
| `p_value`          | Raw p-value from t-test                                                                 |
| `p_value_adj`      | Adjusted p-value (FDR-corrected)                                                        |
| `significant`      | Boolean indicating if gene passes thresholds                                            |

### Interpreting Log2 Fold Change

- `log2_fold_change > 0`: Gene is upregulated in group 1
- `log2_fold_change < 0`: Gene is upregulated in group 2
- `log2_fold_change = 1`: 2-fold increase
- `log2_fold_change = -1`: 2-fold decrease
- `log2_fold_change = 2`: 4-fold increase

## Example Use Cases

### 1. Find Microglia Markers Between Datasets

```javascript
const response = await fetch("/api/v1/deg/between-datasets", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    dataset_id1: "AD093044",
    dataset_id2: "AD093045",
    cell_type: "Microglia",
    min_pct: 0.1,
    logfc_threshold: 0.5,
    p_value_threshold: 0.05,
    top_n: 100,
  }),
});

const data = await response.json();

// data.summary.significant_genes = number of DEGs
// data.genes = list of DEGs sorted by significance
```

### 2. Compare Two Cell Types

```javascript
const response = await fetch("/api/v1/deg/within-dataset", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    dataset_id: "AD093044",
    group1_filters: { cell_type: "Microglia" },
    group2_filters: { cell_type: "Astrocyte" },
    min_pct: 0.2,
    logfc_threshold: 1.0,
    p_value_threshold: 0.01,
    top_n: 50,
  }),
});

const data = await response.json();

// Genes upregulated in Microglia (positive log2FC)
const microglia_markers = data.genes.filter((g) => g.log2_fold_change > 0);

// Genes upregulated in Astrocyte (negative log2FC)
const astrocyte_markers = data.genes.filter((g) => g.log2_fold_change < 0);
```

### 3. Interactive Cell Type Selection

```javascript
// 1. Get available cell types
const cellTypesResponse = await fetch("/api/v1/deg/cell-types/AD093044");
const cellTypes = await cellTypesResponse.json();

// 2. Let user select two cell types
const cellType1 = "Microglia"; // User selection
const cellType2 = "Astrocyte"; // User selection

// 3. Run DEG analysis
const degResponse = await fetch("/api/v1/deg/within-dataset", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    dataset_id: "AD093044",
    group1_filters: { cell_type: cellType1 },
    group2_filters: { cell_type: cellType2 },
    top_n: 100,
  }),
});

const results = await degResponse.json();

// 4. Display results
console.log(`Found ${results.summary.significant_genes} DEGs`);
console.log(
  `Upregulated in ${cellType1}: ${results.summary.upregulated_in_group1}`
);
console.log(
  `Upregulated in ${cellType2}: ${results.summary.upregulated_in_group2}`
);
```

### 4. Volcano Plot Data

```javascript
const response = await fetch("/api/v1/deg/between-datasets", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    dataset_id1: "AD093044",
    dataset_id2: "AD093045",
    cell_type: "Microglia",
  }),
});

const data = await response.json();

// Prepare data for volcano plot
const volcanoData = data.genes.map((gene) => ({
  x: gene.log2_fold_change,
  y: -Math.log10(gene.p_value_adj),
  label: gene.gene,
  significant: gene.significant,
}));

// Plot with Plotly
Plotly.newPlot(
  "volcano",
  [
    {
      x: volcanoData.map((d) => d.x),
      y: volcanoData.map((d) => d.y),
      mode: "markers",
      type: "scatter",
      text: volcanoData.map((d) => d.label),
      marker: {
        color: volcanoData.map((d) => (d.significant ? "red" : "gray")),
      },
    },
  ],
  {
    xaxis: { title: "Log2 Fold Change" },
    yaxis: { title: "-Log10(Adjusted P-value)" },
  }
);
```

## Performance & Caching

### First Run (No Cache)

- **Time**: 5-30 seconds (depending on dataset size and filtering)
- Results are computed and cached automatically

### Subsequent Runs (Cached)

- **Time**: <100ms
- Same parameters return cached results instantly
- Cache key is based on all parameters

### Cache Location

```
h5ad/deg_cache/{md5_hash}.json
```

### Clear Cache

```bash
# Clear all DEG cache
rm -rf h5ad/deg_cache/*

# Clear specific analysis
# (cache file names are MD5 hashes of parameters)
```

## Best Practices

### 1. Parameter Selection

**min_pct** (Minimum Percentage)

- `0.1` (10%): Standard, filters out very lowly expressed genes
- `0.2` (20%): More stringent, focuses on more widely expressed genes
- `0.05` (5%): More permissive, includes more rare genes

**logfc_threshold** (Log2 Fold Change)

- `0.25`: Minimal difference (1.19-fold)
- `0.5`: Moderate difference (1.41-fold)
- `1.0`: Strong difference (2-fold)
- `2.0`: Very strong difference (4-fold)

**p_value_threshold** (Adjusted P-value)

- `0.05`: Standard
- `0.01`: More stringent
- `0.001`: Very stringent

### 2. Cell Number Considerations

- **Minimum**: At least 20 cells per group recommended
- **Optimal**: 100+ cells per group for robust statistics
- **Warning**: Very small groups (<10 cells) may produce unreliable results

### 3. Result Interpretation

1. **Check summary first**:

   - `total_genes_tested`: How many genes passed min_pct filter
   - `significant_genes`: How many passed all thresholds

2. **Focus on top genes**:

   - Results are sorted by significance
   - Use `top_n` parameter to limit results

3. **Consider biological context**:
   - High fold change + low p-value = strong, confident result
   - High fold change + high p-value = variable gene
   - Low fold change + low p-value = consistent but small effect

## Troubleshooting

### "No cells remaining after filtering"

- Check cell type spelling (case-sensitive)
- Use `/deg/cell-types/{id}` to see available options
- Verify metadata column names with `/deg/metadata-columns/{id}`

### "Analysis taking too long"

- First run computes and caches (may take 10-30s)
- Subsequent runs with same parameters are instant
- Consider using `top_n` to limit results

### "Too few significant genes"

- Try lowering `logfc_threshold` (e.g., 0.25)
- Try relaxing `p_value_threshold` (e.g., 0.1)
- Try lowering `min_pct` (e.g., 0.05)
- Check if groups are too similar

### "Too many significant genes"

- Try increasing `logfc_threshold` (e.g., 1.0)
- Try more stringent `p_value_threshold` (e.g., 0.01)
- Try increasing `min_pct` (e.g., 0.2)
- Use `top_n` to focus on most significant

## API Documentation

Interactive API documentation available at:

```
http://localhost:9117/sskind-backend/docs
```

Look for the "DEG Analysis" section for interactive testing.

## Statistical Methods

### T-test

- Two-sample t-test between group means
- Handles sparse matrices efficiently
- NaN values replaced with p=1.0 (no difference)

### Multiple Testing Correction

- Method: Benjamini-Hochberg FDR (False Discovery Rate)
- Controls for multiple comparisons
- Adjusts p-values to reduce false positives

### Gene Filtering

1. **Expression filter**: Gene must be expressed in ≥ min_pct cells in either group
2. **Fold change filter**: |log2FC| ≥ logfc_threshold
3. **Statistical filter**: Adjusted p-value < p_value_threshold

## Future Enhancements

Planned features:

- [ ] Wilcoxon rank-sum test option
- [ ] Batch effect correction
- [ ] Pseudobulk analysis
- [ ] Gene set enrichment analysis (GSEA)
- [ ] Export results to CSV
- [ ] Visualization endpoints (volcano plot, MA plot)
