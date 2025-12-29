# Module Score API Guide

## Overview

The Module Score API allows you to calculate aggregate gene expression scores across multiple genes using the Seurat-style scoring method implemented in Scanpy. This is useful for gene signatures, pathway activity, and cell state characterization.

## Endpoint

**POST** `/api/v1/h5ad/{dataset_id}/module-score`

## What is a Module Score?

A **module score** is calculated using `scanpy.tl.score_genes`, which implements the Seurat scoring approach [Satija et al., 2015]. The score is the **average expression of a set of genes after subtraction by the average expression of a reference set of genes**. The reference set is randomly sampled from the gene pool for each binned expression value.

**Key characteristics:**

- **Relative scoring**: Scores can be negative (lower than reference) or positive (higher than reference)
- **Normalized**: Automatically accounts for overall expression levels
- **Robust**: Uses reference genes matched by expression level to reduce bias
- **Reproducible**: Same approach as Seurat's `AddModuleScore`

**Use cases:**

- Gene signature scoring (e.g., immune response signature)
- Pathway activity scoring
- Cell state characterization
- Multi-gene biomarker evaluation

## Request Format

```json
{
  "gene_list": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "use_raw": false
}
```

### Parameters

| Parameter   | Type          | Required | Default | Description                                      |
| ----------- | ------------- | -------- | ------- | ------------------------------------------------ |
| `gene_list` | array[string] | ✅       | -       | List of gene symbols                             |
| `use_raw`   | boolean       | ❌       | false   | Use raw counts (false = normalized, recommended) |

### How It Works

The scoring method (`scanpy.tl.score_genes`):

1. **Bins cells** by their total expression (25 bins by default)
2. **Selects reference genes** randomly from the gene pool (50 genes by default) for each bin
3. **Calculates score** = mean(expression of gene_list) - mean(expression of reference set)
4. **Normalizes** by expression level, so scores are comparable across cells

**Important notes:**

- Scores can be **negative** (gene set expressed less than reference)
- Scores can be **positive** (gene set expressed more than reference)
- Higher scores indicate stronger signature expression relative to background

## Response Format

```json
{
  "module_score": [0.6, -0.2, 1.4, 2.0, -0.5, ...],
  "cell_ids": ["cell_1", "cell_2", ...],
  "genes_used": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "genes_not_found": [],
  "method": "scanpy_score_genes",
  "n_genes_used": 5,
  "stats": {
    "mean": 0.1234,
    "median": 0.0500,
    "std": 0.5678,
    "min": -2.5,
    "max": 3.2,
    "cells_with_expression": 7436,
    "mean_expressed": 0.2345
  },
  "data_type": "normalized"
}
```

**Note**: Module scores can be negative because they represent expression relative to a reference set. A negative score means the gene set is expressed less than the reference genes for that cell's expression level.

### Response Fields

| Field             | Description                                                                                                |
| ----------------- | ---------------------------------------------------------------------------------------------------------- |
| `module_score`    | Array of scores (one per cell). Can be negative (lower than reference) or positive (higher than reference) |
| `cell_ids`        | Cell identifiers                                                                                           |
| `genes_used`      | Genes found in the dataset                                                                                 |
| `genes_not_found` | Genes not found (if any)                                                                                   |
| `method`          | Always "scanpy_score_genes" (Seurat-style scoring)                                                         |
| `n_genes_used`    | Number of genes used                                                                                       |
| `stats`           | Summary statistics                                                                                         |
| `data_type`       | "normalized" or "raw"                                                                                      |

## Example Use Cases

### 1. Alzheimer's Disease Signature

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/module-score" \
  -H "Content-Type: application/json" \
  -d '{
    "gene_list": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"]
  }'
```

### 2. T Cell Signature

```javascript
const response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"],
  }),
});

const data = await response.json();

console.log(`T cell score calculated for ${data.cell_ids.length} cells`);
console.log(`Mean score: ${data.stats.mean}`);
console.log(
  `Found ${data.n_genes_used}/${
    data.genes_used.length + data.genes_not_found.length
  } genes`
);
```

### 3. Microglia Activation Signature

```javascript
const response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: [
      "AIF1",
      "CD68",
      "ITGAM",
      "PTPRC", // Microglia markers
      "IL1B",
      "TNF",
      "IL6", // Inflammatory
      "CD74",
      "HLA-DRA", // Antigen presentation
    ],
  }),
});

const data = await response.json();
```

### 4. Cell Cycle Signature

```javascript
const g1s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4"];
const g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2"];

// Calculate G1/S score
const g1s_response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: g1s_genes,
  }),
});

// Calculate G2/M score
const g2m_response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: g2m_genes,
  }),
});

const g1s_data = await g1s_response.json();
const g2m_data = await g2m_response.json();

// Compare scores
for (let i = 0; i < g1s_data.cell_ids.length; i++) {
  if (
    g1s_data.module_score[i] > threshold ||
    g2m_data.module_score[i] > threshold
  ) {
    console.log(`Cell ${i} is cycling`);
  }
}
```

## Visualization Examples

### 1. UMAP with Module Score Overlay

```javascript
// Get UMAP coordinates and module score
const umap_response = await fetch("/api/v1/h5ad/AD093044/embedding/umap");
const umap = await umap_response.json();

const score_response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: ["APOE", "APP", "MAPT"],
  }),
});
const scores = await score_response.json();

// Plot with Plotly
Plotly.newPlot(
  "umap-plot",
  [
    {
      x: umap.coordinates.map((c) => c[0]),
      y: umap.coordinates.map((c) => c[1]),
      mode: "markers",
      type: "scattergl",
      marker: {
        color: scores.module_score,
        colorscale: "Viridis",
        showscale: true,
        colorbar: {
          title: "AD Signature Score",
        },
      },
      text: scores.cell_ids,
    },
  ],
  {
    title: "Alzheimer's Disease Signature Score",
    xaxis: { title: "UMAP 1" },
    yaxis: { title: "UMAP 2" },
  }
);
```

### 2. Module Score Distribution

```javascript
const response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  }),
});

const data = await response.json();

// Plot histogram
Plotly.newPlot(
  "histogram",
  [
    {
      x: data.module_score,
      type: "histogram",
      nbinsx: 50,
      marker: { color: "steelblue" },
    },
  ],
  {
    title: "Distribution of AD Signature Scores",
    xaxis: { title: "Module Score" },
    yaxis: { title: "Number of Cells" },
  }
);
```

### 3. Compare Module Scores Between Cell Types

```javascript
// Get module score
const score_response = await fetch("/api/v1/h5ad/AD093044/module-score", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    gene_list: ["APOE", "APP", "MAPT"],
  }),
});

// Get cell type metadata
const metadata_response = await fetch(
  "/api/v1/h5ad/AD093044/metadata?columns=cell_type"
);

const scores = await score_response.json();
const metadata = await metadata_response.json();

// Group by cell type
const cell_types = {};
for (let i = 0; i < scores.cell_ids.length; i++) {
  const cell_type = metadata.data.cell_type[i];
  if (!cell_types[cell_type]) {
    cell_types[cell_type] = [];
  }
  cell_types[cell_type].push(scores.module_score[i]);
}

// Plot box plot
const traces = Object.keys(cell_types).map((cell_type) => ({
  y: cell_types[cell_type],
  type: "box",
  name: cell_type,
}));

Plotly.newPlot("boxplot", traces, {
  title: "AD Signature Score by Cell Type",
  yaxis: { title: "Module Score" },
});
```

## Best Practices

### 1. Gene Selection

**Good practices:**

- Use 5-20 genes for robust signatures
- Include well-validated markers
- Use genes with similar biological function
- Check all genes are found: `genes_not_found` should be empty

**Example signatures:**

```javascript
// Well-curated signatures
const signatures = {
  microglia: ["AIF1", "CX3CR1", "P2RY12", "TMEM119", "HEXB"],
  astrocyte: ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "GJA1"],
  oligodendrocyte: ["MBP", "MOG", "PLP1", "MAG", "MOBP"],
  neuron: ["SNAP25", "SYT1", "SYN1", "RBFOX3", "NEFL"],
};
```

### 2. Normalization

**Always use normalized data (`use_raw: false`)** for:

- ✅ Visualization
- ✅ Cross-sample comparison
- ✅ Module scores
- ✅ Most analyses

**Use raw data (`use_raw: true`)** only for:

- Statistical testing (DEG analysis handles this internally)
- Specific algorithms requiring raw counts

### 3. Interpretation

**Module score values:**

- Scores are **relative to a reference set** of genes
- **Negative scores**: Gene set expressed less than reference (low signature activity)
- **Positive scores**: Gene set expressed more than reference (high signature activity)
- **Zero**: Gene set expressed similarly to reference
- Compare relative scores between cells/groups
- Use percentiles or absolute thresholds for classification

**Example thresholding:**

```javascript
// Calculate percentiles
const sorted = [...data.module_score].sort((a, b) => a - b);
const p75 = sorted[Math.floor(sorted.length * 0.75)];
const p90 = sorted[Math.floor(sorted.length * 0.9)];

// Classify cells (scores can be negative!)
const classifications = data.module_score.map((score) => {
  if (score >= p90) return "high";
  if (score >= p75) return "medium";
  if (score > 0) return "low";
  return "negative"; // Below reference level
});
```

## Troubleshooting

### "None of the genes found"

- Check gene symbols (case-sensitive)
- Use `/h5ad/{id}/genes/search` to find correct names
- Verify gene names match the organism (human vs mouse)

### "Module scores all zero or very low"

- Check if genes are actually expressed
- Try different genes
- Verify data normalization
- Remember: scores are relative to reference, so low scores don't necessarily mean no expression

### Unexpected score distribution

- Check `genes_not_found` - missing genes affect scores
- Scores can be negative (this is normal - means below reference level)
- Compare `data_type` (normalized vs raw)
- Ensure you're using appropriate gene signatures for your cell types

## Performance

- **First calculation**: 2-10 seconds (depending on dataset size)
- **Response size**: ~100KB for 10,000 cells
- **No caching**: Module scores computed on-demand (fast with normalized data)

## API Documentation

Interactive documentation:

```
http://localhost:9117/sskind-backend/docs
```

Look for `POST /h5ad/{dataset_id}/module-score` under "H5AD Data" section.

## Common Gene Signatures

### Neuroscience

```javascript
const signatures = {
  // Alzheimer's disease
  ad_risk: ["APOE", "APP", "MAPT", "PSEN1", "PSEN2", "TREM2"],

  // Neuroinflammation
  neuroinflam: ["IL1B", "TNF", "IL6", "PTGS2", "NOS2"],

  // Synaptic function
  synapse: ["SYN1", "SYT1", "SNAP25", "VAMP2", "STX1A"],

  // Myelin
  myelin: ["MBP", "MOG", "PLP1", "MAG", "MOBP", "CNP"],

  // Oxidative stress
  oxidative: ["SOD1", "SOD2", "CAT", "GPX1", "NQO1"],
};
```

### Cell Type Markers

```javascript
const cell_markers = {
  microglia: ["AIF1", "CX3CR1", "P2RY12", "TMEM119", "HEXB", "CSF1R"],
  astrocyte: ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "GJA1", "ALDH1L1"],
  oligodendrocyte: ["MBP", "MOG", "PLP1", "MAG", "MOBP", "OLIG1"],
  neuron_excitatory: ["SLC17A7", "CAMK2A", "SATB2", "TBR1"],
  neuron_inhibitory: ["GAD1", "GAD2", "SLC32A1", "PVALB"],
  endothelial: ["CLDN5", "FLT1", "VWF", "PECAM1"],
  opc: ["PDGFRA", "CSPG4", "SOX10", "OLIG2"],
};
```

## Advanced: Custom Scoring

For custom scoring algorithms, you can:

1. Get individual gene expressions
2. Implement custom logic in frontend
3. Or extend the backend with new methods

```javascript
// Example: Custom z-score normalization
const gene_responses = await Promise.all(
  [
    fetch("/api/v1/h5ad/AD093044/expression/APOE"),
    fetch("/api/v1/h5ad/AD093044/expression/APP"),
    fetch("/api/v1/h5ad/AD093044/expression/MAPT"),
  ].map((p) => p.then((r) => r.json()))
);

// Z-score normalize each gene
const z_scores = gene_responses.map((gene) => {
  const expr = gene.expression;
  const mean = gene.stats.mean;
  const std = gene.stats.std;
  return expr.map((val) => (val - mean) / std);
});

// Average z-scores
const custom_score = z_scores[0].map(
  (_, i) =>
    z_scores.reduce((sum, scores) => sum + scores[i], 0) / z_scores.length
);
```
