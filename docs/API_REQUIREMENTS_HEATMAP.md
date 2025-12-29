# Heatmap/Dotplot API Requirements

## Overview

This API endpoint generates heatmap or dotplot visualizations showing gene expression across cell types. The endpoint should aggregate expression data for multiple genes and return it in a format suitable for visualization, with optional scaling and clustering.

## Endpoint

**POST** `/api/v1/h5ad/{dataset_id}/heatmap`

## Request Format

```json
{
  "gene_list": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "plot_type": "heatmap",
  "scale_expression": false,
  "cluster_rows": false,
  "cluster_columns": false
}
```

### Parameters

| Parameter          | Type          | Required | Default   | Description                                                |
| ------------------ | ------------- | -------- | --------- | ---------------------------------------------------------- |
| `gene_list`        | array[string] | ✅       | -         | List of gene symbols to visualize                          |
| `plot_type`        | string        | ❌       | "heatmap" | Visualization type: `"heatmap"` or `"dotplot"`             |
| `scale_expression` | boolean       | ❌       | false     | Scale gene expression (z-score normalization per gene)     |
| `cluster_rows`     | boolean       | ❌       | false     | Cluster rows (genes) using hierarchical clustering         |
| `cluster_columns`  | boolean       | ❌       | false     | Cluster columns (cell types) using hierarchical clustering |

### Parameter Details

#### `scale_expression` (Z-score Normalization)

- When `true`: Normalize each gene's expression across cell types using z-score
- Formula: `z = (x - mean) / std`
- This makes expression levels comparable across genes with different baseline expression
- Applied per gene (row-wise normalization)

#### `cluster_rows` (Gene Clustering)

- When `true`: Reorder genes using hierarchical clustering based on expression similarity
- Use distance metric (e.g., Euclidean, correlation) and linkage method (e.g., average, complete)
- Return genes in clustered order (not original order)
- Useful for identifying gene co-expression patterns

#### `cluster_columns` (Cell Type Clustering)

- When `true`: Reorder cell types using hierarchical clustering based on expression similarity
- Use distance metric and linkage method
- Return cell types in clustered order (not original order)
- Useful for identifying similar cell type expression patterns

## Response Format

### For Heatmap

```json
{
  "data": [
    [2.5, 1.8, 0.3, 4.2, ...],  // Gene 1 expression across cell types
    [1.2, 3.5, 2.1, 0.8, ...],  // Gene 2 expression across cell types
    ...
  ],
  "genes": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "cell_types": ["Microglia", "Astrocyte", "Neuron", "Oligodendrocyte", ...],
  "genes_used": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "genes_not_found": [],
  "plot_type": "heatmap",
  "scale_expression": false,
  "cluster_rows": false,
  "cluster_columns": false,
  "stats": {
    "min_expression": 0.0,
    "max_expression": 5.2,
    "mean_expression": 2.1
  }
}
```

### For Dotplot

```json
{
  "data": [
    [
      {"expression": 2.5, "percentage": 85.3},
      {"expression": 1.8, "percentage": 72.1},
      {"expression": 0.3, "percentage": 15.2},
      ...
    ],
    [
      {"expression": 1.2, "percentage": 65.4},
      {"expression": 3.5, "percentage": 91.2},
      ...
    ],
    ...
  ],
  "genes": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "cell_types": ["Microglia", "Astrocyte", "Neuron", "Oligodendrocyte", ...],
  "genes_used": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
  "genes_not_found": [],
  "plot_type": "dotplot",
  "scale_expression": false,
  "cluster_rows": false,
  "cluster_columns": false,
  "stats": {
    "min_expression": 0.0,
    "max_expression": 5.2,
    "mean_expression": 2.1,
    "min_percentage": 0.0,
    "max_percentage": 100.0
  }
}
```

### Response Fields

| Field              | Type          | Description                                                                     |
| ------------------ | ------------- | ------------------------------------------------------------------------------- |
| `data`             | array[array]  | 2D array: `data[gene_index][cell_type_index]` (reordered if clustering applied) |
| `genes`            | array[string] | List of genes (in clustered order if `cluster_rows=true`)                       |
| `cell_types`       | array[string] | List of cell types (in clustered order if `cluster_columns=true`)               |
| `genes_used`       | array[string] | Genes found in the dataset                                                      |
| `genes_not_found`  | array[string] | Genes not found (if any)                                                        |
| `plot_type`        | string        | "heatmap" or "dotplot"                                                          |
| `scale_expression` | boolean       | Whether expression was scaled                                                   |
| `cluster_rows`     | boolean       | Whether rows were clustered                                                     |
| `cluster_columns`  | boolean       | Whether columns were clustered                                                  |
| `stats`            | object        | Summary statistics                                                              |

## Implementation Notes

### Data Aggregation

1. **For each gene and cell type combination:**

   - Filter cells by cell type
   - Calculate average expression (use normalized/log-transformed data)
   - For dotplot: also calculate percentage of cells with expression > 0

2. **Expression Calculation:**

   - Use normalized expression data (not raw counts)
   - Apply log1p transformation: `log(1 + expression)`
   - Handle zero values appropriately

3. **Scaling (`scale_expression=true`):**

   - For each gene (row), calculate z-score across all cell types
   - Formula: `z = (x - mean(x)) / std(x)`
   - This normalizes each gene independently
   - After scaling, expression values will be centered around 0 with std=1

4. **Clustering (`cluster_rows=true` or `cluster_columns=true`):**
   - Use hierarchical clustering (e.g., scipy.cluster.hierarchy)
   - Distance metric: Euclidean or correlation distance
   - Linkage method: average, complete, or ward
   - Reorder data, genes, and cell_types arrays according to clustering dendrogram
   - Return reordered arrays in response

### Example Python Implementation (Pseudocode)

```python
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist, squareform

def generate_heatmap(dataset_id, gene_list, plot_type="heatmap",
                     scale_expression=False, cluster_rows=False,
                     cluster_columns=False):
    # Load dataset
    adata = load_h5ad(dataset_id)

    # Get unique cell types
    cell_types = sorted(adata.obs['cell_type'].unique())

    # Prepare data matrix
    data = []
    genes_used = []
    genes_not_found = []

    for gene in gene_list:
        if gene not in adata.var_names:
            genes_not_found.append(gene)
            continue

        genes_used.append(gene)
        gene_data = []

        for cell_type in cell_types:
            # Filter cells by type
            cell_mask = adata.obs['cell_type'] == cell_type
            gene_expression = adata[cell_mask, gene].X

            # Calculate average (log-transformed)
            avg_expr = np.mean(np.log1p(gene_expression))

            if plot_type == "dotplot":
                # Calculate percentage expressing
                pct_expressing = np.mean(gene_expression > 0) * 100
                gene_data.append({
                    "expression": float(avg_expr),
                    "percentage": float(pct_expressing)
                })
            else:
                gene_data.append(float(avg_expr))

        data.append(gene_data)

    # Convert to numpy array for easier manipulation
    if plot_type == "heatmap":
        data_array = np.array(data)
    else:
        # Extract expression values for scaling/clustering
        expr_array = np.array([[d['expression'] for d in row] for row in data])

    # Apply scaling if requested
    if scale_expression:
        if plot_type == "heatmap":
            # Z-score normalize each row (gene)
            data_array = (data_array - data_array.mean(axis=1, keepdims=True)) / (data_array.std(axis=1, keepdims=True) + 1e-10)
        else:
            expr_array = (expr_array - expr_array.mean(axis=1, keepdims=True)) / (expr_array.std(axis=1, keepdims=True) + 1e-10)
            # Update expression values in data
            for i, row in enumerate(data):
                for j, point in enumerate(row):
                    point['expression'] = float(expr_array[i, j])

    # Apply clustering if requested
    if cluster_rows:
        if plot_type == "heatmap":
            # Cluster genes (rows)
            linkage_matrix = linkage(data_array, method='average', metric='euclidean')
            row_order = leaves_list(linkage_matrix)
            data = [data[i] for i in row_order]
            genes_used = [genes_used[i] for i in row_order]
        else:
            linkage_matrix = linkage(expr_array, method='average', metric='euclidean')
            row_order = leaves_list(linkage_matrix)
            data = [data[i] for i in row_order]
            genes_used = [genes_used[i] for i in row_order]

    if cluster_columns:
        if plot_type == "heatmap":
            # Cluster cell types (columns)
            data_array = np.array(data)
            linkage_matrix = linkage(data_array.T, method='average', metric='euclidean')
            col_order = leaves_list(linkage_matrix)
            data = [[row[i] for i in col_order] for row in data]
            cell_types = [cell_types[i] for i in col_order]
        else:
            expr_array = np.array([[d['expression'] for d in row] for row in data])
            linkage_matrix = linkage(expr_array.T, method='average', metric='euclidean')
            col_order = leaves_list(linkage_matrix)
            data = [[row[i] for i in col_order] for row in data]
            cell_types = [cell_types[i] for i in col_order]

    # Calculate stats
    if plot_type == "heatmap":
        all_values = [v for row in data for v in row]
    else:
        all_values = [v['expression'] for row in data for v in row]

    stats = {
        "min_expression": float(np.min(all_values)),
        "max_expression": float(np.max(all_values)),
        "mean_expression": float(np.mean(all_values))
    }

    if plot_type == "dotplot":
        all_pcts = [v['percentage'] for row in data for v in row]
        stats["min_percentage"] = float(np.min(all_pcts))
        stats["max_percentage"] = float(np.max(all_pcts))

    return {
        "data": data,
        "genes": genes_used,
        "cell_types": cell_types,
        "genes_used": genes_used,
        "genes_not_found": genes_not_found,
        "plot_type": plot_type,
        "scale_expression": scale_expression,
        "cluster_rows": cluster_rows,
        "cluster_columns": cluster_columns,
        "stats": stats
    }
```

## Example Requests

### Basic Heatmap

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/heatmap" \
  -H "Content-Type: application/json" \
  -d '{
    "gene_list": ["APOE", "APP", "MAPT"],
    "plot_type": "heatmap"
  }'
```

### Heatmap with Scaling and Clustering

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/heatmap" \
  -H "Content-Type: application/json" \
  -d '{
    "gene_list": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
    "plot_type": "heatmap",
    "scale_expression": true,
    "cluster_rows": true,
    "cluster_columns": true
  }'
```

### Dotplot with Clustering

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/AD093044/heatmap" \
  -H "Content-Type: application/json" \
  -d '{
    "gene_list": ["APOE", "APP", "MAPT"],
    "plot_type": "dotplot",
    "cluster_rows": true,
    "cluster_columns": false
  }'
```

## Notes

- **Scaling** is applied before clustering (if both are enabled)
- **Clustering** reorders the data, so the `genes` and `cell_types` arrays in the response reflect the clustered order
- Use consistent clustering parameters (distance metric, linkage method) for reproducibility
- For large gene lists (>100 genes), clustering may be slow - consider caching or optimization
