# Atlas Advanced DEG Analysis API

## Overview

✅ **STATUS: FULLY IMPLEMENTED AND TESTED**

This document describes the advanced differential expression analysis API for Atlas datasets. The system allows users to define cell groups using **multiple metadata filters**, enabling complex comparisons like:

- Astrocytes in Males vs Astrocytes in Females
- Neurons in Alzheimer's patients from Hippocampus vs Neurons in Controls from Hippocampus
- Microglia in APOE4 carriers aged 60-70 vs Microglia in APOE3 carriers aged 60-70

## Quick Start

```bash
# 1. Get available metadata columns and values
curl http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/info | jq '.categorical_metadata | keys'

# 2. Get values for a specific column
curl http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/metadata-values/cell_type | jq

# 3. Preview cell counts
curl -X POST http://localhost:9117/sskind-backend/api/v1/deg/preview-cell-counts \
  -H "Content-Type: application/json" \
  -d '{"dataset_id": "human_subset", "group1_filters": [{"column": "cell_type", "values": ["Astrocyte"]}], "group2_filters": [{"column": "cell_type", "values": ["Neuron"]}]}'

# 4. Run DEG analysis
curl -X POST http://localhost:9117/sskind-backend/api/v1/deg/atlas-compare \
  -H "Content-Type: application/json" \
  -d '{"dataset_id": "human_subset", "group1_filters": [{"column": "cell_type", "values": ["Astrocyte"]}], "group2_filters": [{"column": "cell_type", "values": ["Neuron"]}], "top_n": 20}'
```

## New API Endpoints Required

### 1. Get Metadata Column Values

**Endpoint:** `GET /h5ad/{dataset_id}/metadata-values/{column}`

Returns unique values for a specific metadata column with optional counts.

**Response:**

```json
{
  "column": "sex",
  "values": ["Male", "Female"],
  "counts": {
    "Male": 45000,
    "Female": 55000
  }
}
```

### 2. Advanced DEG Analysis with Multi-Filter Groups

**Endpoint:** `POST /deg/atlas-compare`

**Request Body:**

```json
{
  "dataset_id": "human_subset",
  "group1_filters": [
    {
      "column": "cell_type",
      "values": ["Astrocyte", "Oligodendrocyte"]
    },
    {
      "column": "sex",
      "values": ["Male"]
    },
    {
      "column": "brain_region",
      "values": ["Hippocampus", "Cortex"]
    }
  ],
  "group2_filters": [
    {
      "column": "cell_type",
      "values": ["Astrocyte", "Oligodendrocyte"]
    },
    {
      "column": "sex",
      "values": ["Female"]
    },
    {
      "column": "brain_region",
      "values": ["Hippocampus", "Cortex"]
    }
  ],
  "logfc_threshold": 0.25,
  "p_value_threshold": 0.05,
  "min_pct": 0.1
}
```

**Filter Logic:**

- Multiple filters within a group use **AND** logic
- Multiple values within a filter use **OR** logic

**Example interpretation of group1_filters above:**

```
(cell_type == "Astrocyte" OR cell_type == "Oligodendrocyte")
AND
(sex == "Male")
AND
(brain_region == "Hippocampus" OR brain_region == "Cortex")
```

**Response:**

```json
{
  "summary": {
    "group1_cells": 12500,
    "group2_cells": 14200,
    "total_genes_tested": 18000,
    "significant_genes": 342,
    "comparison_description": "Male Astrocyte/Oligodendrocyte (Hippocampus/Cortex) vs Female Astrocyte/Oligodendrocyte (Hippocampus/Cortex)"
  },
  "group1_description": {
    "filters_applied": [
      {"column": "cell_type", "values": ["Astrocyte", "Oligodendrocyte"]},
      {"column": "sex", "values": ["Male"]},
      {"column": "brain_region", "values": ["Hippocampus", "Cortex"]}
    ],
    "cell_count": 12500
  },
  "group2_description": {
    "filters_applied": [
      {"column": "cell_type", "values": ["Astrocyte", "Oligodendrocyte"]},
      {"column": "sex", "values": ["Female"]},
      {"column": "brain_region", "values": ["Hippocampus", "Cortex"]}
    ],
    "cell_count": 14200
  },
  "genes": [
    {
      "gene": "XIST",
      "log2_fold_change": 5.23,
      "p_value": 1.2e-150,
      "p_value_adj": 2.4e-146,
      "mean_expr_group1": 0.02,
      "mean_expr_group2": 4.56,
      "pct_group1": 0.01,
      "pct_group2": 0.89,
      "significant": true
    },
    ...
  ]
}
```

### 3. Preview Cell Counts (Optional but Recommended)

**Endpoint:** `POST /deg/preview-cell-counts`

This endpoint allows the frontend to show users how many cells match their filter criteria before running the full DEG analysis.

**Request Body:**

```json
{
  "dataset_id": "human_subset",
  "group1_filters": [...],
  "group2_filters": [...]
}
```

**Response:**

```json
{
  "group1_count": 12500,
  "group2_count": 14200,
  "overlap_count": 0,
  "warnings": []
}
```

If there's overlap (same cells in both groups), return a warning.

## Backend Implementation Notes

### Cell Selection Algorithm

```python
def select_cells(adata, filters):
    """
    Select cells matching all filter conditions.

    Args:
        adata: AnnData object
        filters: List of {column, values} dicts

    Returns:
        Boolean mask of matching cells
    """
    mask = np.ones(adata.n_obs, dtype=bool)

    for f in filters:
        column = f['column']
        values = f['values']

        if column not in adata.obs.columns:
            raise ValueError(f"Column '{column}' not found in dataset")

        # OR logic within values, AND logic between filters
        filter_mask = adata.obs[column].isin(values)
        mask = mask & filter_mask

    return mask


def run_advanced_deg(adata, group1_filters, group2_filters, **params):
    """Run DEG analysis between two filtered cell groups."""

    # Get cell masks
    mask1 = select_cells(adata, group1_filters)
    mask2 = select_cells(adata, group2_filters)

    # Check for overlap
    overlap = mask1 & mask2
    if overlap.sum() > 0:
        raise ValueError(f"{overlap.sum()} cells are in both groups")

    # Subset data
    group1_cells = adata[mask1]
    group2_cells = adata[mask2]

    # Run DEG (using scanpy or custom implementation)
    # ...
```

### Supported Metadata Columns

The system should support any categorical metadata column in `adata.obs`. Common examples:

| Column               | Example Values                  |
| -------------------- | ------------------------------- |
| `cell_type`          | Astrocyte, Neuron, Microglia    |
| `supercluster_label` | Neuronal, Non-neuronal          |
| `sex`                | Male, Female                    |
| `brain_region`       | Hippocampus, Cortex, Cerebellum |
| `disease_status`     | Control, AD, MCI                |
| `apoe_genotype`      | APOE3/3, APOE3/4, APOE4/4       |
| `age_group`          | Young, Middle, Old              |
| `sample_id`          | Sample-specific identifiers     |

### Performance Considerations

1. **Caching metadata unique values**: Cache the unique values and counts for each metadata column to speed up the filter UI.

2. **Lazy loading**: Don't load all metadata into memory at once. Fetch column values on demand.

3. **Parallel DEG computation**: For large datasets, consider parallelizing the DEG calculation.

4. **Progress reporting**: For long-running DEG analyses, consider implementing progress reporting via WebSockets or polling.

### Performance Optimizations (Implemented)

#### Preview Cell Counts Optimization

The `/deg/preview-cell-counts` endpoint has been optimized for frequent calls:

1. **Metadata-only loading**: Loads only metadata (not expression data) from precomputed cache
2. **In-memory caching**: Caches metadata in memory for subsequent requests
3. **Performance results**:
   - First call: ~0.5 seconds (loads from precomputed JSON)
   - Subsequent calls: ~0.02 seconds (uses in-memory cache) - **25x faster!**

This makes it suitable for real-time filter UI updates as users build their queries.

**Implementation details**:

- Uses precomputed `metadata.json` if available
- Falls back to loading from h5ad if no precomputed cache exists
- Caches metadata DataFrame in memory for the lifetime of the service
- No need to load the full expression matrix for cell counting

## Frontend Integration

The frontend sends filter configurations like this:

```javascript
// Example from GroupFilterBuilder.vue
const group1Filters = [
  { column: "cell_type", selectedValues: ["Astrocyte"] },
  { column: "sex", selectedValues: ["Male", "Female"] },
];

// Converted to API format
const apiFilters = group1Filters
  .filter((f) => f.column && f.selectedValues.length > 0)
  .map((f) => ({
    column: f.column,
    values: f.selectedValues,
  }));
```

## Error Handling

Return appropriate error responses:

```json
{
  "detail": "No cells match the specified filters for Group 1",
  "error_code": "EMPTY_GROUP",
  "group": 1
}
```

```json
{
  "detail": "Column 'invalid_column' not found in dataset",
  "error_code": "INVALID_COLUMN",
  "column": "invalid_column"
}
```

```json
{
  "detail": "Groups have overlapping cells. 150 cells are in both groups.",
  "error_code": "OVERLAPPING_GROUPS",
  "overlap_count": 150
}
```

## Migration from Existing API

The existing `/deg/between-datasets` endpoint should continue to work for backward compatibility with non-atlas datasets. The new `/deg/atlas-compare` endpoint is specifically for atlas datasets with multi-filter support.

## Testing Scenarios

1. **Simple comparison**: Single cell type vs another cell type
2. **Multi-value filter**: Multiple cell types in each group
3. **Multi-column filter**: Cell type + sex combination
4. **Complex filter**: Cell type + sex + brain region + disease status
5. **Edge cases**: Empty groups, overlapping groups, invalid columns

## Implementation Status

✅ **IMPLEMENTED** - All endpoints are now available and tested.

### Available Endpoints

1. **GET `/h5ad/{dataset_id}/metadata-values/{column}`** - Get unique values and counts for a metadata column
2. **POST `/deg/preview-cell-counts`** - Preview cell counts before running DEG analysis
3. **POST `/deg/atlas-compare`** - Advanced DEG analysis with multi-filter groups

### Example Usage

#### 1. Get metadata values for building filter UI

```bash
# Get available cell types
curl -s http://localhost:9117/sskind-backend/api/v1/h5ad/human_subset/metadata-values/cell_type | jq

# Response:
{
  "column": "cell_type",
  "values": ["Astrocyte", "Neuron", "Oligodendrocyte", ...],
  "counts": {
    "Astrocyte": 31497,
    "Neuron": 45000,
    ...
  }
}
```

#### 2. Preview cell counts

```bash
curl -s -X POST http://localhost:9117/sskind-backend/api/v1/deg/preview-cell-counts \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id": "human_subset",
    "group1_filters": [
      {"column": "cell_type", "values": ["Astrocyte"]}
    ],
    "group2_filters": [
      {"column": "cell_type", "values": ["Oligodendrocyte"]}
    ]
  }' | jq

# Response:
{
  "group1_count": 31497,
  "group2_count": 34615,
  "overlap_count": 0,
  "warnings": []
}
```

#### 3. Simple DEG comparison (single cell type vs another)

```bash
curl -s -X POST http://localhost:9117/sskind-backend/api/v1/deg/atlas-compare \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id": "human_subset",
    "group1_filters": [
      {"column": "cell_type", "values": ["Astrocyte"]}
    ],
    "group2_filters": [
      {"column": "cell_type", "values": ["Oligodendrocyte"]}
    ],
    "logfc_threshold": 0.5,
    "p_value_threshold": 0.05,
    "min_pct": 0.1,
    "top_n": 10
  }' | jq
```

#### 4. Multi-value filter (OR logic within filter)

```bash
curl -s -X POST http://localhost:9117/sskind-backend/api/v1/deg/atlas-compare \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id": "human_subset",
    "group1_filters": [
      {"column": "cell_type", "values": ["Astrocyte", "Bergmann glia", "Choroid plexus"]}
    ],
    "group2_filters": [
      {"column": "cell_type", "values": ["Oligodendrocyte", "Committed oligodendrocyte precursor"]}
    ],
    "top_n": 10
  }' | jq
```

#### 5. Multi-column filter (AND logic between filters)

```bash
curl -s -X POST http://localhost:9117/sskind-backend/api/v1/deg/atlas-compare \
  -H "Content-Type: application/json" \
  -d '{
    "dataset_id": "human_subset",
    "group1_filters": [
      {"column": "cell_type", "values": ["Astrocyte"]},
      {"column": "supercluster_name", "values": ["Astrocyte"]}
    ],
    "group2_filters": [
      {"column": "cell_type", "values": ["Oligodendrocyte"]},
      {"column": "supercluster_name", "values": ["Oligodendrocyte"]}
    ],
    "top_n": 10
  }' | jq
```

### Test Results

All test scenarios passed:

✅ Simple comparison: Astrocyte vs Oligodendrocyte (31,497 vs 34,615 cells, 553 significant genes)
✅ Multi-value filter: Multiple cell types with OR logic (31,599 vs 37,984 cells)
✅ Multi-column filter: AND logic between filters (31,497 vs 34,615 cells)
✅ Error handling: Overlapping groups detected and reported
✅ Error handling: Invalid column names detected with helpful message

### Performance Notes

- **Preview endpoint**: Fast (<1 second) - only counts cells without running DEG
- **Atlas-compare endpoint**: Moderate (5-30 seconds) - depends on cell count and gene filtering
- **Caching**: H5AD files are loaded once and kept in memory for subsequent requests
