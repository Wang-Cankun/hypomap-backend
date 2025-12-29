# Spatial Transcriptomics API Guide

## Overview

The ssKIND backend supports comprehensive spatial transcriptomics data visualization and analysis through RESTful API endpoints. This guide covers all spatial-specific features including spatial coordinates, Spatially Variable Genes (SVG), precomputed DEG, deconvolution predictions, and Cell-Cell Communication (CCC) data.

## API Discovery

The backend provides a queryable API metadata endpoint that returns structured information about all available endpoints. This allows frontends to programmatically discover available endpoints, their parameters, response schemas, and examples.

### Get API Metadata

```bash
GET /api/v1/api-metadata
```

Returns comprehensive metadata about all API endpoints organized by category.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/api-metadata"
```

**Response Structure:**

```json
{
  "base_url": "http://localhost:9117/sskind-backend/api/v1",
  "api_version": "1.0.0",
  "categories": {
    "spatial": [...],
    "general_h5ad": [...],
    "analysis_features": [...],
    "deg": [...]
  },
  "quick_reference": {
    "check_features": "/h5ad/{dataset_id}/spatial/info",
    "spatial_coords": "/h5ad/{dataset_id}/spatial/coordinates",
    "spatial_plot": "/h5ad/{dataset_id}/spatial/plot-data?genes=GENE1,GENE2",
    "list_datasets": "/h5ad/datasets",
    "list_features": "/h5ad/analysis-features"
  }
}
```

Each endpoint in the categories includes:

- `name`: Human-readable name
- `method`: HTTP method (GET, POST, etc.)
- `path`: Endpoint path with placeholders
- `description`: What the endpoint does
- `path_params`: Path parameters with types and descriptions
- `query_params`: Query parameters with types, defaults, and descriptions
- `response`: Response type and schema
- `example`: Example URL and response

**Frontend Usage Example:**

```javascript
// Fetch API metadata
const metadata = await fetch("/api/v1/api-metadata").then((r) => r.json());

// Get spatial endpoints
const spatialEndpoints = metadata.categories.spatial;

// Find endpoint by name
const spatialInfoEndpoint = spatialEndpoints.find(
  (e) => e.name === "Get Spatial Info"
);

// Build URL dynamically
const datasetId = "ST024001";
const url =
  metadata.base_url +
  spatialInfoEndpoint.path.replace("{dataset_id}", datasetId);

// Make request
const data = await fetch(url).then((r) => r.json());
```

## Supported Dataset Types

- **Visium (10x Genomics)**: Multi-cell-level spatial data
  - Features: Spatial coordinates, SVG, DEG, Deconvolution (Tangram)
  - **Demo Dataset ID**: `ST024001`
- **Xenium (10x Genomics)**: Single-cell-level spatial data
  - Features: Spatial coordinates, SVG, DEG, CCC (COMMOT)
  - **Demo Dataset ID**: `ST034001`

### Demo Datasets

The following datasets are available for testing:

| Dataset ID | Type   | Cells   | Genes | Features                                                            |
| ---------- | ------ | ------- | ----- | ------------------------------------------------------------------- |
| `ST024001` | Visium | 2,615   | 2,306 | Spatial coords, SVG, DEG (10 groups), Deconvolution (19 cell types) |
| `ST034001` | Xenium | 164,081 | 415   | Spatial coords, SVG, DEG (254 groups), CCC                          |

## Quick Start

### 1. Preprocess Spatial Dataset

```bash
# Preprocess Visium dataset (ST024001)
python scripts/preprocess_h5ad.py --h5ad h5ad/raw/ST024001.h5ad --dataset-id ST024001

# Preprocess Xenium dataset (ST034001)
python preprocess_h5ad.py --h5ad h5ad/raw/ST034001.h5ad --dataset-id ST034001
```

**Note**: The demo datasets `ST024001` (Visium) and `ST034001` (Xenium) are already preprocessed and available for testing.

### 2. Check Available Features

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/info"
```

**Response:**

```json
{
  "has_spatial_coordinates": true,
  "has_svg": true,
  "has_precomputed_deg": true,
  "has_deconvolution": true,
  "has_ccc": false,
  "dataset_type": "visium"
}
```

### 3. Get Spatial Coordinates

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/coordinates"
```

## API Endpoints

### Spatial Information

#### Get Spatial Features Info

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/info
```

Returns information about available spatial features in the dataset.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/info"
```

**Response:**

```json
{
  "has_spatial_coordinates": true,
  "has_svg": true,
  "has_precomputed_deg": true,
  "has_deconvolution": true,
  "has_ccc": false,
  "dataset_type": "visium"
}
```

### Spatial Coordinates

#### Get Spatial Coordinates

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/coordinates
```

Returns x, y spatial coordinates for all spots/cells.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/coordinates"
```

**Response:**

```json
{
  "coordinates": [
    [9575.0, 4911.0],
    [12671.0, 9443.0],
    ...
  ],
  "cell_ids": ["cell_1", "cell_2", ...],
  "spatial_key": "spatial"
}
```

#### Get Spatial Plot Data (Combined) - Full Resolution

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/plot-data?genes=APOE,APP&metadata=annotation,seurat_clusters
```

Returns combined spatial coordinates + gene expression + metadata in one request.

**Important:** This endpoint returns **FULL-RESOLUTION coordinates** (not transformed for image overlay).

- Use this when you **don't need** to overlay on tissue images
- Coordinates are in original full-resolution space (e.g., [9575.0, 4911.0])
- Does **NOT** accept `sample_key` or `image_key` parameters
- For image overlay, use `/spatial/plot-data-with-image` instead

**Parameters:**

- `genes` (optional): Comma-separated list of gene symbols
- `metadata` (optional): Comma-separated list of metadata columns (e.g., `annotation`, `seurat_clusters`, `cell_type`)

**Common Metadata Columns:**

- `annotation` - Cell type annotations
- `seurat_clusters` - Seurat cluster assignments
- `leiden` - Leiden cluster assignments
- `cell_type` - Cell type labels
- Any other column from `adata.obs`

**Example - Get cell types/clusters:**

```bash
# Get spatial coordinates with cell type annotations
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data?metadata=annotation,seurat_clusters"
```

**Response:**

```json
{
  "coordinates": [[9575.0, 4911.0], ...],
  "cell_ids": ["cell_1", ...],
  "metadata": {
    "annotation": ["ctx-deep-layers", "ctx-olfactory", ...],
    "seurat_clusters": ["0", "1", "2", ...]
  }
}
```

**Example - Get coordinates + genes + cell types:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data?genes=APOE&metadata=annotation"
```

**Response:**

```json
{
  "coordinates": [[9575.0, 4911.0], ...],
  "cell_ids": ["cell_1", ...],
  "genes": {
    "APOE": [0.5, 1.2, 0.0, ...]
  },
  "metadata": {
    "annotation": ["ctx-deep-layers", "ctx-olfactory", ...]
  }
}
```

### Spatial Images

#### Get Spatial Image Info

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/image/info?sample_key=slice1
```

Returns metadata about available spatial images including sample keys, image dimensions, and scale factors.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/image/info"
```

**Response:**

```json
{
  "samples": {
    "slice1": {
      "image_keys": ["hires"],
      "image_shapes": { "hires": [600, 600, 3] },
      "scalefactors": {
        "spot_diameter_fullres": 1.0,
        "tissue_hires_scalef": 0.04010695
      },
      "image_paths": { "hires": "spatial_images/slice1/hires.png" }
    }
  },
  "default_sample": "slice1"
}
```

#### Get Spatial Image File

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/image/{sample_key}/{image_key}
```

Returns the PNG image file for the specified sample and image type. The image can be used as a background for overlaying spatial coordinates.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/image/slice1/hires" -o tissue_image.png
```

**Response:** PNG image file (Content-Type: image/png)

#### Get Transformed Spatial Coordinates

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/coordinates/transformed?sample_key=slice1&image_key=hires
```

Returns spatial coordinates transformed to match image space. This is essential for overlaying coordinates on the tissue image.

**The transformation applies:** `image_coord = fullres_coord * tissue_hires_scalef`

**Parameters:**

- `sample_key` (optional): Sample/slice key (uses default if not specified)
- `image_key` (optional): Image type to transform for (default: "hires")

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/coordinates/transformed?image_key=hires"
```

**Response:**

```json
{
  "coordinates": [[384.02, 196.97], ...],
  "cell_ids": ["cell_1", ...],
  "sample_key": "slice1",
  "image_key": "hires",
  "scale_factor": 0.04010695,
  "original_coordinates": [[9575.0, 4911.0], ...],
  "image_shape": [600, 600, 3]
}
```

#### Get Spatial Plot Data with Image (Recommended for Image Overlay)

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/plot-data-with-image?genes=APOE&image_key=hires&sample_key=slice1
```

**This is the recommended endpoint for spatial visualization with image background.**

**Key differences from `/spatial/plot-data`:**

- Returns **TRANSFORMED coordinates** matching image space (e.g., [384.02, 196.97])
- Includes image URL and transformation metadata
- Accepts `sample_key` and `image_key` parameters to select specific image
- Coordinates are scaled to match the image dimensions

Returns everything needed for spatial visualization with image overlay:

- Transformed coordinates (matching image space)
- Image URL and metadata
- Gene expression (if genes specified)
- Metadata (if metadata columns specified)
- Scale factors for coordinate transformation

**Parameters:**

- `sample_key` (optional): Sample/slice key (e.g., "slice1", "slice1.1"). Uses default if not specified
- `image_key` (optional): Image type (default: "hires", can be "lowres" or other available types)
- `genes` (optional): Comma-separated list of genes
- `metadata` (optional): Comma-separated list of metadata columns (e.g., "annotation", "seurat_clusters")

**Note:** The `sample_key` and `image_key` parameters are **only used by this endpoint** to transform coordinates to match the selected image. The `/spatial/plot-data` endpoint does NOT accept these parameters.

**Design guidance:** If you want a single-purpose call for overlay (coords + metadata + image) without gene expression, use `/spatial/plot-coordinates-with-image`. Fetch gene expression separately (per existing expression endpoints) and join on `cell_ids`. This keeps calls small and focused.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data-with-image?genes=APOE&image_key=hires"
```

**Response:**

```json
{
  "coordinates": [[384.02, 196.97], ...],
  "cell_ids": ["cell_1", ...],
  "sample_key": "slice1",
  "image_key": "hires",
  "image_url": "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/image/slice1/hires",
  "image_shape": [600, 600, 3],
  "scale_factor": 0.04010695,
  "original_coordinates": [[9575.0, 4911.0], ...],
  "scalefactors": {
    "spot_diameter_fullres": 1.0,
    "tissue_hires_scalef": 0.04010695
  },
  "genes": {
    "APOE": [0.5, 1.2, ...]
  }
}
```

**Displaying Categories (Cell Type/Cluster):**

To display cell types, clusters, or other categorical metadata for each spot, use the `metadata` parameter:

```bash
# Get spatial data with cell type/cluster annotations (with image overlay)
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data-with-image?metadata=annotation,seurat_clusters&sample_key=slice1.1"

# Get spatial data with cell type/cluster annotations (without image, full resolution)
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data?metadata=annotation,seurat_clusters"
```

**Key Difference:**

- `/spatial/plot-data-with-image`: Coordinates are **transformed** to match image space (use for image overlay)
- `/spatial/plot-data`: Coordinates are **full resolution** (use when no image overlay needed)

#### Get Spatial Plot Coordinates with Image (No Gene Expression)

```bash
GET /api/v1/h5ad/{dataset_id}/spatial/plot-coordinates-with-image?metadata=annotation,seurat_clusters&sample_key=slice1&image_key=hires
```

Single-purpose endpoint for image overlay that returns:

- Transformed coordinates (image space)
- Metadata (if requested)
- Image URL + shape + scale factors
- Original coordinates

Use this when you do NOT want gene expression in the same call. Fetch gene expression separately (per gene or batch endpoints) and join on `cell_ids`.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-coordinates-with-image?metadata=annotation,seurat_clusters&sample_key=slice1&image_key=hires"
```

**Response includes metadata:**

```json
{
  "coordinates": [[384.02, 196.97], ...],
  "cell_ids": ["cell_1", ...],
  "metadata": {
    "annotation": ["ctx-deep-layers", "ctx-olfactory", ...],
    "seurat_clusters": ["0", "1", "2", ...]
  },
  "image_url": "...",
  ...
}
```

**Frontend Usage:**

```javascript
// Get plot data with image and cell type/cluster metadata
const response = await fetch(
  "/api/v1/h5ad/ST024001/spatial/plot-data-with-image?metadata=annotation,seurat_clusters&image_key=hires"
);
const data = await response.json();

// Load image
const img = new Image();
img.src = data.image_url;

// When image loads, overlay coordinates
img.onload = () => {
  const canvas = document.getElementById("spatial-canvas");
  const ctx = canvas.getContext("2d");

  // Set canvas size to match image
  canvas.width = data.image_shape[1]; // width
  canvas.height = data.image_shape[0]; // height

  // Draw image
  ctx.drawImage(img, 0, 0);

  // Overlay coordinates (already transformed to image space)
  data.coordinates.forEach((coord, i) => {
    const x = coord[0];
    const y = coord[1];

    // Option 1: Color by gene expression
    if (data.genes && data.genes.APOE) {
      const expression = data.genes.APOE[i];
      const color = getColorForExpression(expression);
      ctx.fillStyle = color;
      ctx.fillRect(x - 2, y - 2, 4, 4);
    }

    // Option 2: Color by cell type/cluster (if metadata provided)
    if (data.metadata && data.metadata.annotation) {
      const cellType = data.metadata.annotation[i];
      const color = getColorForCellType(cellType);
      ctx.fillStyle = color;
      ctx.fillRect(x - 2, y - 2, 4, 4);
    }
  });
};
```

**Example: Display Cell Types/Clusters on Spatial Plot**

```javascript
// Get spatial data with cell type annotations
const response = await fetch(
  "/api/v1/h5ad/ST024001/spatial/plot-data-with-image?metadata=annotation"
);
const data = await response.json();

// Load image
const img = new Image();
img.src = data.image_url;

img.onload = () => {
  const canvas = document.getElementById("spatial-canvas");
  const ctx = canvas.getContext("2d");
  canvas.width = data.image_shape[1];
  canvas.height = data.image_shape[0];

  // Draw image background
  ctx.drawImage(img, 0, 0);

  // Get unique cell types for color mapping
  const cellTypes = [...new Set(data.metadata.annotation)];
  const colorMap = createColorMap(cellTypes);

  // Overlay coordinates colored by cell type
  data.coordinates.forEach((coord, i) => {
    const x = coord[0];
    const y = coord[1];
    const cellType = data.metadata.annotation[i];
    const color = colorMap[cellType];

    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x, y, 3, 0, 2 * Math.PI);
    ctx.fill();
  });
};
```

### Spatially Variable Genes (SVG)

#### List All SVG

```bash
GET /api/v1/h5ad/{dataset_id}/svg/list?top_n=50&min_score=2.0
```

Returns all spatially variable genes with scores, ranks, and statistics.

**Parameters:**

- `top_n` (optional): Return only top N genes by rank
- `min_score` (optional): Minimum gft_score threshold

**Getting ALL SVG genes:**

- Omit `top_n` to return the full list (payload can be large).
- For ST024001 there are ~2306 SVG entries.
- The `/svg/top` endpoint caps `n` at 1000; use `/svg/list` when you need everything.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/svg/list?top_n=20"
```

**Response:**

```json
{
  "genes": [
    {
      "gene": "gm42418",
      "gft_score": 4.23,
      "svg_rank": 2.0,
      "cutoff_gft_score": 1.0,
      "pvalue": 1.33e-77,
      "fdr": 3.27e-76
    },
    ...
  ],
  "total_genes": 2306,
  "filtered_count": 20,
  "available_columns": ["gft_score", "svg_rank", "cutoff_gft_score", "pvalue", "fdr"]
}
```

#### Get Top N SVG

```bash
GET /api/v1/h5ad/{dataset_id}/svg/top?n=50
```

Convenience endpoint to get top N spatially variable genes.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/svg/top?n=10"
```

#### Get SVG for Specific Gene

```bash
GET /api/v1/h5ad/{dataset_id}/svg/{gene}
```

Returns SVG statistics for a single gene.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/svg/nrgn"
```

**Response:**

```json
{
  "gene": "nrgn",
  "gft_score": 3.91,
  "svg_rank": 3.0,
  "cutoff_gft_score": 1.0,
  "pvalue": 5.95e-94,
  "fdr": 2.03e-92
}
```

### Precomputed DEG

#### List DEG Groups

```bash
GET /api/v1/h5ad/{dataset_id}/deg/precomputed/groups
```

Returns list of available comparison groups in precomputed DEG results.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/deg/precomputed/groups"
```

**Response:**

```json
{
  "groups": ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"],
  "method": "wilcoxon"
}
```

#### Get Precomputed DEG Results

```bash
GET /api/v1/h5ad/{dataset_id}/deg/precomputed?group=0
```

Returns differentially expressed genes from precomputed analysis stored in the h5ad file.

**Parameters:**

- `group` (optional): Filter by specific group name

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/deg/precomputed?group=0"
```

**Response:**

```json
{
  "groups": ["0", "1", ...],
  "params": {
    "method": "wilcoxon",
    "use_raw": false
  },
  "genes": {
    "0": [
      {
        "gene": "gene1",
        "score": 0.85,
        "logfoldchange": 1.23,
        "pval": 0.001,
        "pval_adj": 0.005
      },
      ...
    ]
  },
  "method": "wilcoxon"
}
```

### Deconvolution (Visium Only)

#### Get Deconvolution Predictions

```bash
GET /api/v1/h5ad/{dataset_id}/deconvolution/predictions
```

Returns cell type predictions from Tangram deconvolution analysis.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/deconvolution/predictions"
```

**Response:**

```json
{
  "predictions": {
    "01 IT-ET Glut": [1.37, 1.37, 1.37, ...],
    "02 IT Glut": [0.85, 0.82, 0.88, ...],
    ...
  },
  "cell_ids": ["cell_1", "cell_2", ...],
  "cell_types": ["01 IT-ET Glut", "02 IT Glut", ...]
}
```

#### Get Deconvolution Plot Data

```bash
GET /api/v1/h5ad/{dataset_id}/deconvolution/plot-data
```

Returns spatial coordinates + deconvolution predictions for visualization.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/deconvolution/plot-data"
```

**Response:**

```json
{
  "coordinates": [[9575.0, 4911.0], ...],
  "cell_ids": ["cell_1", ...],
  "metadata": {
    "01 IT-ET Glut": [1.37, ...],
    "02 IT Glut": [0.85, ...],
    ...
  }
}
```

### Cell-Cell Communication (Xenium Only)

#### Get CCC Interactions

```bash
GET /api/v1/h5ad/{dataset_id}/ccc/interactions
```

Returns cell-cell communication interaction data from COMMOT analysis.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST034001/ccc/interactions"
```

**Response:**

```json
{
  "interactions": [...],
  "ligands": ["TAC1", ...],
  "receptors": ["TACR1", ...],
  "cell_types": [...],
  "data": {
    "info": {
      "df_ligrec": [
        {
          "ligand": "TAC1",
          "receptor": "TACR1",
          "pathway": "TAC"
        }
      ],
      "distance_threshold": 200
    },
    ...
  }
}
```

#### Get CCC Network

```bash
GET /api/v1/h5ad/{dataset_id}/ccc/network
```

Returns simplified network representation of cell-cell communication.

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST034001/ccc/network"
```

## Analysis Features Management

### List All Analysis Features

```bash
GET /api/v1/h5ad/analysis-features?dataset_type=visium&has_svg=true
```

Returns list of all datasets with their available analysis features.

**Query Parameters:**

- `dataset_type` (optional): Filter by `visium` or `xenium`
- `has_spatial` (optional): Filter by spatial coordinates availability
- `has_svg` (optional): Filter by SVG availability
- `has_deconvolution` (optional): Filter by deconvolution availability
- `has_ccc` (optional): Filter by CCC availability

**Example:**

```bash
# List all Visium datasets
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/analysis-features?dataset_type=visium"

# List datasets with CCC
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/analysis-features?has_ccc=true"
```

**Response:**

```json
[
  {
    "dataset_id": "ST024001",
    "has_umap": true,
    "has_tsne": false,
    "has_pca": true,
    "has_spatial_coordinates": true,
    "has_svg": true,
    "has_precomputed_deg": true,
    "has_deconvolution": true,
    "has_ccc": false,
    "dataset_type": "visium",
    "n_cells": 2615,
    "n_genes": 2306,
    "n_deg_groups": 10,
    "n_deconv_cell_types": 19,
    "id": 1,
    "created_at": "2025-12-05T02:51:19",
    "updated_at": null
  }
]
```

### Get Features for Specific Dataset

```bash
GET /api/v1/h5ad/analysis-features/{dataset_id}
```

**Example:**

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/analysis-features/ST024001"
```

### Sync Analysis Features

```bash
POST /api/v1/h5ad/analysis-features/{dataset_id}/sync
```

Automatically detects and updates analysis features by reading the h5ad file and cache.

**Example:**

```bash
curl -X POST "http://localhost:9117/sskind-backend/api/v1/h5ad/analysis-features/ST024001/sync"
```

## Frontend Integration Examples

### 1. Spatial Plot with Image Background (Recommended)

```javascript
// Get plot data with image (coordinates already transformed)
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data-with-image?genes=APOE&image_key=hires"
);
const data = await response.json();

// Load image
const img = new Image();
img.src = data.image_url;

img.onload = () => {
  // Create canvas
  const canvas = document.getElementById("spatial-canvas");
  const ctx = canvas.getContext("2d");
  canvas.width = data.image_shape[1];
  canvas.height = data.image_shape[0];

  // Draw image background
  ctx.drawImage(img, 0, 0);

  // Overlay coordinates (already in image space)
  data.coordinates.forEach((coord, i) => {
    const x = coord[0];
    const y = coord[1];
    const expr = data.genes.APOE[i];

    // Color by expression
    const color = getColorForExpression(expr);
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x, y, 3, 0, 2 * Math.PI);
    ctx.fill();
  });
};
```

### 2. Basic Spatial Plot (Without Image)

```javascript
// Get spatial coordinates (full resolution)
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/coordinates"
);
const data = await response.json();

// Plot with ECharts or similar
const x = data.coordinates.map((c) => c[0]);
const y = data.coordinates.map((c) => c[1]);

// Render scatter plot with spatial coordinates
```

### 3. Spatial Plot with Cell Type/Cluster Categories

```javascript
// Get spatial data with cell type/cluster annotations
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data-with-image?metadata=annotation,seurat_clusters"
);
const data = await response.json();

// Load image
const img = new Image();
img.src = data.image_url;

img.onload = () => {
  const canvas = document.getElementById("spatial-canvas");
  const ctx = canvas.getContext("2d");
  canvas.width = data.image_shape[1];
  canvas.height = data.image_shape[0];

  // Draw image background
  ctx.drawImage(img, 0, 0);

  // Create color map for cell types
  const cellTypes = [...new Set(data.metadata.annotation)];
  const colors = generateColors(cellTypes.length);
  const colorMap = {};
  cellTypes.forEach((type, i) => {
    colorMap[type] = colors[i];
  });

  // Overlay coordinates colored by cell type
  data.coordinates.forEach((coord, i) => {
    const x = coord[0];
    const y = coord[1];
    const cellType = data.metadata.annotation[i];
    const color = colorMap[cellType];

    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x, y, 3, 0, 2 * Math.PI);
    ctx.fill();
  });
};
```

### 4. Spatial Plot with Gene Expression on Image

```javascript
// Get transformed coordinates + gene expression + image
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data-with-image?genes=APOE"
);
const data = await response.json();

// Use image as background and overlay colored points
// Coordinates are already transformed to match image space
```

### 5. Combined: Cell Types + Gene Expression

```javascript
// Get both cell types and gene expression
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data-with-image?genes=APOE&metadata=annotation"
);
const data = await response.json();

// You can color by either:
// - Cell type: data.metadata.annotation[i]
// - Gene expression: data.genes.APOE[i]
```

### 3. Top Spatially Variable Genes

```javascript
// Get top 20 SVG
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/svg/top?n=20"
);
const svgData = await response.json();

// Display in table or list
svgData.genes.forEach((gene) => {
  console.log(`${gene.gene}: rank=${gene.svg_rank}, score=${gene.gft_score}`);
});
```

### 4. Deconvolution Visualization

```javascript
// Get deconvolution predictions
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/deconvolution/plot-data"
);
const data = await response.json();

// Visualize each cell type prediction on spatial coordinates
Object.keys(data.metadata).forEach((cellType) => {
  const predictions = data.metadata[cellType];
  // Create overlay for this cell type
});
```

### 5. Check Available Features

```javascript
// Check what features are available
const response = await fetch(
  "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/info"
);
const info = await response.json();

if (info.has_deconvolution) {
  // Load deconvolution data
}
if (info.has_ccc) {
  // Load CCC data
}
```

## Data Structure

### Spatial Coordinates

Stored in `ad.obsm['spatial']`, `ad.obsm['spatial_xy']`, or `ad.obsm['X_spatial']`

- Format: `[[x1, y1], [x2, y2], ...]`
- Units: Pixel coordinates or physical coordinates

### SVG Data

Stored in `ad.var` columns:

- `gft_score`: Global Fourier Transform score
- `svg_rank`: Rank by spatial variability
- `cutoff_gft_score`: Threshold score
- `pvalue`: Statistical significance
- `fdr`: False Discovery Rate

### Precomputed DEG

Stored in `ad.uns["rank_genes_groups"]`

- Contains: gene names, scores, log fold changes, p-values
- Multiple comparison groups supported

### Deconvolution

Stored in `ad.obsm["tangram_ct_pred"]`

- Matrix: cells × cell types
- Values: Prediction probabilities

### CCC Data

Stored in `ad.uns['commot_user_database']` or `ad.uns['commot']`

- Contains: ligand-receptor interactions, cell types, distance thresholds

## Preprocessing

The preprocessing script automatically extracts and caches all spatial features:

```bash
python scripts/preprocess_h5ad.py --h5ad h5ad/raw/ST024001.h5ad --dataset-id ST024001
```

This creates:

```
h5ad/precomputed/ST024001/
├── spatial_coordinates.json
├── svg_data.json
├── deg_precomputed.json
├── deconvolution.json      # Visium only
├── ccc_interactions.json   # Xenium only
└── spatial_info.json
```

## Performance

- **Spatial coordinates**: ~15ms (cached), ~150ms (first access)
- **SVG data**: ~20ms (cached), ~200ms (first access)
- **DEG results**: ~50ms (cached), ~500ms (first access)
- **Deconvolution**: ~30ms (cached), ~300ms (first access)
- **CCC data**: ~100ms (cached), ~1000ms (first access)

## Error Handling

All endpoints return appropriate HTTP status codes:

- `200`: Success
- `404`: Dataset not found or feature not available
- `500`: Server error

Example error response:

```json
{
  "detail": "Spatial coordinates not found in dataset ST024001"
}
```

## Best Practices

1. **Check features first**: Use `/spatial/info` to verify available features
2. **Use combined endpoints**: `/spatial/plot-data` is more efficient than multiple calls
3. **Cache on frontend**: Spatial coordinates don't change, cache them
4. **Filter SVG results**: Use `top_n` or `min_score` to limit results
5. **Sync after preprocessing**: Call `/analysis-features/{id}/sync` after preprocessing
6. **Gene symbols are case-sensitive**: The demo Visium dataset `ST024001` uses lower-case mouse symbols (e.g., `fth1`, `mt3`, `h2afz`; genes like `APOE/Plp1` are not present). If a requested gene is missing, expression arrays will be empty.

## Common Use Cases

### Find all Visium datasets with deconvolution

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/analysis-features?dataset_type=visium&has_deconvolution=true"
```

### Get top 10 SVG for visualization

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/svg/top?n=10"
```

### Visualize gene expression on spatial coordinates

```bash
curl "http://localhost:9117/sskind-backend/api/v1/h5ad/ST024001/spatial/plot-data?genes=APOE,APP"
```

## Related Documentation

- **H5AD API Guide**: `H5AD_API_GUIDE.md` - General h5ad endpoints
- **H5AD Architecture**: `H5AD_ARCHITECTURE_SUMMARY.md` - System architecture
- **API Summary**: `API_SUMMARY.md` - All API endpoints overview
- **Interactive Docs**: http://localhost:9117/sskind-backend/docs

## Troubleshooting

### Feature not available

- Check if dataset was preprocessed: `ls h5ad/precomputed/{dataset_id}/`
- Verify feature exists in h5ad file
- Check `/spatial/info` endpoint

### Slow performance

- Ensure preprocessing completed successfully
- Check cache files exist in `h5ad/precomputed/`
- Use combined endpoints instead of multiple calls

### Missing data

- Run preprocessing: `python preprocess_h5ad.py --h5ad h5ad/raw/{dataset_id}.h5ad`
- Sync features: `POST /analysis-features/{dataset_id}/sync`
