"""
API Metadata endpoint - Returns structured API documentation for frontend discovery
"""
from fastapi import APIRouter
from typing import Dict, List, Any
from app.config import settings

router = APIRouter()


@router.get("/api-metadata")
async def get_api_metadata() -> Dict[str, Any]:
    """
    Get structured API metadata for frontend discovery
    
    Returns comprehensive information about all available endpoints including:
    - Endpoint paths and methods
    - Parameters and their types
    - Response schemas
    - Examples
    - Feature categories
    """
    base_url = f"http://localhost:{settings.PORT}{settings.GLOBAL_PREFIX}{settings.API_PREFIX}"
    
    return {
        "base_url": base_url,
        "api_version": settings.APP_VERSION,
        "categories": {
            "spatial": get_spatial_endpoints(base_url),
            "general_h5ad": get_general_h5ad_endpoints(base_url),
            "analysis_features": get_analysis_features_endpoints(base_url),
            "deg": get_deg_endpoints(base_url)
        },
        "quick_reference": {
            "check_features": f"{base_url}/h5ad/{{dataset_id}}/spatial/info",
            "spatial_coords": f"{base_url}/h5ad/{{dataset_id}}/spatial/coordinates",
            "spatial_plot": f"{base_url}/h5ad/{{dataset_id}}/spatial/plot-data?genes=GENE1,GENE2",
            "list_datasets": f"{base_url}/h5ad/datasets",
            "list_features": f"{base_url}/h5ad/analysis-features"
        },
        "demo_datasets": {
            "visium": {
                "dataset_id": "ST024001",
                "type": "visium",
                "cells": 2615,
                "genes": 2306,
                "features": ["spatial_coordinates", "svg", "deg", "deconvolution"],
                "deg_groups": 10,
                "deconv_cell_types": 19
            },
            "xenium": {
                "dataset_id": "ST034001",
                "type": "xenium",
                "cells": 164081,
                "genes": 415,
                "features": ["spatial_coordinates", "svg", "deg", "ccc"],
                "deg_groups": 254
            }
        }
    }


def get_spatial_endpoints(base_url: str) -> List[Dict[str, Any]]:
    """Get spatial transcriptomics endpoints metadata"""
    return [
        {
            "name": "Get Spatial Info",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/info",
            "description": "Get information about available spatial features in the dataset",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier (e.g., ST024001)"
                }
            },
            "query_params": {},
            "response": {
                "type": "SpatialInfoResponse",
                "schema": {
                    "has_spatial_coordinates": "boolean",
                    "has_svg": "boolean",
                    "has_precomputed_deg": "boolean",
                    "has_deconvolution": "boolean",
                    "has_ccc": "boolean",
                    "dataset_type": "string | null (visium, xenium, unknown)"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/info",
                "response": {
                    "has_spatial_coordinates": True,
                    "has_svg": True,
                    "has_precomputed_deg": True,
                    "has_deconvolution": True,
                    "has_ccc": False,
                    "dataset_type": "visium"
                }
            }
        },
        {
            "name": "Get Spatial Coordinates",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/coordinates",
            "description": "Get spatial coordinates (x, y) for all spots/cells",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {},
            "response": {
                "type": "SpatialCoordinatesResponse",
                "schema": {
                    "coordinates": "array[array[float]]",
                    "cell_ids": "array[string]",
                    "spatial_key": "string | null"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/coordinates",
                "response": {
                    "coordinates": [[9575.0, 4911.0], [12671.0, 9443.0]],
                    "cell_ids": ["cell_1", "cell_2"],
                    "spatial_key": "spatial"
                }
            }
        },
        {
            "name": "Get Spatial Plot Data",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/plot-data",
            "description": "Get combined spatial plot data (coordinates + gene expression + metadata)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "genes": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated list of gene symbols",
                    "example": "APOE,APP"
                },
                "metadata": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated list of metadata columns (e.g., annotation, seurat_clusters). Returns full-resolution coordinates (no image transform).",
                    "example": "annotation,seurat_clusters"
                }
            },
            "response": {
                "type": "SpatialPlotDataResponse",
                "schema": {
                    "coordinates": "array[array[float]]",
                    "cell_ids": "array[string]",
                    "genes": "object<string, array[float]> | null",
                    "metadata": "object<string, array[any]> | null"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/plot-data?genes=APOE&metadata=annotation",
                "response": {
                    "coordinates": [[9575.0, 4911.0]],
                    "cell_ids": ["cell_1"],
                    "genes": {"APOE": [0.5, 1.2]},
                    "metadata": {"annotation": ["Neuron", "Astrocyte"]}
                }
            }
        },
        {
            "name": "Get Spatial Plot Coordinates with Image",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/plot-coordinates-with-image",
            "description": "Get transformed coordinates + metadata + image info (no gene expression). Use for image overlay.",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "sample_key": {
                    "type": "string",
                    "required": False,
                    "description": "Sample/slice key (default sample if omitted)"
                },
                "image_key": {
                    "type": "string",
                    "required": False,
                    "default": "hires",
                    "description": "Image type"
                },
                "metadata": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated metadata columns (e.g., annotation,seurat_clusters)"
                }
            },
            "response": {
                "type": "object",
                "schema": {
                    "coordinates": "array[array[float]] (transformed to image space)",
                    "cell_ids": "array[string]",
                    "sample_key": "string",
                    "image_key": "string",
                    "image_url": "string",
                    "image_shape": "array[int]",
                    "scale_factor": "float",
                    "original_coordinates": "array[array[float]]",
                    "scalefactors": "object",
                    "metadata": "object<string, array[any]> | null"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/plot-coordinates-with-image?metadata=annotation,seurat_clusters&sample_key=slice1"
            }
        },
        {
            "name": "Get Top SVG",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/svg/top",
            "description": "Get top N spatially variable genes",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "n": {
                    "type": "integer",
                    "required": False,
                    "default": 50,
                    "min": 1,
                    "max": 1000,
                    "description": "Number of top genes to return"
                }
            },
            "response": {
                "type": "SVGResponse",
                "schema": {
                    "genes": "array[SVGGeneResult]",
                    "total_genes": "integer",
                    "filtered_count": "integer",
                    "available_columns": "array[string]"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/svg/top?n=10",
                "response": {
                    "genes": [
                        {
                            "gene": "gm42418",
                            "gft_score": 4.23,
                            "svg_rank": 2.0,
                            "pvalue": 1.33e-77,
                            "fdr": 3.27e-76
                        }
                    ],
                    "total_genes": 2306,
                    "filtered_count": 10
                }
            }
        },
        {
            "name": "Get SVG List",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/svg/list",
            "description": "Get all spatially variable genes with optional filters",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "top_n": {
                    "type": "integer",
                    "required": False,
                    "description": "Return only top N genes by rank"
                },
                "min_score": {
                    "type": "float",
                    "required": False,
                    "description": "Minimum gft_score threshold"
                }
            },
            "response": {
                "type": "SVGResponse"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/svg/list?top_n=20&min_score=2.0"
            }
        },
        {
            "name": "Get SVG for Gene",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/svg/{gene}",
            "description": "Get SVG statistics for a specific gene",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                },
                "gene": {
                    "type": "string",
                    "required": True,
                    "description": "Gene symbol"
                }
            },
            "query_params": {},
            "response": {
                "type": "SVGGeneResult",
                "schema": {
                    "gene": "string",
                    "gft_score": "float | null",
                    "svg_rank": "float | null",
                    "cutoff_gft_score": "float | null",
                    "pvalue": "float | null",
                    "fdr": "float | null"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/svg/nrgn"
            }
        },
        {
            "name": "Get Precomputed DEG Groups",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/deg/precomputed/groups",
            "description": "List available comparison groups in precomputed DEG results",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {},
            "response": {
                "type": "object",
                "schema": {
                    "groups": "array[string]",
                    "method": "string"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/deg/precomputed/groups",
                "response": {
                    "groups": ["0", "1", "2"],
                    "method": "wilcoxon"
                }
            }
        },
        {
            "name": "Get Precomputed DEG",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/deg/precomputed",
            "description": "Get precomputed DEG results from h5ad file",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "group": {
                    "type": "string",
                    "required": False,
                    "description": "Filter by specific group name"
                }
            },
            "response": {
                "type": "PrecomputedDEGResponse"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/deg/precomputed?group=0"
            }
        },
        {
            "name": "Get Deconvolution Predictions",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/deconvolution/predictions",
            "description": "Get cell type predictions from Tangram deconvolution (Visium only)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {},
            "response": {
                "type": "DeconvolutionResponse",
                "schema": {
                    "predictions": "object<string, array[float]>",
                    "cell_ids": "array[string]",
                    "cell_types": "array[string]"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/deconvolution/predictions"
            },
            "dataset_types": ["visium"]
        },
        {
            "name": "Get Deconvolution Plot Data",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/deconvolution/plot-data",
            "description": "Get spatial coordinates + deconvolution predictions",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {},
            "response": {
                "type": "SpatialPlotDataResponse"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/deconvolution/plot-data"
            },
            "dataset_types": ["visium"]
        },
        {
            "name": "Get CCC Interactions",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/ccc/interactions",
            "description": "Get cell-cell communication interaction data (Xenium only)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {},
            "response": {
                "type": "CCCInteractionsResponse",
                "schema": {
                    "interactions": "array",
                    "ligands": "array[string]",
                    "receptors": "array[string]",
                    "cell_types": "array[string]",
                    "data": "object"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST034001/ccc/interactions"
            },
            "dataset_types": ["xenium"]
        },
        {
            "name": "Get CCC Network",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/ccc/network",
            "description": "Get simplified network representation of CCC",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {},
            "response": {
                "type": "object",
                "schema": {
                    "ligands": "array[string]",
                    "receptors": "array[string]",
                    "cell_types": "array[string]",
                    "interactions": "array"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST034001/ccc/network"
            },
            "dataset_types": ["xenium"]
        },
        {
            "name": "Get Spatial Image Info",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/image/info",
            "description": "Get spatial image metadata and available samples/slices",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "sample_key": {
                    "type": "string",
                    "required": False,
                    "description": "Optional specific sample/slice key to filter"
                }
            },
            "response": {
                "type": "SpatialImageInfoResponse",
                "schema": {
                    "samples": "object<string, SpatialImageSampleInfo>",
                    "default_sample": "string"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/image/info"
            }
        },
        {
            "name": "Get Spatial Image",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/image/{sample_key}/{image_key}",
            "description": "Get spatial tissue image file (PNG)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                },
                "sample_key": {
                    "type": "string",
                    "required": True,
                    "description": "Sample/slice key (e.g., 'slice1' for Visium)"
                },
                "image_key": {
                    "type": "string",
                    "required": True,
                    "default": "hires",
                    "description": "Image type (hires, lowres, etc.)"
                }
            },
            "query_params": {},
            "response": {
                "type": "FileResponse",
                "content_type": "image/png",
                "description": "PNG image file"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/image/slice1/hires"
            }
        },
        {
            "name": "Get Transformed Spatial Coordinates",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/coordinates/transformed",
            "description": "Get spatial coordinates transformed to match image space",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "sample_key": {
                    "type": "string",
                    "required": False,
                    "description": "Sample/slice key (uses default if not specified)"
                },
                "image_key": {
                    "type": "string",
                    "required": False,
                    "default": "hires",
                    "description": "Image type to transform for"
                }
            },
            "response": {
                "type": "SpatialCoordinatesTransformedResponse",
                "schema": {
                    "coordinates": "array[array[float]]",
                    "cell_ids": "array[string]",
                    "sample_key": "string",
                    "image_key": "string",
                    "scale_factor": "float",
                    "original_coordinates": "array[array[float]]",
                    "image_shape": "array[int] | null"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/coordinates/transformed?image_key=hires"
            }
        },
        {
            "name": "Get Spatial Plot Data with Image",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/spatial/plot-data-with-image",
            "description": "Get combined spatial plot data with image information (recommended for visualization)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True,
                    "description": "Dataset identifier"
                }
            },
            "query_params": {
                "sample_key": {
                    "type": "string",
                    "required": False,
                    "description": "Sample/slice key"
                },
                "image_key": {
                    "type": "string",
                    "required": False,
                    "default": "hires",
                    "description": "Image type"
                },
                "genes": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated list of genes"
                },
                "metadata": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated list of metadata columns"
                }
            },
            "response": {
                "type": "object",
                "schema": {
                    "coordinates": "array[array[float]]",
                    "cell_ids": "array[string]",
                    "image_url": "string",
                    "image_shape": "array[int]",
                    "scale_factor": "float",
                    "scalefactors": "object",
                    "genes": "object<string, array[float]> | null",
                    "metadata": "object<string, array[any]> | null"
                }
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/spatial/plot-data-with-image?genes=APOE&image_key=hires"
            }
        }
    ]


def get_general_h5ad_endpoints(base_url: str) -> List[Dict[str, Any]]:
    """Get general h5ad endpoints metadata"""
    return [
        {
            "name": "List Datasets",
            "method": "GET",
            "path": "/h5ad/datasets",
            "description": "List all available h5ad datasets",
            "path_params": {},
            "query_params": {},
            "response": {
                "type": "array[DatasetInfo]"
            },
            "example": {
                "url": f"{base_url}/h5ad/datasets"
            }
        },
        {
            "name": "Get Dataset Info",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/info",
            "description": "Get detailed information about a specific dataset",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True
                }
            },
            "query_params": {},
            "response": {
                "type": "DatasetInfo"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/info"
            }
        },
        {
            "name": "Get Embedding",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/embedding/{embedding_type}",
            "description": "Get embedding coordinates (UMAP, tSNE, PCA)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True
                },
                "embedding_type": {
                    "type": "string",
                    "required": True,
                    "enum": ["umap", "tsne", "pca"]
                }
            },
            "query_params": {},
            "response": {
                "type": "EmbeddingResponse"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/embedding/umap"
            }
        },
        {
            "name": "Get Gene Expression",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/expression/{gene}",
            "description": "Get gene expression for a specific gene",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True
                },
                "gene": {
                    "type": "string",
                    "required": True
                }
            },
            "query_params": {},
            "response": {
                "type": "GeneExpressionResponse"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/expression/APOE"
            }
        },
        {
            "name": "Get Plot Data",
            "method": "GET",
            "path": "/h5ad/{dataset_id}/plot-data",
            "description": "Get combined data for plotting (embedding + genes + metadata)",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True
                }
            },
            "query_params": {
                "embedding": {
                    "type": "string",
                    "required": False,
                    "default": "umap",
                    "enum": ["umap", "tsne", "pca"]
                },
                "genes": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated gene symbols"
                },
                "metadata": {
                    "type": "string",
                    "required": False,
                    "description": "Comma-separated metadata columns"
                }
            },
            "response": {
                "type": "PlotDataResponse"
            },
            "example": {
                "url": f"{base_url}/h5ad/ST024001/plot-data?embedding=umap&genes=APOE&metadata=cell_type"
            }
        }
    ]


def get_analysis_features_endpoints(base_url: str) -> List[Dict[str, Any]]:
    """Get analysis features management endpoints metadata"""
    return [
        {
            "name": "List Analysis Features",
            "method": "GET",
            "path": "/h5ad/analysis-features",
            "description": "List all datasets with their available analysis features",
            "path_params": {},
            "query_params": {
                "dataset_type": {
                    "type": "string",
                    "required": False,
                    "enum": ["visium", "xenium"],
                    "description": "Filter by dataset type"
                },
                "has_spatial": {
                    "type": "boolean",
                    "required": False,
                    "description": "Filter by spatial coordinates availability"
                },
                "has_svg": {
                    "type": "boolean",
                    "required": False,
                    "description": "Filter by SVG availability"
                },
                "has_deconvolution": {
                    "type": "boolean",
                    "required": False,
                    "description": "Filter by deconvolution availability"
                },
                "has_ccc": {
                    "type": "boolean",
                    "required": False,
                    "description": "Filter by CCC availability"
                }
            },
            "response": {
                "type": "array[H5ADAnalysisFeatures]"
            },
            "example": {
                "url": f"{base_url}/h5ad/analysis-features?dataset_type=visium&has_deconvolution=true"
            }
        },
        {
            "name": "Get Analysis Features",
            "method": "GET",
            "path": "/h5ad/analysis-features/{dataset_id}",
            "description": "Get analysis features for a specific dataset",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True
                }
            },
            "query_params": {},
            "response": {
                "type": "H5ADAnalysisFeatures"
            },
            "example": {
                "url": f"{base_url}/h5ad/analysis-features/ST024001"
            }
        },
        {
            "name": "Sync Analysis Features",
            "method": "POST",
            "path": "/h5ad/analysis-features/{dataset_id}/sync",
            "description": "Automatically detect and update analysis features from h5ad file",
            "path_params": {
                "dataset_id": {
                    "type": "string",
                    "required": True
                }
            },
            "query_params": {},
            "response": {
                "type": "H5ADAnalysisFeatures"
            },
            "example": {
                "url": f"{base_url}/h5ad/analysis-features/ST024001/sync"
            }
        }
    ]


def get_deg_endpoints(base_url: str) -> List[Dict[str, Any]]:
    """Get DEG analysis endpoints metadata"""
    return [
        {
            "name": "DEG Between Datasets",
            "method": "POST",
            "path": "/deg/between-datasets",
            "description": "Perform DEG analysis between two datasets",
            "path_params": {},
            "query_params": {},
            "request_body": {
                "type": "DEGBetweenDatasetsRequest",
                "schema": {
                    "dataset_id1": "string",
                    "dataset_id2": "string",
                    "cell_type": "string | null",
                    "min_pct": "float (default: 0.1)",
                    "logfc_threshold": "float (default: 0.25)",
                    "p_value_threshold": "float (default: 0.05)",
                    "top_n": "integer | null"
                }
            },
            "response": {
                "type": "DEGAnalysisResult"
            },
            "example": {
                "url": f"{base_url}/deg/between-datasets",
                "body": {
                    "dataset_id1": "AD093044",
                    "dataset_id2": "AD093045",
                    "cell_type": "Neuron",
                    "top_n": 100
                }
            }
        },
        {
            "name": "DEG Within Dataset",
            "method": "POST",
            "path": "/deg/within-dataset",
            "description": "Perform DEG analysis within a single dataset",
            "path_params": {},
            "query_params": {},
            "request_body": {
                "type": "DEGWithinDatasetRequest",
                "schema": {
                    "dataset_id": "string",
                    "group1_filters": "object<string, string>",
                    "group2_filters": "object<string, string>",
                    "min_pct": "float",
                    "logfc_threshold": "float",
                    "p_value_threshold": "float",
                    "top_n": "integer | null"
                }
            },
            "response": {
                "type": "DEGAnalysisResult"
            },
            "example": {
                "url": f"{base_url}/deg/within-dataset",
                "body": {
                    "dataset_id": "AD093044",
                    "group1_filters": {"cell_type": "Neuron"},
                    "group2_filters": {"cell_type": "Astrocyte"}
                }
            }
        }
    ]

