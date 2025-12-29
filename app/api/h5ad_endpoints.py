"""
API endpoints for H5AD data visualization
"""
from fastapi import APIRouter, HTTPException, Query, Depends, Response
from fastapi.responses import FileResponse
from typing import List, Optional
import logging
from pathlib import Path
from sqlalchemy.orm import Session

from app.services.h5ad_service import H5ADService
from app.database import get_db
from app import models
from app.schemas import (
    EmbeddingResponse,
    GeneExpressionResponse,
    MetadataResponse,
    MetadataValuesResponse,
    PlotDataResponse,
    GeneSearchResult,
    ModuleScoreResponse,
    ModuleScoreRequest,
    HeatmapRequest,
    HeatmapResponse,
    SpatialCoordinatesResponse,
    SpatialPlotDataResponse,
    SVGResponse,
    SVGGeneResult,
    PrecomputedDEGResponse,
    DeconvolutionResponse,
    CCCInteractionsResponse,
    SpatialInfoResponse,
    SpatialImageInfoResponse,
    SpatialCoordinatesTransformedResponse,
    H5ADAnalysisFeatures,
    H5ADAnalysisFeaturesCreate,
    H5ADAnalysisFeaturesUpdate
)

logger = logging.getLogger(__name__)

router = APIRouter()
h5ad_service = H5ADService()


@router.get("/datasets")
async def list_h5ad_datasets():
    """
    List all available h5ad datasets
    """
    # TODO: Query from database when h5ad_files table is populated
    # For now, scan the h5ad/raw directory
    import os
    from pathlib import Path
    
    h5ad_dir = Path("h5ad/raw")
    if not h5ad_dir.exists():
        return []
    
    datasets = []
    for file in h5ad_dir.glob("*.h5ad"):
        dataset_id = file.stem
        try:
            info = h5ad_service.get_dataset_info(dataset_id)
            datasets.append(info)
        except Exception as e:
            logger.error(f"Error loading info for {dataset_id}: {e}")
    
    return datasets


@router.get("/{dataset_id}/info")
async def get_dataset_info(dataset_id: str):
    """
    Get detailed information about a specific dataset
    """
    try:
        return h5ad_service.get_dataset_info(dataset_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error getting dataset info: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/embedding/{embedding_type}", response_model=EmbeddingResponse)
async def get_embedding(
    dataset_id: str,
    embedding_type: str
):
    """
    Get embedding coordinates (UMAP, tSNE, PCA)
    
    Args:
        dataset_id: Dataset identifier
        embedding_type: Type of embedding (umap, tsne, pca)
    """
    try:
        data = h5ad_service.get_embedding(dataset_id, embedding_type)
        return EmbeddingResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting embedding: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/expression/{gene}", response_model=GeneExpressionResponse)
async def get_gene_expression(
    dataset_id: str,
    gene: str
):
    """
    Get gene expression for a specific gene
    
    Args:
        dataset_id: Dataset identifier
        gene: Gene symbol
    """
    try:
        data = h5ad_service.get_gene_expression(dataset_id, gene)
        return GeneExpressionResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting gene expression: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/metadata", response_model=MetadataResponse)
async def get_metadata(
    dataset_id: str,
    columns: Optional[str] = Query(None, description="Comma-separated list of metadata columns")
):
    """
    Get cell metadata
    
    Args:
        dataset_id: Dataset identifier
        columns: Optional comma-separated list of specific columns to return
    """
    try:
        column_list = columns.split(",") if columns else None
        data = h5ad_service.get_metadata(dataset_id, column_list)
        return MetadataResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error getting metadata: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/metadata-values/{column}", response_model=MetadataValuesResponse)
async def get_metadata_values(
    dataset_id: str,
    column: str
):
    """
    Get unique values and counts for a specific metadata column
    
    This endpoint is useful for building filter UIs in the frontend.
    Returns all unique values for the specified column along with cell counts.
    
    Args:
        dataset_id: Dataset identifier
        column: Metadata column name
        
    Example:
        GET /h5ad/human_subset/metadata-values/cell_type
        
    Returns:
        {
            "column": "cell_type",
            "values": ["Astrocyte", "Neuron", "Microglia", ...],
            "counts": {
                "Astrocyte": 15000,
                "Neuron": 45000,
                "Microglia": 8000,
                ...
            }
        }
    """
    try:
        data = h5ad_service.get_metadata_column_values(dataset_id, column)
        return MetadataValuesResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting metadata values: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/plot-data", response_model=PlotDataResponse)
async def get_plot_data(
    dataset_id: str,
    embedding: str = Query("umap", description="Embedding type (umap, tsne, pca)"),
    genes: Optional[str] = Query(None, description="Comma-separated list of genes"),
    metadata: Optional[str] = Query(None, description="Comma-separated list of metadata columns")
):
    """
    Get combined data for plotting (optimized endpoint)
    
    This endpoint returns everything needed for visualization in one request:
    - Embedding coordinates
    - Gene expression (if genes specified)
    - Metadata (if metadata columns specified)
    
    Args:
        dataset_id: Dataset identifier
        embedding: Embedding type (umap, tsne, pca)
        genes: Optional comma-separated gene symbols
        metadata: Optional comma-separated metadata columns
    """
    try:
        gene_list = genes.split(",") if genes else None
        metadata_list = metadata.split(",") if metadata else None
        
        data = h5ad_service.get_plot_data(
            dataset_id=dataset_id,
            embedding=embedding,
            genes=gene_list,
            metadata_cols=metadata_list
        )
        
        return PlotDataResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting plot data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/genes/search", response_model=List[GeneSearchResult])
async def search_genes(
    dataset_id: str,
    q: str = Query(..., description="Search query"),
    limit: int = Query(50, ge=1, le=500, description="Maximum number of results")
):
    """
    Search for genes by symbol
    
    Args:
        dataset_id: Dataset identifier
        q: Search query (gene symbol or partial match)
        limit: Maximum number of results (1-500)
    """
    try:
        results = h5ad_service.search_genes(dataset_id, q, limit)
        return [GeneSearchResult(**r) for r in results]
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error searching genes: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{dataset_id}/module-score", response_model=ModuleScoreResponse)
async def calculate_module_score(
    dataset_id: str,
    request: ModuleScoreRequest
):
    """
    Calculate module score for a list of genes using scanpy.tl.score_genes
    
    **Module Score**: Seurat-style scoring method that calculates the average expression
    of a set of genes after subtraction by the average expression of a reference set of genes.
    The reference set is randomly sampled from the gene pool for each binned expression value.
    
    This reproduces the approach in Seurat [Satija et al., 2015] and has been implemented
    for Scanpy by Davide Cittaro.
    
    **Use Cases:**
    - Gene signature scoring (e.g., immune signature, cell cycle signature)
    - Pathway activity scoring
    - Multi-gene biomarker evaluation
    - Cell state characterization
    
    **Example:**
    ```json
    {
      "gene_list": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"],
      "use_raw": false
    }
    ```
    
    **Parameters:**
    - **gene_list**: List of gene symbols (e.g., ["APOE", "APP", "MAPT"])
    - **use_raw**: Use raw counts if True, normalized if False (default: False, recommended)
    
    **Returns:**
    - Module score for each cell (can be negative, as it's relative to reference set)
    - Which genes were found/not found
    - Statistics (mean, median, std, etc.)
    - Method used: "scanpy_score_genes"
    """
    try:
        result = h5ad_service.calculate_module_score(
            dataset_id=dataset_id,
            gene_list=request.gene_list,
            use_raw=request.use_raw
        )
        return ModuleScoreResponse(**result)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error calculating module score: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{dataset_id}/heatmap", response_model=HeatmapResponse)
async def generate_heatmap(
    dataset_id: str,
    request: HeatmapRequest
):
    """
    Generate heatmap or dotplot data for genes across cell types
    
    **Heatmap**: Shows average expression values (log-transformed) for each gene across each cell type
    
    **Dotplot**: Shows both expression level and percentage of cells expressing each gene
    
    **Example (Heatmap):**
    ```json
    {
      "gene_list": ["APOE", "APP", "MAPT", "PSEN1", "PSEN2"],
      "plot_type": "heatmap"
    }
    ```
    
    **Example (Dotplot):**
    ```json
    {
      "gene_list": ["APOE", "APP", "MAPT"],
      "plot_type": "dotplot"
    }
    ```
    
    **Parameters:**
    - **gene_list**: List of gene symbols to visualize
    - **plot_type**: "heatmap" or "dotplot" (default: "heatmap")
    - **scale_expression**: Z-score normalize each gene across cell types (default: false)
    - **cluster_rows**: Cluster genes using hierarchical clustering (default: false)
    - **cluster_columns**: Cluster cell types using hierarchical clustering (default: false)
    
    **Returns:**
    - 2D data matrix: `data[gene_index][cell_type_index]`
    - For heatmap: array of expression values (log-transformed)
    - For dotplot: array of objects with `expression` and `percentage`
    - List of genes and cell types
    - Which genes were found/not found
    - Summary statistics
    
    **Note**: Requires 'cell_type' column in metadata
    """
    try:
        result = h5ad_service.generate_heatmap(
            dataset_id=dataset_id,
            gene_list=request.gene_list,
            plot_type=request.plot_type,
            scale_expression=request.scale_expression,
            cluster_rows=request.cluster_rows,
            cluster_columns=request.cluster_columns,
            group_by=request.group_by
        )
        return HeatmapResponse(**result)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error generating heatmap: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Spatial Transcriptomics Endpoints

@router.get("/{dataset_id}/spatial/info", response_model=SpatialInfoResponse)
async def get_spatial_info(dataset_id: str):
    """
    Get information about available spatial features in the dataset
    
    Returns information about what spatial features are available:
    - Spatial coordinates
    - Spatially Variable Genes (SVG)
    - Precomputed DEG
    - Deconvolution predictions (Visium)
    - Cell-Cell Communication (CCC) data (Xenium)
    - Dataset type (visium, xenium, unknown)
    - Spatial image availability
    """
    try:
        info = h5ad_service.detect_spatial_features(dataset_id)
        
        # Check for spatial image
        try:
            image_info = h5ad_service.get_spatial_image_info(dataset_id)
            info["has_spatial_image"] = len(image_info.get("samples", {})) > 0
        except:
            info["has_spatial_image"] = False
        
        return SpatialInfoResponse(**info)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error getting spatial info: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/spatial/coordinates", response_model=SpatialCoordinatesResponse)
async def get_spatial_coordinates(dataset_id: str):
    """
    Get spatial coordinates (x, y) for spots/cells
    
    Returns spatial coordinates that can be used to visualize data in physical space.
    Checks multiple possible keys: spatial, spatial_xy, X_spatial
    """
    try:
        data = h5ad_service.get_spatial_coordinates(dataset_id)
        return SpatialCoordinatesResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting spatial coordinates: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/spatial/plot-data", response_model=SpatialPlotDataResponse)
async def get_spatial_plot_data(
    dataset_id: str,
    genes: Optional[str] = Query(None, description="Comma-separated list of genes"),
    metadata: Optional[str] = Query(None, description="Comma-separated list of metadata columns")
):
    """
    Get combined spatial plot data (coordinates + gene expression + metadata)
    
    **Note:** This endpoint returns FULL-RESOLUTION coordinates (not transformed for image overlay).
    For image overlay visualization, use `/spatial/plot-data-with-image` instead.
    
    This endpoint returns everything needed for spatial visualization in one request:
    - Spatial coordinates (x, y) - FULL RESOLUTION (not transformed)
    - Gene expression (if genes specified)
    - Metadata (if metadata columns specified)
    
    **Key differences from `/spatial/plot-data-with-image`:**
    - Returns full-resolution coordinates (e.g., [9575.0, 4911.0])
    - Does NOT include image URL or transformation info
    - Does NOT accept sample_key or image_key parameters
    - Use this when you don't need image overlay
    
    Args:
        dataset_id: Dataset identifier
        genes: Optional comma-separated gene symbols
        metadata: Optional comma-separated metadata columns
    """
    try:
        gene_list = genes.split(",") if genes else None
        metadata_list = metadata.split(",") if metadata else None
        
        data = h5ad_service.get_spatial_plot_data(
            dataset_id=dataset_id,
            genes=gene_list,
            metadata_cols=metadata_list
        )
        
        return SpatialPlotDataResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting spatial plot data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/svg/list", response_model=SVGResponse)
async def get_svg_list(
    dataset_id: str,
    top_n: Optional[int] = Query(None, ge=1, description="Return only top N genes by rank"),
    min_score: Optional[float] = Query(None, description="Minimum gft_score threshold")
):
    """
    Get Spatially Variable Genes (SVG) list
    
    Returns all spatially variable genes with their scores, ranks, and statistics.
    Extracts from ad.var columns: gft_score, svg_rank, cutoff_gft_score, pvalue, fdr
    
    Args:
        dataset_id: Dataset identifier
        top_n: Optional, return only top N genes by rank
        min_score: Optional, minimum gft_score threshold
    """
    try:
        data = h5ad_service.get_svg_data(dataset_id, top_n=top_n, min_score=min_score)
        return SVGResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting SVG data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/svg/top", response_model=SVGResponse)
async def get_svg_top(
    dataset_id: str,
    n: int = Query(50, ge=1, le=1000, description="Number of top genes to return")
):
    """
    Get top N spatially variable genes
    
    Convenience endpoint to get top N SVG genes sorted by rank.
    """
    try:
        data = h5ad_service.get_svg_data(dataset_id, top_n=n)
        return SVGResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting top SVG genes: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/svg/{gene}", response_model=SVGGeneResult)
async def get_svg_gene(dataset_id: str, gene: str):
    """
    Get SVG statistics for a specific gene
    
    Returns SVG data (gft_score, svg_rank, pvalue, fdr) for a single gene.
    """
    try:
        data = h5ad_service.get_svg_data(dataset_id)
        # Find the gene
        for gene_data in data["genes"]:
            if gene_data["gene"] == gene:
                return SVGGeneResult(**gene_data)
        
        raise HTTPException(status_code=404, detail=f"Gene {gene} not found in SVG data")
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting SVG gene data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/deg/precomputed/groups")
async def get_precomputed_deg_groups(dataset_id: str):
    """
    List available comparison groups in precomputed DEG results
    
    Returns list of group names available in ad.uns["rank_genes_groups"]
    """
    try:
        data = h5ad_service.get_precomputed_deg(dataset_id)
        return {"groups": data["groups"], "method": data["method"]}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting DEG groups: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/deg/precomputed", response_model=PrecomputedDEGResponse)
async def get_precomputed_deg(
    dataset_id: str,
    group: Optional[str] = Query(None, description="Filter by specific group name")
):
    """
    Get precomputed DEG results from ad.uns["rank_genes_groups"]
    
    Returns differentially expressed genes that were precomputed and stored in the h5ad file.
    This is different from the dynamic DEG analysis endpoints - this retrieves existing results.
    
    Args:
        dataset_id: Dataset identifier
        group: Optional group name to filter (e.g., '0', '1', etc.)
    """
    try:
        data = h5ad_service.get_precomputed_deg(dataset_id, group=group)
        return PrecomputedDEGResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting precomputed DEG: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/deconvolution/predictions", response_model=DeconvolutionResponse)
async def get_deconvolution_predictions(dataset_id: str):
    """
    Get deconvolution (Tangram) cell type predictions
    
    Returns cell type predictions from Tangram deconvolution analysis.
    Extracts from ad.obsm["tangram_ct_pred"]
    
    Available for Visium datasets.
    """
    try:
        data = h5ad_service.get_deconvolution_predictions(dataset_id)
        return DeconvolutionResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting deconvolution predictions: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/deconvolution/plot-data", response_model=SpatialPlotDataResponse)
async def get_deconvolution_plot_data(dataset_id: str):
    """
    Get spatial coordinates + deconvolution predictions
    
    Combined endpoint for visualizing deconvolution results on spatial coordinates.
    """
    try:
        # Get spatial coordinates
        spatial_data = h5ad_service.get_spatial_coordinates(dataset_id)
        
        # Get deconvolution predictions
        deconv_data = h5ad_service.get_deconvolution_predictions(dataset_id)
        
        # Combine into plot data format
        result = {
            "coordinates": spatial_data["coordinates"],
            "cell_ids": spatial_data["cell_ids"],
            "metadata": deconv_data["predictions"]  # Use predictions as metadata
        }
        
        return SpatialPlotDataResponse(**result)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting deconvolution plot data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/ccc/interactions", response_model=CCCInteractionsResponse)
async def get_ccc_interactions(dataset_id: str):
    """
    Get Cell-Cell Communication (CCC) interaction data
    
    Returns CCC interaction data from COMMOT analysis.
    Extracts from ad.uns['commot_user_database'] or ad.uns['commot']
    
    Available for Xenium datasets.
    """
    try:
        data = h5ad_service.get_ccc_interactions(dataset_id)
        return CCCInteractionsResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting CCC interactions: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/ccc/network")
async def get_ccc_network(dataset_id: str):
    """
    Get CCC network structure
    
    Returns a simplified network representation of cell-cell communication interactions.
    """
    try:
        data = h5ad_service.get_ccc_interactions(dataset_id)
        
        # Extract network structure
        network = {
            "ligands": data.get("ligands", []),
            "receptors": data.get("receptors", []),
            "cell_types": data.get("cell_types", []),
            "interactions": data.get("interactions", [])
        }
        
        return network
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting CCC network: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Spatial Image Endpoints

@router.get("/{dataset_id}/spatial/image/info", response_model=SpatialImageInfoResponse)
async def get_spatial_image_info(
    dataset_id: str,
    sample_key: Optional[str] = Query(None, description="Optional specific sample/slice key")
):
    """
    Get spatial image metadata and available samples
    
    Returns information about available spatial images including:
    - Available sample keys (slices/regions)
    - Image dimensions for each image type
    - Scale factors for coordinate transformation
    - Paths to cached image files
    
    Args:
        dataset_id: Dataset identifier
        sample_key: Optional specific sample key to filter
    """
    try:
        info = h5ad_service.get_spatial_image_info(dataset_id, sample_key)
        return SpatialImageInfoResponse(**info)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting spatial image info: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/spatial/image/{sample_key}/{image_key}")
async def get_spatial_image(
    dataset_id: str,
    sample_key: str,
    image_key: str = "hires"
):
    """
    Get spatial tissue image file
    
    Returns the PNG image file for the specified sample and image type.
    The image can be used as a background for overlaying spatial coordinates.
    
    Args:
        dataset_id: Dataset identifier
        sample_key: Sample/slice key (e.g., 'slice1' for Visium, region name for Xenium)
        image_key: Image type (default: 'hires', can be 'lowres' or other available types)
    """
    try:
        # Get image info to find the path
        image_info = h5ad_service.get_spatial_image_info(dataset_id, sample_key)
        
        if sample_key not in image_info["samples"]:
            raise HTTPException(
                status_code=404,
                detail=f"Sample key '{sample_key}' not found. Available: {list(image_info['samples'].keys())}"
            )
        
        sample_info = image_info["samples"][sample_key]
        
        if image_key not in sample_info.get("image_paths", {}):
            raise HTTPException(
                status_code=404,
                detail=f"Image key '{image_key}' not found. Available: {sample_info.get('image_keys', [])}"
            )
        
        # Get full path from relative path stored in cache
        rel_path = sample_info["image_paths"][image_key]
        precomputed_dir = h5ad_service.get_precomputed_dir(dataset_id)
        image_path = precomputed_dir / rel_path
        
        if not image_path.exists():
            raise HTTPException(status_code=404, detail=f"Image file not found: {image_path}")
        
        return FileResponse(
            path=str(image_path),
            media_type="image/png",
            filename=f"{dataset_id}_{sample_key}_{image_key}.png"
        )
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting spatial image: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/spatial/coordinates/transformed", response_model=SpatialCoordinatesTransformedResponse)
async def get_spatial_coordinates_transformed(
    dataset_id: str,
    sample_key: Optional[str] = Query(None, description="Sample/slice key (uses default if not specified)"),
    image_key: str = Query("hires", description="Image type to transform for (hires, lowres, etc.)")
):
    """
    Get spatial coordinates transformed to match image space
    
    Returns coordinates that are scaled to match the spatial image dimensions.
    This is essential for overlaying coordinates on the tissue image.
    
    The transformation applies: image_coord = fullres_coord * tissue_hires_scalef
    
    Args:
        dataset_id: Dataset identifier
        sample_key: Optional sample key (uses default if None)
        image_key: Image key to transform for (default: "hires")
    
    Returns:
        - coordinates: Transformed coordinates in image space
        - original_coordinates: Original full-resolution coordinates
        - scale_factor: Applied scale factor
        - image_shape: Image dimensions [height, width, channels]
    """
    try:
        data = h5ad_service.get_spatial_coordinates_transformed(
            dataset_id=dataset_id,
            sample_key=sample_key,
            image_key=image_key
        )
        return SpatialCoordinatesTransformedResponse(**data)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting transformed coordinates: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/spatial/plot-data-with-image")
async def get_spatial_plot_data_with_image(
    dataset_id: str,
    sample_key: Optional[str] = Query(None, description="Sample/slice key"),
    image_key: str = Query("hires", description="Image type"),
    genes: Optional[str] = Query(None, description="Comma-separated list of genes"),
    metadata: Optional[str] = Query(None, description="Comma-separated list of metadata columns")
):
    """
    Get combined spatial plot data with image information
    
    Returns everything needed for spatial visualization with image overlay:
    - Transformed coordinates (matching image space)
    - Image metadata and URL
    - Gene expression (if genes specified)
    - Metadata (if metadata columns specified)
    
    This is the recommended endpoint for spatial visualization with image background.
    """
    try:
        # Get transformed coordinates
        coords_data = h5ad_service.get_spatial_coordinates_transformed(
            dataset_id=dataset_id,
            sample_key=sample_key,
            image_key=image_key
        )
        
        # Get image info
        image_info = h5ad_service.get_spatial_image_info(dataset_id, coords_data["sample_key"])
        sample_info = image_info["samples"][coords_data["sample_key"]]
        
        # Build image URL (use full path with base URL)
        from app.config import settings
        base_url = f"http://localhost:{settings.PORT}{settings.GLOBAL_PREFIX}{settings.API_PREFIX}"
        image_url = f"{base_url}/h5ad/{dataset_id}/spatial/image/{coords_data['sample_key']}/{image_key}"
        
        result = {
            "coordinates": coords_data["coordinates"],
            "cell_ids": coords_data["cell_ids"],
            "sample_key": coords_data["sample_key"],
            "image_key": image_key,
            "image_url": image_url,
            "image_shape": coords_data["image_shape"],
            "scale_factor": coords_data["scale_factor"],
            "original_coordinates": coords_data["original_coordinates"],
            "scalefactors": sample_info.get("scalefactors", {})
        }
        
        # Get gene expression if requested
        if genes:
            gene_list = genes.split(",")
            result["genes"] = {}
            for gene in gene_list:
                try:
                    gene_data = h5ad_service.get_gene_expression(dataset_id, gene)
                    result["genes"][gene] = gene_data["expression"]
                except ValueError as e:
                    logger.warning(f"Gene {gene} not found: {e}")
                    # Return empty list to satisfy schema and signal missing gene
                    result["genes"][gene] = []
        
        # Get metadata if requested
        if metadata:
            metadata_list = metadata.split(",")
            metadata_data = h5ad_service.get_metadata(dataset_id, metadata_list)
            result["metadata"] = metadata_data["data"]
        
        return result
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting spatial plot data with image: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}/spatial/plot-coordinates-with-image")
async def get_spatial_plot_coordinates_with_image(
    dataset_id: str,
    sample_key: Optional[str] = Query(None, description="Sample/slice key"),
    image_key: str = Query("hires", description="Image type"),
    metadata: Optional[str] = Query(None, description="Comma-separated list of metadata columns")
):
    """
    Get transformed spatial coordinates with image info (no gene expression).

    Use this for image overlay when you only need coordinates + metadata.
    Fetch gene expression separately (e.g., /expression endpoints) to keep
    responsibilities single-purpose.
    """
    try:
        coords_data = h5ad_service.get_spatial_coordinates_transformed(
            dataset_id=dataset_id,
            sample_key=sample_key,
            image_key=image_key
        )

        image_info = h5ad_service.get_spatial_image_info(dataset_id, coords_data["sample_key"])
        sample_info = image_info["samples"][coords_data["sample_key"]]

        from app.config import settings
        base_url = f"http://localhost:{settings.PORT}{settings.GLOBAL_PREFIX}{settings.API_PREFIX}"
        image_url = f"{base_url}/h5ad/{dataset_id}/spatial/image/{coords_data['sample_key']}/{image_key}"

        result = {
            "coordinates": coords_data["coordinates"],
            "cell_ids": coords_data["cell_ids"],
            "sample_key": coords_data["sample_key"],
            "image_key": image_key,
            "image_url": image_url,
            "image_shape": coords_data["image_shape"],
            "scale_factor": coords_data["scale_factor"],
            "original_coordinates": coords_data["original_coordinates"],
            "scalefactors": sample_info.get("scalefactors", {})
        }

        if metadata:
            metadata_list = metadata.split(",")
            metadata_data = h5ad_service.get_metadata(dataset_id, metadata_list)
            result["metadata"] = metadata_data["data"]

        return result
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error getting spatial plot coordinates with image: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Analysis Features Management Endpoints

@router.get("/analysis-features", response_model=List[H5ADAnalysisFeatures])
async def list_analysis_features(
    db: Session = Depends(get_db),
    dataset_type: Optional[str] = Query(None, description="Filter by dataset type (visium, xenium)"),
    has_spatial: Optional[bool] = Query(None, description="Filter by spatial coordinates availability"),
    has_svg: Optional[bool] = Query(None, description="Filter by SVG availability"),
    has_deconvolution: Optional[bool] = Query(None, description="Filter by deconvolution availability"),
    has_ccc: Optional[bool] = Query(None, description="Filter by CCC availability")
):
    """
    List all datasets with their available analysis features
    
    Returns a list of all h5ad datasets and what analysis results are available for each.
    Can filter by various features.
    """
    query = db.query(models.H5ADAnalysisFeatures)
    
    if dataset_type:
        query = query.filter(models.H5ADAnalysisFeatures.dataset_type == dataset_type)
    if has_spatial is not None:
        query = query.filter(models.H5ADAnalysisFeatures.has_spatial_coordinates == has_spatial)
    if has_svg is not None:
        query = query.filter(models.H5ADAnalysisFeatures.has_svg == has_svg)
    if has_deconvolution is not None:
        query = query.filter(models.H5ADAnalysisFeatures.has_deconvolution == has_deconvolution)
    if has_ccc is not None:
        query = query.filter(models.H5ADAnalysisFeatures.has_ccc == has_ccc)
    
    features = query.all()
    return features


@router.get("/analysis-features/{dataset_id}", response_model=H5ADAnalysisFeatures)
async def get_analysis_features(dataset_id: str, db: Session = Depends(get_db)):
    """
    Get analysis features for a specific dataset
    
    Returns what analysis results are available for the specified dataset.
    """
    features = db.query(models.H5ADAnalysisFeatures).filter(
        models.H5ADAnalysisFeatures.dataset_id == dataset_id
    ).first()
    
    if not features:
        raise HTTPException(status_code=404, detail=f"Analysis features not found for dataset {dataset_id}")
    
    return features


@router.post("/analysis-features", response_model=H5ADAnalysisFeatures)
async def create_analysis_features(
    features: H5ADAnalysisFeaturesCreate,
    db: Session = Depends(get_db)
):
    """
    Create or update analysis features for a dataset
    
    This endpoint is typically called automatically during preprocessing,
    but can be used to manually register analysis features.
    """
    # Check if exists
    existing = db.query(models.H5ADAnalysisFeatures).filter(
        models.H5ADAnalysisFeatures.dataset_id == features.dataset_id
    ).first()
    
    if existing:
        # Update existing
        for key, value in features.model_dump(exclude_unset=True).items():
            setattr(existing, key, value)
        db.commit()
        db.refresh(existing)
        return existing
    else:
        # Create new
        db_features = models.H5ADAnalysisFeatures(**features.model_dump())
        db.add(db_features)
        db.commit()
        db.refresh(db_features)
        return db_features


@router.put("/analysis-features/{dataset_id}", response_model=H5ADAnalysisFeatures)
async def update_analysis_features(
    dataset_id: str,
    features: H5ADAnalysisFeaturesUpdate,
    db: Session = Depends(get_db)
):
    """
    Update analysis features for a dataset
    """
    db_features = db.query(models.H5ADAnalysisFeatures).filter(
        models.H5ADAnalysisFeatures.dataset_id == dataset_id
    ).first()
    
    if not db_features:
        raise HTTPException(status_code=404, detail=f"Analysis features not found for dataset {dataset_id}")
    
    # Update only provided fields
    for key, value in features.model_dump(exclude_unset=True).items():
        setattr(db_features, key, value)
    
    db.commit()
    db.refresh(db_features)
    return db_features


@router.post("/analysis-features/{dataset_id}/sync")
async def sync_analysis_features(dataset_id: str, db: Session = Depends(get_db)):
    """
    Sync analysis features from h5ad file to database
    
    Automatically detects and updates analysis features by reading the h5ad file
    and precomputed cache. Useful for updating features after preprocessing.
    """
    try:
        # Get spatial info from service
        spatial_info = h5ad_service.detect_spatial_features(dataset_id)
        
        # Get dataset info
        dataset_info = h5ad_service.get_dataset_info(dataset_id)
        
        # Get additional counts
        n_deg_groups = None
        n_deconv_cell_types = None
        
        try:
            deg_data = h5ad_service.get_precomputed_deg(dataset_id)
            n_deg_groups = len(deg_data.get("groups", []))
        except:
            pass
        
        try:
            deconv_data = h5ad_service.get_deconvolution_predictions(dataset_id)
            n_deconv_cell_types = len(deconv_data.get("cell_types", []))
        except:
            pass
        
        # Check embeddings
        embeddings = dataset_info.get("embeddings", [])
        
        # Create/update features
        features_data = {
            "dataset_id": dataset_id,
            "has_umap": "umap" in embeddings,
            "has_tsne": "tsne" in embeddings,
            "has_pca": "pca" in embeddings,
            "has_spatial_coordinates": spatial_info.get("has_spatial_coordinates", False),
            "has_svg": spatial_info.get("has_svg", False),
            "has_precomputed_deg": spatial_info.get("has_precomputed_deg", False),
            "has_deconvolution": spatial_info.get("has_deconvolution", False),
            "has_ccc": spatial_info.get("has_ccc", False),
            "dataset_type": spatial_info.get("dataset_type"),
            "n_cells": dataset_info.get("n_cells"),
            "n_genes": dataset_info.get("n_genes"),
            "n_deg_groups": n_deg_groups,
            "n_deconv_cell_types": n_deconv_cell_types
        }
        
        # Check if exists
        existing = db.query(models.H5ADAnalysisFeatures).filter(
            models.H5ADAnalysisFeatures.dataset_id == dataset_id
        ).first()
        
        if existing:
            # Update
            for key, value in features_data.items():
                setattr(existing, key, value)
            db.commit()
            db.refresh(existing)
            return existing
        else:
            # Create
            db_features = models.H5ADAnalysisFeatures(**features_data)
            db.add(db_features)
            db.commit()
            db.refresh(db_features)
            return db_features
            
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error syncing analysis features: {e}")
        raise HTTPException(status_code=500, detail=str(e))

