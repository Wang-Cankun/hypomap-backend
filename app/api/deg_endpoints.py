"""
API endpoints for Differential Expression Gene (DEG) Analysis
"""
from fastapi import APIRouter, HTTPException, Body
import logging

from app.services.deg_service import DEGService
from app.schemas import (
    DEGBetweenDatasetsRequest,
    DEGWithinDatasetRequest,
    DEGAnalysisResult,
    AtlasDEGRequest,
    AtlasDEGResult,
    PreviewCellCountsRequest,
    PreviewCellCountsResponse
)

logger = logging.getLogger(__name__)

router = APIRouter()
deg_service = DEGService()


@router.post("/between-datasets", response_model=DEGAnalysisResult)
async def analyze_deg_between_datasets(request: DEGBetweenDatasetsRequest):
    """
    Perform DEG analysis between two datasets
    
    **Use Cases:**
    - Compare different conditions (e.g., AD093044 vs AD093045)
    - Filter by cell type (e.g., only "Neuron" cells)
    - Filter by metadata (e.g., specific cluster or condition)
    
    **Example:**
    ```json
    {
      "dataset_id1": "AD093044",
      "dataset_id2": "AD093045",
      "cell_type": "Neuron",
      "min_pct": 0.1,
      "logfc_threshold": 0.25,
      "p_value_threshold": 0.05,
      "top_n": 100
    }
    ```
    
    **Parameters:**
    - **dataset_id1**: First dataset (group 1)
    - **dataset_id2**: Second dataset (group 2)
    - **cell_type**: Optional, filter by cell type
    - **metadata_filters1/2**: Optional, additional metadata filters
    - **min_pct**: Minimum fraction of cells expressing gene (default: 0.1)
    - **logfc_threshold**: Minimum log2 fold change (default: 0.25)
    - **p_value_threshold**: Adjusted p-value threshold (default: 0.05)
    - **top_n**: Return only top N genes (optional)
    
    **Returns:**
    - List of genes with fold changes, p-values, and significance
    - Summary statistics
    """
    try:
        result = deg_service.analyze_deg_between_datasets(
            dataset_id1=request.dataset_id1,
            dataset_id2=request.dataset_id2,
            cell_type=request.cell_type,
            cell_type2=request.cell_type2,
            cell_types_group1=request.cell_types_group1,
            cell_types_group2=request.cell_types_group2,
            metadata_filters1=request.metadata_filters1,
            metadata_filters2=request.metadata_filters2,
            min_pct=request.min_pct,
            logfc_threshold=request.logfc_threshold,
            p_value_threshold=request.p_value_threshold,
            top_n=request.top_n
        )
        return result
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error in DEG analysis: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"DEG analysis failed: {str(e)}")


@router.post("/within-dataset", response_model=DEGAnalysisResult)
async def analyze_deg_within_dataset(request: DEGWithinDatasetRequest):
    """
    Perform DEG analysis within a single dataset between two groups
    
    **Use Cases:**
    - Compare cell types within same dataset
    - Compare conditions within same dataset
    - Compare clusters within same dataset
    
    **Example:**
    ```json
    {
      "dataset_id": "AD093044",
      "group1_filters": {"cell_type": "Neuron"},
      "group2_filters": {"cell_type": "Astrocyte"},
      "min_pct": 0.1,
      "logfc_threshold": 0.5,
      "p_value_threshold": 0.01,
      "top_n": 50
    }
    ```
    
    **Parameters:**
    - **dataset_id**: Dataset to analyze
    - **group1_filters**: Metadata filters for group 1 (e.g., {"cell_type": "Neuron"})
    - **group2_filters**: Metadata filters for group 2 (e.g., {"cell_type": "Astrocyte"})
    - **min_pct**: Minimum fraction of cells expressing gene
    - **logfc_threshold**: Minimum log2 fold change
    - **p_value_threshold**: Adjusted p-value threshold
    - **top_n**: Return only top N genes
    
    **Returns:**
    - List of genes with fold changes, p-values, and significance
    - Summary statistics
    """
    try:
        result = deg_service.analyze_deg_within_dataset(
            dataset_id=request.dataset_id,
            group1_filters=request.group1_filters,
            group2_filters=request.group2_filters,
            min_pct=request.min_pct,
            logfc_threshold=request.logfc_threshold,
            p_value_threshold=request.p_value_threshold,
            top_n=request.top_n
        )
        return result
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error in DEG analysis: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"DEG analysis failed: {str(e)}")


@router.get("/cell-types/{dataset_id}")
async def get_available_cell_types(dataset_id: str):
    """
    Get available cell types in a dataset for filtering
    
    **Example:**
    ```
    GET /deg/cell-types/AD093044
    ```
    
    **Returns:**
    List of unique cell types available in the dataset
    """
    try:
        import anndata
        from pathlib import Path
        
        h5ad_path = Path("h5ad/raw") / f"{dataset_id}.h5ad"
        if not h5ad_path.exists():
            raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
        
        adata = anndata.read_h5ad(h5ad_path)
        
        if 'cell_type' not in adata.obs.columns:
            return {"cell_types": [], "message": "No 'cell_type' column found"}
        
        cell_types = adata.obs['cell_type'].unique().tolist()
        cell_type_counts = adata.obs['cell_type'].value_counts().to_dict()
        
        return {
            "dataset_id": dataset_id,
            "cell_types": sorted(cell_types),
            "counts": cell_type_counts
        }
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error getting cell types: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/metadata-columns/{dataset_id}")
async def get_metadata_columns(dataset_id: str):
    """
    Get available metadata columns for filtering
    
    **Example:**
    ```
    GET /deg/metadata-columns/AD093044
    ```
    
    **Returns:**
    List of metadata columns and their unique values
    """
    try:
        import anndata
        from pathlib import Path
        
        h5ad_path = Path("h5ad/raw") / f"{dataset_id}.h5ad"
        if not h5ad_path.exists():
            raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
        
        adata = anndata.read_h5ad(h5ad_path)
        
        metadata = {}
        for col in adata.obs.columns:
            unique_values = adata.obs[col].unique().tolist()
            # Limit to reasonable number of unique values
            if len(unique_values) <= 100:
                metadata[col] = {
                    "unique_values": sorted([str(v) for v in unique_values]),
                    "n_unique": len(unique_values),
                    "counts": adata.obs[col].value_counts().head(20).to_dict()
                }
        
        return {
            "dataset_id": dataset_id,
            "metadata_columns": metadata
        }
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")
    except Exception as e:
        logger.error(f"Error getting metadata columns: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/preview-cell-counts", response_model=PreviewCellCountsResponse)
async def preview_cell_counts(request: PreviewCellCountsRequest):
    """
    Preview cell counts for filter criteria without running full DEG analysis.
    
    This endpoint allows the frontend to show users how many cells match their 
    filter criteria before running the full DEG analysis.
    
    **Example:**
    ```json
    {
      "dataset_id": "human_subset",
      "group1_filters": [
        {"column": "cell_type", "values": ["Astrocyte"]},
        {"column": "supercluster_name", "values": ["Astrocyte"]}
      ],
      "group2_filters": [
        {"column": "cell_type", "values": ["Neuron"]},
        {"column": "supercluster_name", "values": ["Neuron"]}
      ]
    }
    ```
    
    **Returns:**
    - Cell counts for each group
    - Overlap count (should be 0)
    - Warnings if groups are too small or overlapping
    """
    try:
        # Convert Pydantic models to dicts
        group1_filters = [f.dict() for f in request.group1_filters]
        group2_filters = [f.dict() for f in request.group2_filters]
        
        result = deg_service.preview_cell_counts(
            dataset_id=request.dataset_id,
            group1_filters=group1_filters,
            group2_filters=group2_filters
        )
        return PreviewCellCountsResponse(**result)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {request.dataset_id} not found")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error previewing cell counts: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/atlas-compare", response_model=AtlasDEGResult)
async def atlas_compare(request: AtlasDEGRequest):
    """
    Perform advanced DEG analysis with multi-filter groups for Atlas datasets.
    
    This endpoint enables complex comparisons using multiple metadata filters:
    - Astrocytes in Males vs Astrocytes in Females
    - Neurons in AD patients from Hippocampus vs Neurons in Controls from Hippocampus
    - Multiple cell types with specific conditions
    
    **Filter Logic:**
    - Multiple filters within a group use AND logic
    - Multiple values within a filter use OR logic
    
    **Example:**
    ```json
    {
      "dataset_id": "human_subset",
      "group1_filters": [
        {"column": "cell_type", "values": ["Astrocyte", "Oligodendrocyte"]},
        {"column": "supercluster_name", "values": ["Astrocyte"]}
      ],
      "group2_filters": [
        {"column": "cell_type", "values": ["Neuron"]},
        {"column": "supercluster_name", "values": ["Neuron"]}
      ],
      "logfc_threshold": 0.25,
      "p_value_threshold": 0.05,
      "min_pct": 0.1
    }
    ```
    
    **Parameters:**
    - **dataset_id**: Atlas dataset identifier (e.g., "human_subset")
    - **group1_filters**: List of filters for group 1
    - **group2_filters**: List of filters for group 2
    - **logfc_threshold**: Minimum log2 fold change (default: 0.25)
    - **p_value_threshold**: Adjusted p-value threshold (default: 0.05)
    - **min_pct**: Minimum fraction of cells expressing gene (default: 0.1)
    - **top_n**: Optional limit on number of genes returned
    
    **Returns:**
    - Summary statistics (cell counts, significant genes)
    - Group descriptions with applied filters
    - List of DEG results sorted by significance
    
    **Error Codes:**
    - `EMPTY_GROUP`: No cells match filters for a group
    - `INVALID_COLUMN`: Column not found in dataset
    - `OVERLAPPING_GROUPS`: Groups have overlapping cells
    """
    try:
        # Convert Pydantic models to dicts
        group1_filters = [f.dict() for f in request.group1_filters]
        group2_filters = [f.dict() for f in request.group2_filters]
        
        result = deg_service.atlas_compare(
            dataset_id=request.dataset_id,
            group1_filters=group1_filters,
            group2_filters=group2_filters,
            logfc_threshold=request.logfc_threshold,
            p_value_threshold=request.p_value_threshold,
            min_pct=request.min_pct,
            top_n=request.top_n
        )
        return AtlasDEGResult(**result)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset {request.dataset_id} not found")
    except ValueError as e:
        error_msg = str(e)
        # Provide specific error codes for common cases
        if "No cells match" in error_msg:
            if "Group 1" in error_msg:
                raise HTTPException(status_code=400, detail={"detail": error_msg, "error_code": "EMPTY_GROUP", "group": 1})
            else:
                raise HTTPException(status_code=400, detail={"detail": error_msg, "error_code": "EMPTY_GROUP", "group": 2})
        elif "overlapping" in error_msg.lower():
            # Extract overlap count if present
            import re
            match = re.search(r'(\d+) cells', error_msg)
            overlap_count = int(match.group(1)) if match else 0
            raise HTTPException(status_code=400, detail={"detail": error_msg, "error_code": "OVERLAPPING_GROUPS", "overlap_count": overlap_count})
        elif "not found in dataset" in error_msg:
            # Extract column name if present
            match = re.search(r"'([^']+)' not found", error_msg)
            column = match.group(1) if match else None
            raise HTTPException(status_code=400, detail={"detail": error_msg, "error_code": "INVALID_COLUMN", "column": column})
        else:
            raise HTTPException(status_code=400, detail=error_msg)
    except Exception as e:
        logger.error(f"Error in atlas DEG comparison: {e}")
        raise HTTPException(status_code=500, detail=str(e))

