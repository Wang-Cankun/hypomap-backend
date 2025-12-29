"""
Upload endpoints for HypoMap - handles h5ad, regulon, and CellChat file uploads
"""
import os
import shutil
from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse
from typing import Optional
import json

router = APIRouter()

# Get data directories from environment or use defaults
H5AD_RAW_DIR = os.getenv("H5AD_RAW_DIR", "h5ad/raw")
H5AD_PRECOMPUTED_DIR = os.getenv("H5AD_PRECOMPUTED_DIR", "h5ad/precomputed")

# Ensure directories exist
os.makedirs(H5AD_RAW_DIR, exist_ok=True)
os.makedirs(H5AD_PRECOMPUTED_DIR, exist_ok=True)


def run_preprocessing(dataset_id: str):
    """Background task to preprocess uploaded h5ad file"""
    try:
        import subprocess
        result = subprocess.run(
            ["python", "scripts/preprocess_h5ad.py", dataset_id],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        )

        # Update status file
        status_file = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id, "preprocess_status.json")
        os.makedirs(os.path.dirname(status_file), exist_ok=True)

        status = {
            "status": "completed" if result.returncode == 0 else "failed",
            "message": result.stdout if result.returncode == 0 else result.stderr
        }

        with open(status_file, "w") as f:
            json.dump(status, f)

    except Exception as e:
        status_file = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id, "preprocess_status.json")
        os.makedirs(os.path.dirname(status_file), exist_ok=True)
        with open(status_file, "w") as f:
            json.dump({"status": "failed", "message": str(e)}, f)


@router.post("/upload/h5ad")
async def upload_h5ad(
    file: UploadFile = File(...),
    dataset_id: Optional[str] = None
):
    """Upload an h5ad file"""
    if not file.filename.endswith('.h5ad'):
        raise HTTPException(status_code=400, detail="File must be .h5ad format")

    # Generate dataset_id from filename if not provided
    if not dataset_id:
        dataset_id = os.path.splitext(file.filename)[0]

    # Save file
    file_path = os.path.join(H5AD_RAW_DIR, f"{dataset_id}.h5ad")

    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to save file: {str(e)}")

    return {
        "dataset_id": dataset_id,
        "filename": file.filename,
        "path": file_path,
        "message": "File uploaded successfully"
    }


@router.post("/upload/regulon/{dataset_id}")
async def upload_regulon(
    dataset_id: str,
    file: UploadFile = File(...),
    cluster: Optional[str] = None
):
    """Upload regulon CSV file for a dataset"""
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="File must be .csv format")

    # Create regulon directory
    regulon_dir = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id, "regulon")
    os.makedirs(regulon_dir, exist_ok=True)

    # Determine filename
    if cluster:
        filename = f"{cluster}_regulon_network.csv"
    else:
        filename = file.filename

    file_path = os.path.join(regulon_dir, filename)

    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to save file: {str(e)}")

    return {
        "dataset_id": dataset_id,
        "filename": filename,
        "path": file_path,
        "message": "Regulon file uploaded successfully"
    }


@router.post("/upload/cellchat/{dataset_id}")
async def upload_cellchat(
    dataset_id: str,
    file: UploadFile = File(...)
):
    """Upload CellChat communications CSV file"""
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="File must be .csv format")

    # Create CCC directory
    ccc_dir = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id, "ccc")
    os.makedirs(ccc_dir, exist_ok=True)

    file_path = os.path.join(ccc_dir, "communications.csv")

    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to save file: {str(e)}")

    return {
        "dataset_id": dataset_id,
        "filename": "communications.csv",
        "path": file_path,
        "message": "CellChat file uploaded successfully"
    }


@router.post("/upload/preprocess/{dataset_id}")
async def start_preprocessing(
    dataset_id: str,
    background_tasks: BackgroundTasks
):
    """Start preprocessing for an uploaded dataset"""
    # Check if h5ad file exists
    h5ad_path = os.path.join(H5AD_RAW_DIR, f"{dataset_id}.h5ad")
    if not os.path.exists(h5ad_path):
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")

    # Set initial status
    status_file = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id, "preprocess_status.json")
    os.makedirs(os.path.dirname(status_file), exist_ok=True)
    with open(status_file, "w") as f:
        json.dump({"status": "processing", "message": "Preprocessing started"}, f)

    # Start background task
    background_tasks.add_task(run_preprocessing, dataset_id)

    return {
        "dataset_id": dataset_id,
        "status": "processing",
        "message": "Preprocessing started in background"
    }


@router.get("/upload/preprocess/status/{dataset_id}")
async def get_preprocessing_status(dataset_id: str):
    """Get preprocessing status for a dataset"""
    status_file = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id, "preprocess_status.json")

    if not os.path.exists(status_file):
        return {"status": "unknown", "message": "No preprocessing status found"}

    with open(status_file, "r") as f:
        return json.load(f)


@router.get("/upload/datasets")
async def list_uploaded_datasets():
    """List all uploaded datasets with their status"""
    datasets = []

    # Check raw h5ad files
    if os.path.exists(H5AD_RAW_DIR):
        for filename in os.listdir(H5AD_RAW_DIR):
            if filename.endswith('.h5ad'):
                dataset_id = os.path.splitext(filename)[0]
                precomputed_dir = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id)

                # Check what's available
                has_regulon = os.path.exists(os.path.join(precomputed_dir, "regulon"))
                has_ccc = os.path.exists(os.path.join(precomputed_dir, "ccc"))
                info_file = os.path.join(precomputed_dir, "info.json")

                info = {}
                if os.path.exists(info_file):
                    with open(info_file, "r") as f:
                        info = json.load(f)

                datasets.append({
                    "dataset_id": dataset_id,
                    "h5ad_exists": True,
                    "preprocessed": os.path.exists(precomputed_dir),
                    "has_regulon": has_regulon,
                    "has_ccc": has_ccc,
                    "n_cells": info.get("n_cells"),
                    "n_genes": info.get("n_genes")
                })

    return {"datasets": datasets, "total": len(datasets)}


@router.delete("/upload/{dataset_id}")
async def delete_dataset(dataset_id: str):
    """Delete an uploaded dataset and its precomputed data"""
    h5ad_path = os.path.join(H5AD_RAW_DIR, f"{dataset_id}.h5ad")
    precomputed_dir = os.path.join(H5AD_PRECOMPUTED_DIR, dataset_id)

    deleted = []

    if os.path.exists(h5ad_path):
        os.remove(h5ad_path)
        deleted.append("h5ad")

    if os.path.exists(precomputed_dir):
        shutil.rmtree(precomputed_dir)
        deleted.append("precomputed")

    if not deleted:
        raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found")

    return {
        "dataset_id": dataset_id,
        "deleted": deleted,
        "message": "Dataset deleted successfully"
    }
