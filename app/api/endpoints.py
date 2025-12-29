from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from datetime import datetime

from app.database import get_db
from app import crud, schemas, models

router = APIRouter()

# Hello World endpoints
@router.get("/", response_model=schemas.HelloWorld)
def hello_world():
    """Get hello world message"""
    return schemas.HelloWorld(
        message="Hello, World!",
        timestamp=datetime.now()
    )

@router.post("/hello", response_model=schemas.HelloWorld)
def create_hello_world(hello_data: schemas.HelloWorld):
    """Create a hello world message"""
    return hello_data

# Message endpoints
@router.post("/messages/", response_model=schemas.Message, status_code=status.HTTP_201_CREATED)
def create_message(message: schemas.MessageCreate, db: Session = Depends(get_db)):
    """Create a new message"""
    return crud.create_message(db=db, message=message)

@router.get("/messages/", response_model=List[schemas.Message])
def read_messages(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """Get all messages with pagination"""
    messages = crud.get_messages(db, skip=skip, limit=limit)
    return messages

@router.get("/messages/{message_id}", response_model=schemas.Message)
def read_message(message_id: int, db: Session = Depends(get_db)):
    """Get a specific message by ID"""
    db_message = crud.get_message(db, message_id=message_id)
    if db_message is None:
        raise HTTPException(status_code=404, detail="Message not found")
    return db_message

@router.put("/messages/{message_id}", response_model=schemas.Message)
def update_message(message_id: int, message: schemas.MessageUpdate, db: Session = Depends(get_db)):
    """Update a message"""
    db_message = crud.update_message(db, message_id=message_id, message=message)
    if db_message is None:
        raise HTTPException(status_code=404, detail="Message not found")
    return db_message

@router.delete("/messages/{message_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_message(message_id: int, db: Session = Depends(get_db)):
    """Delete a message"""
    success = crud.delete_message(db, message_id=message_id)
    if not success:
        raise HTTPException(status_code=404, detail="Message not found")
    return None

# ScRNA-seq Dataset endpoints
@router.post("/datasets/", response_model=schemas.ScRNAseqDataset, status_code=status.HTTP_201_CREATED)
def create_dataset(dataset: schemas.ScRNAseqDatasetCreate, db: Session = Depends(get_db)):
    """Create a new scRNA-seq dataset"""
    return crud.create_scrna_dataset(db=db, dataset=dataset)

@router.get("/datasets/", response_model=List[schemas.ScRNAseqDataset])
def read_datasets(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """Get all scRNA-seq datasets with pagination"""
    datasets = crud.get_scrna_datasets(db, skip=skip, limit=limit)
    return datasets

@router.get("/datasets/{dataset_id}")
def read_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """
    Get a specific dataset by ID - checks scrna_datasets, spatial_datasets, and h5ad datasets
    Returns the appropriate dataset type or 404 if not found
    """
    # First check scrna_datasets
    db_dataset = crud.get_scrna_dataset(db, dataset_id=dataset_id)
    if db_dataset:
        return schemas.ScRNAseqDataset.model_validate(db_dataset)
    
    # Then check spatial_datasets
    spatial_dataset = crud.get_spatial_dataset(db, dataset_id=dataset_id)
    if spatial_dataset:
        return schemas.SpatialDataset.model_validate(spatial_dataset)
    
    # Finally check if it's an h5ad dataset (check if precomputed directory exists)
    from pathlib import Path
    h5ad_precomputed_path = Path(f"h5ad/precomputed/{dataset_id}")
    if h5ad_precomputed_path.exists() and (h5ad_precomputed_path / "info.json").exists():
        # Return a minimal dataset info for h5ad datasets
        import json
        try:
            with open(h5ad_precomputed_path / "info.json", 'r') as f:
                info = json.load(f)
            return {
                "dataset_id": dataset_id,
                "type": "h5ad",
                "n_cells": info.get("n_cells"),
                "n_genes": info.get("n_genes"),
                "has_umap": "umap" in info.get("embeddings", []),
                "has_tsne": "tsne" in info.get("embeddings", []),
                "has_pca": "pca" in info.get("embeddings", []),
                "metadata_columns": info.get("metadata_columns", []),
                "available_embeddings": info.get("embeddings", [])
            }
        except Exception as e:
            # If we can't read info.json, still return basic info
            return {
                "dataset_id": dataset_id,
                "type": "h5ad",
                "available": True
            }
    
    raise HTTPException(status_code=404, detail=f"Dataset {dataset_id} not found in scrna_datasets, spatial_datasets, or h5ad datasets")

@router.put("/datasets/{dataset_id}", response_model=schemas.ScRNAseqDataset)
def update_dataset(dataset_id: str, dataset: schemas.ScRNAseqDatasetUpdate, db: Session = Depends(get_db)):
    """Update a scRNA-seq dataset"""
    db_dataset = crud.update_scrna_dataset(db, dataset_id=dataset_id, dataset=dataset)
    if db_dataset is None:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return db_dataset

@router.delete("/datasets/{dataset_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """Delete a scRNA-seq dataset"""
    success = crud.delete_scrna_dataset(db, dataset_id=dataset_id)
    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return None

# Search endpoints for scRNA-seq datasets
@router.get("/datasets/search/disease/{disease}", response_model=List[schemas.ScRNAseqDataset])
def search_datasets_by_disease(
    disease: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search scRNA-seq datasets by disease"""
    datasets = crud.search_scrna_datasets_by_disease(db, disease=disease, skip=skip, limit=limit)
    return datasets

@router.get("/datasets/search/species/{species}", response_model=List[schemas.ScRNAseqDataset])
def search_datasets_by_species(
    species: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search scRNA-seq datasets by species"""
    datasets = crud.search_scrna_datasets_by_species(db, species=species, skip=skip, limit=limit)
    return datasets

@router.get("/datasets/search/brain-region/{brain_region}", response_model=List[schemas.ScRNAseqDataset])
def search_datasets_by_brain_region(
    brain_region: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search scRNA-seq datasets by brain region"""
    datasets = crud.search_scrna_datasets_by_brain_region(db, brain_region=brain_region, skip=skip, limit=limit)
    return datasets

# Statistics endpoint
@router.get("/datasets/stats/")
def get_dataset_stats(db: Session = Depends(get_db)):
    """Get statistics about the scRNA-seq datasets"""
    total_datasets = db.query(models.ScRNAseqDataset).count()
    
    # Get unique values for key fields
    diseases = db.query(models.ScRNAseqDataset.disease).distinct().all()
    species = db.query(models.ScRNAseqDataset.species).distinct().all()
    brain_regions = db.query(models.ScRNAseqDataset.brain_region).distinct().all()
    
    return {
        "total_datasets": total_datasets,
        "unique_diseases": [d[0] for d in diseases if d[0]],
        "unique_species": [s[0] for s in species if s[0]],
        "unique_brain_regions": [b[0] for b in brain_regions if b[0]]
    }

# Papers endpoints
@router.post("/papers/", response_model=schemas.ScRNAPaper, status_code=status.HTTP_201_CREATED)
def create_paper(paper: schemas.ScRNAPaperCreate, db: Session = Depends(get_db)):
    return crud.create_paper(db, paper)

@router.get("/papers/", response_model=List[schemas.ScRNAPaper])
def read_papers(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.get_papers(db, skip=skip, limit=limit)

@router.get("/papers/{paper_id}", response_model=schemas.ScRNAPaper)
def read_paper(paper_id: str, db: Session = Depends(get_db)):
    db_paper = crud.get_paper(db, paper_id)
    if db_paper is None:
        raise HTTPException(status_code=404, detail="Paper not found")
    return db_paper

@router.put("/papers/{paper_id}", response_model=schemas.ScRNAPaper)
def update_paper(paper_id: str, paper: schemas.ScRNAPaperUpdate, db: Session = Depends(get_db)):
    db_paper = crud.update_paper(db, paper_id, paper)
    if db_paper is None:
        raise HTTPException(status_code=404, detail="Paper not found")
    return db_paper

@router.delete("/papers/{paper_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_paper(paper_id: str, db: Session = Depends(get_db)):
    success = crud.delete_paper(db, paper_id)
    if not success:
        raise HTTPException(status_code=404, detail="Paper not found")
    return None

@router.get("/papers/search/disease/{disease}", response_model=List[schemas.ScRNAPaper])
def search_papers_disease(disease: str, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.search_papers_by_disease(db, disease, skip=skip, limit=limit)

@router.get("/papers/search/species/{species}", response_model=List[schemas.ScRNAPaper])
def search_papers_species(species: str, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.search_papers_by_species(db, species, skip=skip, limit=limit)

@router.get("/papers/search/brain-region/{brain_region}", response_model=List[schemas.ScRNAPaper])
def search_papers_brain_region(brain_region: str, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.search_papers_by_brain_region(db, brain_region, skip=skip, limit=limit)

# Spatial Dataset endpoints
@router.post("/spatial-datasets/", response_model=schemas.SpatialDataset, status_code=status.HTTP_201_CREATED)
def create_spatial_dataset(dataset: schemas.SpatialDatasetCreate, db: Session = Depends(get_db)):
    """Create a new spatial dataset"""
    return crud.create_spatial_dataset(db=db, dataset=dataset)

@router.get("/spatial-datasets/", response_model=List[schemas.SpatialDataset])
def read_spatial_datasets(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """Get all spatial datasets with pagination"""
    datasets = crud.get_spatial_datasets(db, skip=skip, limit=limit)
    return datasets

@router.get("/spatial-datasets/{dataset_id}", response_model=schemas.SpatialDataset)
def read_spatial_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """Get a specific spatial dataset by dataset ID"""
    db_dataset = crud.get_spatial_dataset(db, dataset_id=dataset_id)
    if db_dataset is None:
        raise HTTPException(status_code=404, detail="Spatial dataset not found")
    return db_dataset

@router.put("/spatial-datasets/{dataset_id}", response_model=schemas.SpatialDataset)
def update_spatial_dataset(dataset_id: str, dataset: schemas.SpatialDatasetUpdate, db: Session = Depends(get_db)):
    """Update a spatial dataset"""
    db_dataset = crud.update_spatial_dataset(db, dataset_id=dataset_id, dataset=dataset)
    if db_dataset is None:
        raise HTTPException(status_code=404, detail="Spatial dataset not found")
    return db_dataset

@router.delete("/spatial-datasets/{dataset_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_spatial_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """Delete a spatial dataset"""
    success = crud.delete_spatial_dataset(db, dataset_id=dataset_id)
    if not success:
        raise HTTPException(status_code=404, detail="Spatial dataset not found")
    return None

# Search endpoints for spatial datasets
@router.get("/spatial-datasets/search/disease/{disease}", response_model=List[schemas.SpatialDataset])
def search_spatial_datasets_by_disease(
    disease: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial datasets by disease"""
    datasets = crud.search_spatial_datasets_by_disease(db, disease=disease, skip=skip, limit=limit)
    return datasets

@router.get("/spatial-datasets/search/species/{species}", response_model=List[schemas.SpatialDataset])
def search_spatial_datasets_by_species(
    species: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial datasets by species"""
    datasets = crud.search_spatial_datasets_by_species(db, species=species, skip=skip, limit=limit)
    return datasets

@router.get("/spatial-datasets/search/brain-region/{brain_region}", response_model=List[schemas.SpatialDataset])
def search_spatial_datasets_by_brain_region(
    brain_region: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial datasets by brain region"""
    datasets = crud.search_spatial_datasets_by_brain_region(db, brain_region=brain_region, skip=skip, limit=limit)
    return datasets

@router.get("/spatial-datasets/search/methodology/{methodology}", response_model=List[schemas.SpatialDataset])
def search_spatial_datasets_by_methodology(
    methodology: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial datasets by methodology (e.g., Visium, MERFISH, Xenium)"""
    datasets = crud.search_spatial_datasets_by_methodology(db, methodology=methodology, skip=skip, limit=limit)
    return datasets

@router.get("/spatial-datasets/search/study/{data_id}", response_model=List[schemas.SpatialDataset])
def search_spatial_datasets_by_study(
    data_id: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Get all datasets from a specific study/project (data_id like ST024)"""
    datasets = crud.search_spatial_datasets_by_data_id(db, data_id=data_id, skip=skip, limit=limit)
    return datasets

# Statistics endpoint for spatial datasets
@router.get("/spatial-datasets/stats/")
def get_spatial_dataset_stats(db: Session = Depends(get_db)):
    """Get statistics about the spatial datasets"""
    total_datasets = db.query(models.SpatialDataset).count()
    
    # Get unique values for key fields
    diseases = db.query(models.SpatialDataset.disease).distinct().all()
    species = db.query(models.SpatialDataset.species).distinct().all()
    brain_regions = db.query(models.SpatialDataset.brain_region).distinct().all()
    methodologies = db.query(models.SpatialDataset.methodology).distinct().all()
    
    return {
        "total_datasets": total_datasets,
        "unique_diseases": [d[0] for d in diseases if d[0]],
        "unique_species": [s[0] for s in species if s[0]],
        "unique_brain_regions": [b[0] for b in brain_regions if b[0]],
        "unique_methodologies": [m[0] for m in methodologies if m[0]]
    }

# Spatial Papers endpoints
@router.post("/spatial-papers/", response_model=schemas.SpatialPaper, status_code=status.HTTP_201_CREATED)
def create_spatial_paper(paper: schemas.SpatialPaperCreate, db: Session = Depends(get_db)):
    """Create a new spatial paper"""
    return crud.create_spatial_paper(db=db, paper=paper)

@router.get("/spatial-papers/", response_model=List[schemas.SpatialPaper])
def read_spatial_papers(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """Get all spatial papers with pagination"""
    papers = crud.get_spatial_papers(db, skip=skip, limit=limit)
    return papers

@router.get("/spatial-papers/{data_id}", response_model=schemas.SpatialPaper)
def read_spatial_paper(data_id: str, db: Session = Depends(get_db)):
    """Get a specific spatial paper by data ID"""
    db_paper = crud.get_spatial_paper(db, data_id=data_id)
    if db_paper is None:
        raise HTTPException(status_code=404, detail="Spatial paper not found")
    return db_paper

@router.put("/spatial-papers/{data_id}", response_model=schemas.SpatialPaper)
def update_spatial_paper(data_id: str, paper: schemas.SpatialPaperUpdate, db: Session = Depends(get_db)):
    """Update a spatial paper"""
    db_paper = crud.update_spatial_paper(db, data_id=data_id, paper=paper)
    if db_paper is None:
        raise HTTPException(status_code=404, detail="Spatial paper not found")
    return db_paper

@router.delete("/spatial-papers/{data_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_spatial_paper(data_id: str, db: Session = Depends(get_db)):
    """Delete a spatial paper"""
    success = crud.delete_spatial_paper(db, data_id=data_id)
    if not success:
        raise HTTPException(status_code=404, detail="Spatial paper not found")
    return None

# Search endpoints for spatial papers
@router.get("/spatial-papers/search/disease/{disease}", response_model=List[schemas.SpatialPaper])
def search_spatial_papers_by_disease(
    disease: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial papers by disease"""
    papers = crud.search_spatial_papers_by_disease(db, disease=disease, skip=skip, limit=limit)
    return papers

@router.get("/spatial-papers/search/species/{species}", response_model=List[schemas.SpatialPaper])
def search_spatial_papers_by_species(
    species: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial papers by species"""
    papers = crud.search_spatial_papers_by_species(db, species=species, skip=skip, limit=limit)
    return papers

@router.get("/spatial-papers/search/brain-region/{brain_region}", response_model=List[schemas.SpatialPaper])
def search_spatial_papers_by_brain_region(
    brain_region: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial papers by brain region"""
    papers = crud.search_spatial_papers_by_brain_region(db, brain_region=brain_region, skip=skip, limit=limit)
    return papers

@router.get("/spatial-papers/search/methodology/{methodology}", response_model=List[schemas.SpatialPaper])
def search_spatial_papers_by_methodology(
    methodology: str, 
    skip: int = 0, 
    limit: int = 100, 
    db: Session = Depends(get_db)
):
    """Search spatial papers by methodology (e.g., Visium, MERFISH, Xenium)"""
    papers = crud.search_spatial_papers_by_methodology(db, methodology=methodology, skip=skip, limit=limit)
    return papers

# Statistics endpoint for spatial papers
@router.get("/spatial-papers/stats/")
def get_spatial_paper_stats(db: Session = Depends(get_db)):
    """Get statistics about the spatial papers"""
    total_papers = db.query(models.SpatialPaper).count()
    
    # Get unique values for key fields
    diseases = db.query(models.SpatialPaper.disease).distinct().all()
    species = db.query(models.SpatialPaper.species).distinct().all()
    brain_regions = db.query(models.SpatialPaper.brain_region).distinct().all()
    methodologies = db.query(models.SpatialPaper.methodology).distinct().all()
    journals = db.query(models.SpatialPaper.journal).distinct().all()
    
    return {
        "total_papers": total_papers,
        "unique_diseases": [d[0] for d in diseases if d[0]],
        "unique_species": [s[0] for s in species if s[0]],
        "unique_brain_regions": [b[0] for b in brain_regions if b[0]],
        "unique_methodologies": [m[0] for m in methodologies if m[0]],
        "unique_journals": [j[0] for j in journals if j[0]]
    } 