from sqlalchemy.orm import Session
from app import models, schemas
from typing import List, Optional

# Create a new message
def create_message(db: Session, message: schemas.MessageCreate) -> models.Message:
    db_message = models.Message(**message.model_dump())
    db.add(db_message)
    db.commit()
    db.refresh(db_message)
    return db_message

# Get a message by ID
def get_message(db: Session, message_id: int) -> Optional[models.Message]:
    return db.query(models.Message).filter(models.Message.id == message_id).first()

# Get all messages
def get_messages(db: Session, skip: int = 0, limit: int = 100) -> List[models.Message]:
    return db.query(models.Message).offset(skip).limit(limit).all()

# Update a message
def update_message(db: Session, message_id: int, message: schemas.MessageUpdate) -> Optional[models.Message]:
    db_message = get_message(db, message_id)
    if db_message:
        update_data = message.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(db_message, field, value)
        db.commit()
        db.refresh(db_message)
    return db_message

# Delete a message
def delete_message(db: Session, message_id: int) -> bool:
    db_message = get_message(db, message_id)
    if db_message:
        db.delete(db_message)
        db.commit()
        return True
    return False

# ScRNA-seq Dataset CRUD operations
def create_scrna_dataset(db: Session, dataset: schemas.ScRNAseqDatasetCreate) -> models.ScRNAseqDataset:
    db_dataset = models.ScRNAseqDataset(**dataset.model_dump())
    db.add(db_dataset)
    db.commit()
    db.refresh(db_dataset)
    return db_dataset

def get_scrna_dataset(db: Session, dataset_id: str) -> Optional[models.ScRNAseqDataset]:
    return db.query(models.ScRNAseqDataset).filter(models.ScRNAseqDataset.dataset_id == dataset_id).first()

def get_scrna_datasets(db: Session, skip: int = 0, limit: int = 100) -> List[models.ScRNAseqDataset]:
    return db.query(models.ScRNAseqDataset).offset(skip).limit(limit).all()

def update_scrna_dataset(db: Session, dataset_id: str, dataset: schemas.ScRNAseqDatasetUpdate) -> Optional[models.ScRNAseqDataset]:
    db_dataset = get_scrna_dataset(db, dataset_id)
    if db_dataset:
        update_data = dataset.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(db_dataset, field, value)
        db.commit()
        db.refresh(db_dataset)
    return db_dataset

def delete_scrna_dataset(db: Session, dataset_id: str) -> bool:
    db_dataset = get_scrna_dataset(db, dataset_id)
    if db_dataset:
        db.delete(db_dataset)
        db.commit()
        return True
    return False

# Search functions for scRNA-seq datasets
def search_scrna_datasets_by_disease(db: Session, disease: str, skip: int = 0, limit: int = 100) -> List[models.ScRNAseqDataset]:
    return db.query(models.ScRNAseqDataset).filter(
        models.ScRNAseqDataset.disease.ilike(f"%{disease}%")
    ).offset(skip).limit(limit).all()

def search_scrna_datasets_by_species(db: Session, species: str, skip: int = 0, limit: int = 100) -> List[models.ScRNAseqDataset]:
    return db.query(models.ScRNAseqDataset).filter(
        models.ScRNAseqDataset.species.ilike(f"%{species}%")
    ).offset(skip).limit(limit).all()

def search_scrna_datasets_by_brain_region(db: Session, brain_region: str, skip: int = 0, limit: int = 100) -> List[models.ScRNAseqDataset]:
    return db.query(models.ScRNAseqDataset).filter(
        models.ScRNAseqDataset.brain_region.ilike(f"%{brain_region}%")
    ).offset(skip).limit(limit).all()

# Papers CRUD operations
def create_paper(db: Session, paper: schemas.ScRNAPaperCreate) -> models.ScRNAPaper:
    db_paper = models.ScRNAPaper(**paper.model_dump())
    db.add(db_paper)
    db.commit()
    db.refresh(db_paper)
    return db_paper

def get_paper(db: Session, paper_id: str) -> Optional[models.ScRNAPaper]:
    return db.query(models.ScRNAPaper).filter(models.ScRNAPaper.paper_id == paper_id).first()

def get_papers(db: Session, skip: int = 0, limit: int = 100) -> List[models.ScRNAPaper]:
    return db.query(models.ScRNAPaper).offset(skip).limit(limit).all()

def update_paper(db: Session, paper_id: str, paper: schemas.ScRNAPaperUpdate) -> Optional[models.ScRNAPaper]:
    db_paper = get_paper(db, paper_id)
    if db_paper:
        update_data = paper.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(db_paper, field, value)
        db.commit()
        db.refresh(db_paper)
    return db_paper

def delete_paper(db: Session, paper_id: str) -> bool:
    db_paper = get_paper(db, paper_id)
    if db_paper:
        db.delete(db_paper)
        db.commit()
        return True
    return False

# Paper search helpers
def search_papers_by_disease(db: Session, disease: str, skip: int = 0, limit: int = 100) -> List[models.ScRNAPaper]:
    return db.query(models.ScRNAPaper).filter(models.ScRNAPaper.disease.ilike(f"%{disease}%")).offset(skip).limit(limit).all()

def search_papers_by_species(db: Session, species: str, skip: int = 0, limit: int = 100) -> List[models.ScRNAPaper]:
    return db.query(models.ScRNAPaper).filter(models.ScRNAPaper.species.ilike(f"%{species}%")).offset(skip).limit(limit).all()

def search_papers_by_brain_region(db: Session, brain_region: str, skip: int = 0, limit: int = 100) -> List[models.ScRNAPaper]:
    return db.query(models.ScRNAPaper).filter(models.ScRNAPaper.brain_region.ilike(f"%{brain_region}%")).offset(skip).limit(limit).all()

# Spatial Dataset CRUD operations
def create_spatial_dataset(db: Session, dataset: schemas.SpatialDatasetCreate) -> models.SpatialDataset:
    db_dataset = models.SpatialDataset(**dataset.model_dump())
    db.add(db_dataset)
    db.commit()
    db.refresh(db_dataset)
    return db_dataset

def get_spatial_dataset(db: Session, dataset_id: str) -> Optional[models.SpatialDataset]:
    return db.query(models.SpatialDataset).filter(models.SpatialDataset.dataset_id == dataset_id).first()

def get_spatial_datasets(db: Session, skip: int = 0, limit: int = 100) -> List[models.SpatialDataset]:
    return db.query(models.SpatialDataset).offset(skip).limit(limit).all()

def update_spatial_dataset(db: Session, dataset_id: str, dataset: schemas.SpatialDatasetUpdate) -> Optional[models.SpatialDataset]:
    db_dataset = get_spatial_dataset(db, dataset_id)
    if db_dataset:
        update_data = dataset.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(db_dataset, field, value)
        db.commit()
        db.refresh(db_dataset)
    return db_dataset

def delete_spatial_dataset(db: Session, dataset_id: str) -> bool:
    db_dataset = get_spatial_dataset(db, dataset_id)
    if db_dataset:
        db.delete(db_dataset)
        db.commit()
        return True
    return False

# Search functions for spatial datasets
def search_spatial_datasets_by_disease(db: Session, disease: str, skip: int = 0, limit: int = 100) -> List[models.SpatialDataset]:
    return db.query(models.SpatialDataset).filter(
        models.SpatialDataset.disease.ilike(f"%{disease}%")
    ).offset(skip).limit(limit).all()

def search_spatial_datasets_by_species(db: Session, species: str, skip: int = 0, limit: int = 100) -> List[models.SpatialDataset]:
    return db.query(models.SpatialDataset).filter(
        models.SpatialDataset.species.ilike(f"%{species}%")
    ).offset(skip).limit(limit).all()

def search_spatial_datasets_by_brain_region(db: Session, brain_region: str, skip: int = 0, limit: int = 100) -> List[models.SpatialDataset]:
    return db.query(models.SpatialDataset).filter(
        models.SpatialDataset.brain_region.ilike(f"%{brain_region}%")
    ).offset(skip).limit(limit).all()

def search_spatial_datasets_by_methodology(db: Session, methodology: str, skip: int = 0, limit: int = 100) -> List[models.SpatialDataset]:
    return db.query(models.SpatialDataset).filter(
        models.SpatialDataset.methodology.ilike(f"%{methodology}%")
    ).offset(skip).limit(limit).all()

def search_spatial_datasets_by_data_id(db: Session, data_id: str, skip: int = 0, limit: int = 100) -> List[models.SpatialDataset]:
    """Get all datasets from a specific study/project (data_id)"""
    return db.query(models.SpatialDataset).filter(
        models.SpatialDataset.data_id == data_id
    ).offset(skip).limit(limit).all()

# Spatial Paper CRUD operations
def create_spatial_paper(db: Session, paper: schemas.SpatialPaperCreate) -> models.SpatialPaper:
    db_paper = models.SpatialPaper(**paper.model_dump())
    db.add(db_paper)
    db.commit()
    db.refresh(db_paper)
    return db_paper

def get_spatial_paper(db: Session, data_id: str) -> Optional[models.SpatialPaper]:
    return db.query(models.SpatialPaper).filter(models.SpatialPaper.data_id == data_id).first()

def get_spatial_papers(db: Session, skip: int = 0, limit: int = 100) -> List[models.SpatialPaper]:
    return db.query(models.SpatialPaper).offset(skip).limit(limit).all()

def update_spatial_paper(db: Session, data_id: str, paper: schemas.SpatialPaperUpdate) -> Optional[models.SpatialPaper]:
    db_paper = get_spatial_paper(db, data_id)
    if db_paper:
        update_data = paper.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(db_paper, field, value)
        db.commit()
        db.refresh(db_paper)
    return db_paper

def delete_spatial_paper(db: Session, data_id: str) -> bool:
    db_paper = get_spatial_paper(db, data_id)
    if db_paper:
        db.delete(db_paper)
        db.commit()
        return True
    return False

# Search functions for spatial papers
def search_spatial_papers_by_disease(db: Session, disease: str, skip: int = 0, limit: int = 100) -> List[models.SpatialPaper]:
    return db.query(models.SpatialPaper).filter(
        models.SpatialPaper.disease.ilike(f"%{disease}%")
    ).offset(skip).limit(limit).all()

def search_spatial_papers_by_species(db: Session, species: str, skip: int = 0, limit: int = 100) -> List[models.SpatialPaper]:
    return db.query(models.SpatialPaper).filter(
        models.SpatialPaper.species.ilike(f"%{species}%")
    ).offset(skip).limit(limit).all()

def search_spatial_papers_by_brain_region(db: Session, brain_region: str, skip: int = 0, limit: int = 100) -> List[models.SpatialPaper]:
    return db.query(models.SpatialPaper).filter(
        models.SpatialPaper.brain_region.ilike(f"%{brain_region}%")
    ).offset(skip).limit(limit).all()

def search_spatial_papers_by_methodology(db: Session, methodology: str, skip: int = 0, limit: int = 100) -> List[models.SpatialPaper]:
    return db.query(models.SpatialPaper).filter(
        models.SpatialPaper.methodology.ilike(f"%{methodology}%")
    ).offset(skip).limit(limit).all() 