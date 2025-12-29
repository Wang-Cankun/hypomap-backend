from sqlalchemy import Column, Integer, String, DateTime, Text, Boolean, JSON
from sqlalchemy.sql import func
from app.database import Base

class Message(Base):
    __tablename__ = "messages"
    
    id = Column(Integer, primary_key=True, index=True)
    content = Column(String(500), nullable=False)
    author = Column(String(100), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())

class ScRNAseqDataset(Base):
    __tablename__ = "scrna_datasets"
    
    id = Column(Integer, primary_key=True, index=True)
    dataset_id = Column(String(50), unique=True, index=True, nullable=False)
    public_dataset_id = Column(Text)
    normation = Column(String(100))
    pubmed_id = Column(Text)
    disease = Column(String(100))
    status = Column(String(100))
    control = Column(String(100))
    species = Column(String(100))
    brain_region = Column(String(200))
    treatment = Column(String(200))
    sex = Column(String(50))
    stage = Column(String(100))
    age = Column(String(100))
    zhu_age = Column(String(100))
    n_cells = Column(Integer)
    protocol = Column(Text)
    methodology = Column(Text)
    model = Column(String(200))
    created_at = Column(DateTime(timezone=True), server_default=func.now())

class ScRNAPaper(Base):
    __tablename__ = "scrna_papers"

    id = Column(Integer, primary_key=True, index=True)
    paper_id = Column(String(50), index=True, nullable=False)
    public_data_id = Column(Text)
    disease = Column(String(100))
    pubmed_id = Column(Text)
    protocol = Column(String(200))
    methodology = Column(String(200))
    title = Column(Text)
    author = Column(Text)
    citation = Column(Text)
    abstract = Column(Text)
    doi = Column(String(200))
    date_published = Column(String(50))
    brain_region = Column(String(200))
    species = Column(String(100))
    cell_numbers = Column(String(100))
    created_at = Column(DateTime(timezone=True), server_default=func.now())

class SpatialDataset(Base):
    __tablename__ = "spatial_datasets"
    
    id = Column(Integer, primary_key=True, index=True)
    data_id = Column(String(50), index=True)  # Study/project ID (can have duplicates)
    dataset_id = Column(String(50), unique=True, index=True, nullable=False)  # Unique sample ID
    data_path = Column(Text)
    public_dataset_id = Column(Text)
    slide_id = Column(String(100))
    disease = Column(String(100))
    status = Column(String(100))
    species = Column(String(100))
    brain_region = Column(String(200))
    treatment = Column(String(200))
    sex = Column(String(50))
    stage = Column(String(100))
    age = Column(String(100))
    n_cells = Column(Integer)
    protocol = Column(Text)
    methodology = Column(Text)
    model = Column(String(200))
    created_at = Column(DateTime(timezone=True), server_default=func.now())

class SpatialPaper(Base):
    __tablename__ = "spatial_papers"

    id = Column(Integer, primary_key=True, index=True)
    data_id = Column(String(50), index=True)
    public_data_id = Column(Text)
    disease = Column(String(100))
    species = Column(String(100))
    n_samples = Column(Integer)
    pmid = Column(Text)
    library_prep_protocol = Column(String(200))
    methodology = Column(String(200))
    title = Column(Text)
    journal = Column(String(200))
    author = Column(Text)
    citation = Column(Text)
    abstract = Column(Text)
    doi = Column(String(200))
    publish_date = Column(String(50))
    brain_region = Column(String(200))
    drug = Column(String(200))
    single_cell = Column(String(50))
    scrna_index = Column(String(100))
    n_sample_of_scrna = Column(Integer)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

class H5ADAnalysisFeatures(Base):
    __tablename__ = "h5ad_analysis_features"
    
    id = Column(Integer, primary_key=True, index=True)
    dataset_id = Column(String(100), unique=True, index=True, nullable=False)
    
    # Basic features
    has_umap = Column(Boolean, default=False)
    has_tsne = Column(Boolean, default=False)
    has_pca = Column(Boolean, default=False)
    
    # Spatial features
    has_spatial_coordinates = Column(Boolean, default=False)
    has_svg = Column(Boolean, default=False)
    has_precomputed_deg = Column(Boolean, default=False)
    has_deconvolution = Column(Boolean, default=False)
    has_ccc = Column(Boolean, default=False)
    
    # Dataset type
    dataset_type = Column(String(50), nullable=True)  # 'visium', 'xenium', 'unknown', None
    
    # Additional metadata
    n_cells = Column(Integer, nullable=True)
    n_genes = Column(Integer, nullable=True)
    n_deg_groups = Column(Integer, nullable=True)  # Number of DEG comparison groups
    n_deconv_cell_types = Column(Integer, nullable=True)  # Number of deconvolution cell types
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())

class H5ADFile(Base):
    __tablename__ = "h5ad_files"
    
    id = Column(Integer, primary_key=True, index=True)
    dataset_id = Column(String(50), unique=True, index=True, nullable=False)
    file_path = Column(String(500), nullable=False)
    n_cells = Column(Integer)
    n_genes = Column(Integer)
    has_umap = Column(Boolean, default=False)
    has_tsne = Column(Boolean, default=False)
    has_pca = Column(Boolean, default=False)
    available_embeddings = Column(JSON)  # ["umap", "tsne", "pca"]
    available_metadata_columns = Column(JSON)  # ["cell_type", "cluster", ...]
    precomputed_path = Column(String(500))  # Path to precomputed data directory
    is_preprocessed = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now()) 