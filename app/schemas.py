from pydantic import BaseModel, field_validator
from datetime import datetime
from typing import Optional, List, Dict, Any

# Base Message Schema
class MessageBase(BaseModel):
    content: str
    author: str

# Create Message Schema
class MessageCreate(MessageBase):
    pass

# Update Message Schema
class MessageUpdate(BaseModel):
    content: Optional[str] = None
    author: Optional[str] = None

# Message Response Schema
class Message(MessageBase):
    id: int
    created_at: datetime
    updated_at: Optional[datetime] = None
    
    class Config:
        from_attributes = True

# Hello World Schema
class HelloWorld(BaseModel):
    message: str
    timestamp: datetime

# ScRNA-seq Dataset Schemas (updated headers)
class ScRNAseqDatasetBase(BaseModel):
    dataset_id: str
    public_dataset_id: Optional[str] = None
    normation: Optional[str] = None
    pubmed_id: Optional[str] = None
    disease: Optional[str] = None
    status: Optional[str] = None
    control: Optional[str] = None
    species: Optional[str] = None
    brain_region: Optional[str] = None
    treatment: Optional[str] = None
    sex: Optional[str] = None
    stage: Optional[str] = None
    age: Optional[str] = None
    zhu_age: Optional[str] = None
    n_cells: Optional[int] = None
    protocol: Optional[str] = None
    methodology: Optional[str] = None
    model: Optional[str] = None

class ScRNAseqDatasetCreate(ScRNAseqDatasetBase):
    pass

class ScRNAseqDatasetUpdate(BaseModel):
    dataset_id: Optional[str] = None
    public_dataset_id: Optional[str] = None
    normation: Optional[str] = None
    pubmed_id: Optional[str] = None
    disease: Optional[str] = None
    status: Optional[str] = None
    control: Optional[str] = None
    species: Optional[str] = None
    brain_region: Optional[str] = None
    treatment: Optional[str] = None
    sex: Optional[str] = None
    stage: Optional[str] = None
    age: Optional[str] = None
    zhu_age: Optional[str] = None
    n_cells: Optional[int] = None
    protocol: Optional[str] = None
    methodology: Optional[str] = None
    model: Optional[str] = None

class ScRNAseqDataset(ScRNAseqDatasetBase):
    id: int
    created_at: datetime
    
    class Config:
        from_attributes = True

# ScRNA-seq Paper Schemas
class ScRNAPaperBase(BaseModel):
    paper_id: str
    public_data_id: Optional[str] = None
    disease: Optional[str] = None
    pubmed_id: Optional[str] = None
    protocol: Optional[str] = None
    methodology: Optional[str] = None
    title: Optional[str] = None
    author: Optional[str] = None
    citation: Optional[str] = None
    abstract: Optional[str] = None
    doi: Optional[str] = None
    date_published: Optional[str] = None
    brain_region: Optional[str] = None
    species: Optional[str] = None
    cell_numbers: Optional[str] = None

class ScRNAPaperCreate(ScRNAPaperBase):
    pass

class ScRNAPaperUpdate(BaseModel):
    paper_id: Optional[str] = None
    public_data_id: Optional[str] = None
    disease: Optional[str] = None
    pubmed_id: Optional[str] = None
    protocol: Optional[str] = None
    methodology: Optional[str] = None
    title: Optional[str] = None
    author: Optional[str] = None
    citation: Optional[str] = None
    abstract: Optional[str] = None
    doi: Optional[str] = None
    date_published: Optional[str] = None
    brain_region: Optional[str] = None
    species: Optional[str] = None
    cell_numbers: Optional[str] = None

class ScRNAPaper(ScRNAPaperBase):
    id: int
    created_at: datetime

    class Config:
        from_attributes = True

# Spatial Dataset Schemas
class SpatialDatasetBase(BaseModel):
    data_id: Optional[str] = None  # Study/project ID (can have duplicates)
    dataset_id: str  # Unique sample ID
    data_path: Optional[str] = None
    public_dataset_id: Optional[str] = None
    slide_id: Optional[str] = None
    disease: Optional[str] = None
    status: Optional[str] = None
    species: Optional[str] = None
    brain_region: Optional[str] = None
    treatment: Optional[str] = None
    sex: Optional[str] = None
    stage: Optional[str] = None
    age: Optional[str] = None
    n_cells: Optional[int] = None
    protocol: Optional[str] = None
    methodology: Optional[str] = None
    model: Optional[str] = None

class SpatialDatasetCreate(SpatialDatasetBase):
    pass

class SpatialDatasetUpdate(BaseModel):
    data_id: Optional[str] = None
    dataset_id: Optional[str] = None
    data_path: Optional[str] = None
    public_dataset_id: Optional[str] = None
    slide_id: Optional[str] = None
    disease: Optional[str] = None
    status: Optional[str] = None
    species: Optional[str] = None
    brain_region: Optional[str] = None
    treatment: Optional[str] = None
    sex: Optional[str] = None
    stage: Optional[str] = None
    age: Optional[str] = None
    n_cells: Optional[int] = None
    protocol: Optional[str] = None
    methodology: Optional[str] = None
    model: Optional[str] = None

class SpatialDataset(SpatialDatasetBase):
    id: int
    created_at: datetime
    
    class Config:
        from_attributes = True

# Spatial Paper Schemas
class SpatialPaperBase(BaseModel):
    data_id: Optional[str] = None
    public_data_id: Optional[str] = None
    disease: Optional[str] = None
    species: Optional[str] = None
    n_samples: Optional[int] = None
    pmid: Optional[str] = None
    library_prep_protocol: Optional[str] = None
    methodology: Optional[str] = None
    title: Optional[str] = None
    journal: Optional[str] = None
    author: Optional[str] = None
    citation: Optional[str] = None
    abstract: Optional[str] = None
    doi: Optional[str] = None
    publish_date: Optional[str] = None
    brain_region: Optional[str] = None
    drug: Optional[str] = None
    single_cell: Optional[str] = None
    scrna_index: Optional[str] = None
    n_sample_of_scrna: Optional[int] = None

class SpatialPaperCreate(SpatialPaperBase):
    pass

class SpatialPaperUpdate(BaseModel):
    data_id: Optional[str] = None
    public_data_id: Optional[str] = None
    disease: Optional[str] = None
    species: Optional[str] = None
    n_samples: Optional[int] = None
    pmid: Optional[str] = None
    library_prep_protocol: Optional[str] = None
    methodology: Optional[str] = None
    title: Optional[str] = None
    journal: Optional[str] = None
    author: Optional[str] = None
    citation: Optional[str] = None
    abstract: Optional[str] = None
    doi: Optional[str] = None
    publish_date: Optional[str] = None
    brain_region: Optional[str] = None
    drug: Optional[str] = None
    single_cell: Optional[str] = None
    scrna_index: Optional[str] = None
    n_sample_of_scrna: Optional[int] = None

class SpatialPaper(SpatialPaperBase):
    id: int
    created_at: datetime

    class Config:
        from_attributes = True

# H5AD Schemas
class H5ADFileBase(BaseModel):
    dataset_id: str
    file_path: str
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    has_umap: bool = False
    has_tsne: bool = False
    has_pca: bool = False
    available_embeddings: Optional[List[str]] = None
    available_metadata_columns: Optional[List[str]] = None
    precomputed_path: Optional[str] = None
    is_preprocessed: bool = False

class H5ADFileCreate(H5ADFileBase):
    pass

class H5ADFile(H5ADFileBase):
    id: int
    created_at: datetime
    updated_at: Optional[datetime] = None

    class Config:
        from_attributes = True

# API Response Schemas
class EmbeddingResponse(BaseModel):
    embedding_type: str
    coordinates: List[List[float]]
    cell_ids: List[str]

class GeneExpressionResponse(BaseModel):
    gene: str
    expression: List[float]
    cell_ids: List[str]
    stats: dict

class MetadataResponse(BaseModel):
    columns: List[str]
    data: dict
    cell_ids: List[str]

class PlotDataResponse(BaseModel):
    coordinates: List[List[float]]
    genes: Optional[dict] = None
    metadata: Optional[dict] = None
    cell_ids: List[str]

class GeneSearchResult(BaseModel):
    symbol: str
    mean_expression: float
    expressed_cells: int

class ModuleScoreResponse(BaseModel):
    module_score: List[float]
    cell_ids: List[str]
    genes_used: List[str]
    genes_not_found: List[str]
    method: str
    n_genes_used: int
    stats: dict
    data_type: str

class ModuleScoreRequest(BaseModel):
    gene_list: List[str]
    use_raw: bool = False

# Heatmap/Dotplot Schemas
class DotplotDataPoint(BaseModel):
    expression: float
    percentage: float

class HeatmapRequest(BaseModel):
    gene_list: List[str]
    plot_type: str = "heatmap"  # "heatmap" or "dotplot"
    scale_expression: bool = False  # Z-score normalization per gene
    cluster_rows: bool = False  # Cluster genes (rows)
    cluster_columns: bool = False  # Cluster cell types (columns)
    group_by: Optional[str] = "cell_type"  # Metadata column to group by (e.g., annotation, seurat_clusters)

class HeatmapResponse(BaseModel):
    data: List[List[Any]]  # Can be List[List[float]] for heatmap or List[List[DotplotDataPoint]] for dotplot
    genes: List[str]
    cell_types: List[str]
    genes_used: List[str]
    genes_not_found: List[str]
    plot_type: str
    scale_expression: bool
    cluster_rows: bool
    cluster_columns: bool
    stats: dict

# DEG Analysis Schemas
class DEGGeneResult(BaseModel):
    gene: str
    log2_fold_change: float
    mean_expr_group1: float
    mean_expr_group2: float
    pct_group1: float
    pct_group2: float
    p_value: float
    p_value_adj: float
    significant: bool

class DEGGroupInfo(BaseModel):
    dataset_id: Optional[str] = None
    n_cells: int
    cell_type: Optional[str] = None
    cell_types: Optional[List[str]] = None
    filters: Optional[Dict[str, str]] = None

class DEGComparisonInfo(BaseModel):
    group1: DEGGroupInfo
    group2: DEGGroupInfo
    dataset_id: Optional[str] = None  # For within-dataset comparison

class DEGParameters(BaseModel):
    min_pct: float
    logfc_threshold: float
    p_value_threshold: float

class DEGSummary(BaseModel):
    total_genes_tested: int
    significant_genes: int
    upregulated_in_group1: int
    upregulated_in_group2: int

class DEGAnalysisResult(BaseModel):
    comparison: dict
    parameters: DEGParameters
    summary: DEGSummary
    genes: List[DEGGeneResult]

# Request schemas for DEG
class DEGBetweenDatasetsRequest(BaseModel):
    dataset_id1: str
    dataset_id2: str
    cell_type: Optional[str] = None
    cell_type2: Optional[str] = None
    cell_types_group1: Optional[List[str]] = None
    cell_types_group2: Optional[List[str]] = None
    metadata_filters1: Optional[Dict[str, str]] = None
    metadata_filters2: Optional[Dict[str, str]] = None
    min_pct: float = 0.1
    logfc_threshold: float = 0.25
    p_value_threshold: float = 0.05
    top_n: Optional[int] = None

class DEGWithinDatasetRequest(BaseModel):
    dataset_id: str
    group1_filters: Dict[str, str]
    group2_filters: Dict[str, str]
    min_pct: float = 0.1
    logfc_threshold: float = 0.25
    p_value_threshold: float = 0.05
    top_n: Optional[int] = None

# Advanced Atlas DEG Analysis Schemas
class MetadataFilter(BaseModel):
    column: str
    values: List[str]

class AtlasDEGRequest(BaseModel):
    dataset_id: str
    group1_filters: List[MetadataFilter]
    group2_filters: List[MetadataFilter]
    logfc_threshold: float = 0.25
    p_value_threshold: float = 0.05
    min_pct: float = 0.1
    top_n: Optional[int] = None

class GroupDescription(BaseModel):
    filters_applied: List[MetadataFilter]
    cell_count: int

class AtlasDEGSummary(BaseModel):
    group1_cells: int
    group2_cells: int
    total_genes_tested: int
    significant_genes: int
    comparison_description: str

class AtlasDEGResult(BaseModel):
    summary: AtlasDEGSummary
    group1_description: GroupDescription
    group2_description: GroupDescription
    genes: List[DEGGeneResult]

class PreviewCellCountsRequest(BaseModel):
    dataset_id: str
    group1_filters: List[MetadataFilter]
    group2_filters: List[MetadataFilter]

class PreviewCellCountsResponse(BaseModel):
    group1_count: int
    group2_count: int
    overlap_count: int
    warnings: List[str] = []

class MetadataValuesResponse(BaseModel):
    column: str
    values: List[str]
    counts: Dict[str, int]

# Spatial Transcriptomics Schemas
class SpatialCoordinatesResponse(BaseModel):
    coordinates: List[List[float]]
    cell_ids: List[str]
    spatial_key: Optional[str] = None

class SpatialPlotDataResponse(BaseModel):
    coordinates: List[List[float]]
    cell_ids: List[str]
    genes: Optional[Dict[str, List[float]]] = None
    metadata: Optional[Dict[str, List[Any]]] = None

class SVGGeneResult(BaseModel):
    gene: str
    gft_score: Optional[float] = None
    svg_rank: Optional[float] = None
    cutoff_gft_score: Optional[float] = None
    pvalue: Optional[float] = None
    fdr: Optional[float] = None
    
    @field_validator('cutoff_gft_score', mode='before')
    @classmethod
    def validate_cutoff_gft_score(cls, v):
        """Handle boolean/string values in cutoff_gft_score"""
        if v is None:
            return None
        if isinstance(v, bool):
            return float(v)
        if isinstance(v, str):
            if v.lower() in ('true', '1'):
                return 1.0
            if v.lower() in ('false', '0'):
                return 0.0
            try:
                return float(v)
            except (ValueError, TypeError):
                return None
        try:
            return float(v)
        except (ValueError, TypeError):
            return None

class SVGResponse(BaseModel):
    genes: List[SVGGeneResult]
    total_genes: int
    filtered_count: int
    available_columns: List[str]

class PrecomputedDEGGroupResult(BaseModel):
    gene: str
    score: Optional[float] = None
    logfoldchange: Optional[float] = None
    pval: Optional[float] = None
    pval_adj: Optional[float] = None

class PrecomputedDEGResponse(BaseModel):
    groups: List[str]
    params: Dict[str, Any]
    genes: Dict[str, List[PrecomputedDEGGroupResult]]
    method: str

class DeconvolutionResponse(BaseModel):
    predictions: Dict[str, List[float]]
    cell_ids: List[str]
    cell_types: List[str]

class CCCInteractionsResponse(BaseModel):
    interactions: List[Any]
    ligands: List[str]
    receptors: List[str]
    cell_types: List[str]
    data: Dict[str, Any]

class SpatialInfoResponse(BaseModel):
    has_spatial_coordinates: bool
    has_svg: bool
    has_precomputed_deg: bool
    has_deconvolution: bool
    has_ccc: bool
    dataset_type: Optional[str] = None
    has_spatial_image: Optional[bool] = False

class SpatialImageSampleInfo(BaseModel):
    image_keys: List[str]
    image_shapes: Dict[str, List[int]]
    scalefactors: Dict[str, float]
    image_paths: Dict[str, str]

class SpatialImageInfoResponse(BaseModel):
    samples: Dict[str, SpatialImageSampleInfo]
    default_sample: str

class SpatialCoordinatesTransformedResponse(BaseModel):
    coordinates: List[List[float]]
    cell_ids: List[str]
    sample_key: str
    image_key: str
    scale_factor: float
    original_coordinates: List[List[float]]
    image_shape: Optional[List[int]] = None

# H5AD Analysis Features Schemas
class H5ADAnalysisFeaturesBase(BaseModel):
    dataset_id: str
    has_umap: bool = False
    has_tsne: bool = False
    has_pca: bool = False
    has_spatial_coordinates: bool = False
    has_svg: bool = False
    has_precomputed_deg: bool = False
    has_deconvolution: bool = False
    has_ccc: bool = False
    dataset_type: Optional[str] = None
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    n_deg_groups: Optional[int] = None
    n_deconv_cell_types: Optional[int] = None

class H5ADAnalysisFeaturesCreate(H5ADAnalysisFeaturesBase):
    pass

class H5ADAnalysisFeaturesUpdate(BaseModel):
    has_umap: Optional[bool] = None
    has_tsne: Optional[bool] = None
    has_pca: Optional[bool] = None
    has_spatial_coordinates: Optional[bool] = None
    has_svg: Optional[bool] = None
    has_precomputed_deg: Optional[bool] = None
    has_deconvolution: Optional[bool] = None
    has_ccc: Optional[bool] = None
    dataset_type: Optional[str] = None
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    n_deg_groups: Optional[int] = None
    n_deconv_cell_types: Optional[int] = None

class H5ADAnalysisFeatures(H5ADAnalysisFeaturesBase):
    id: int
    created_at: datetime
    updated_at: Optional[datetime] = None
    
    class Config:
        from_attributes = True 