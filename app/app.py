from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.database import engine
from app import models
from app.api.endpoints import router
from app.api.h5ad_endpoints import router as h5ad_router
from app.api.deg_endpoints import router as deg_router
from app.api.api_metadata import router as api_metadata_router
from app.api.upload_endpoints import router as upload_router
from app.api.ai_endpoints import router as ai_router
from app.config import settings

# Create database tables
models.Base.metadata.create_all(bind=engine)

# Create FastAPI app
app = FastAPI(
    title=settings.APP_TITLE,
    description=settings.APP_DESCRIPTION,
    version=settings.APP_VERSION,
    docs_url=f"{settings.GLOBAL_PREFIX}{settings.DOCS_URL}",
    redoc_url=f"{settings.GLOBAL_PREFIX}{settings.REDOC_URL}"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=settings.CORS_CREDENTIALS,
    allow_methods=settings.CORS_METHODS,
    allow_headers=settings.CORS_HEADERS,
)

# Include API routes with global prefix
app.include_router(router, prefix=f"{settings.GLOBAL_PREFIX}{settings.API_PREFIX}")
app.include_router(h5ad_router, prefix=f"{settings.GLOBAL_PREFIX}{settings.API_PREFIX}/h5ad", tags=["H5AD Data"])
app.include_router(deg_router, prefix=f"{settings.GLOBAL_PREFIX}{settings.API_PREFIX}/deg", tags=["DEG Analysis"])
app.include_router(api_metadata_router, prefix=f"{settings.GLOBAL_PREFIX}{settings.API_PREFIX}", tags=["API Metadata"])
app.include_router(upload_router, prefix=f"{settings.GLOBAL_PREFIX}{settings.API_PREFIX}", tags=["Upload"])
app.include_router(ai_router, prefix=f"{settings.GLOBAL_PREFIX}{settings.API_PREFIX}", tags=["AI Analysis"])

@app.get(f"{settings.GLOBAL_PREFIX}/health")
def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "message": "API is running"} 