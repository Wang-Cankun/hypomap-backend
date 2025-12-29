import os
from typing import List
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

class Settings:
    # Application settings
    PORT: int = int(os.getenv("PORT", 9120))
    HOST: str = os.getenv("HOST", "0.0.0.0")
    DEBUG: bool = os.getenv("DEBUG", "True").lower() == "true"
    RELOAD: bool = os.getenv("RELOAD", "True").lower() == "true"

    # Database settings
    DATABASE_URL: str = os.getenv("DATABASE_URL", "sqlite:///./app.db")

    # API settings
    GLOBAL_PREFIX: str = os.getenv("GLOBAL_PREFIX", "/hypomap-backend")
    API_PREFIX: str = os.getenv("API_PREFIX", "/api/v1")
    DOCS_URL: str = os.getenv("DOCS_URL", "/docs")
    REDOC_URL: str = os.getenv("REDOC_URL", "/redoc")

    # CORS settings
    CORS_ORIGINS: List[str] = os.getenv("CORS_ORIGINS", '["*"]').replace('"', '').strip('[]').split(',')
    CORS_CREDENTIALS: bool = os.getenv("CORS_CREDENTIALS", "True").lower() == "true"
    CORS_METHODS: List[str] = os.getenv("CORS_METHODS", '["*"]').replace('"', '').strip('[]').split(',')
    CORS_HEADERS: List[str] = os.getenv("CORS_HEADERS", '["*"]').replace('"', '').strip('[]').split(',')

    # Application metadata
    APP_TITLE: str = "HypoMap Server"
    APP_DESCRIPTION: str = "Single-cell RNA-seq analysis platform with AI-powered hypothesis generation"
    APP_VERSION: str = "1.0.0"

    # Data directories
    H5AD_RAW_DIR: str = os.getenv("H5AD_RAW_DIR", "h5ad/raw")
    H5AD_PRECOMPUTED_DIR: str = os.getenv("H5AD_PRECOMPUTED_DIR", "h5ad/precomputed")

    # AI settings
    ANTHROPIC_API_KEY: str = os.getenv("ANTHROPIC_API_KEY", "")

# Create settings instance
settings = Settings() 