#!/usr/bin/env python3
"""
Simple script to run the FastAPI application
"""
import uvicorn
from app.config import settings

if __name__ == "__main__":
    print("🚀 Starting Modern FastAPI with SQLite...")
    print(f"🌐 Server: http://{settings.HOST}:{settings.PORT}")
    print(f"📚 API Documentation: http://{settings.HOST}:{settings.PORT}{settings.GLOBAL_PREFIX}{settings.DOCS_URL}")
    print(f"🔍 Alternative Docs: http://{settings.HOST}:{settings.PORT}{settings.GLOBAL_PREFIX}{settings.REDOC_URL}")
    print(f"💚 Health Check: http://{settings.HOST}:{settings.PORT}{settings.GLOBAL_PREFIX}/health")
    print(f"🔗 API Base URL: http://{settings.HOST}:{settings.PORT}{settings.GLOBAL_PREFIX}{settings.API_PREFIX}")
    print("=" * 50)
    
    uvicorn.run(
        "app.app:app",  # Use import string format for proper reload
        host=settings.HOST,
        port=settings.PORT,
        reload=settings.RELOAD,
        log_level="info"
    ) 