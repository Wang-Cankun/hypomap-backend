from app.config import settings

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app.app:app",  # Use import string format for proper reload
        host=settings.HOST, 
        port=settings.PORT,
        reload=settings.RELOAD
    )