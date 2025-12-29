#!/usr/bin/env python3
"""
Simple script to import spatial papers from CSV into SQLite database
Run this script from the project root directory
"""
import csv
import os
import sys
import argparse
from typing import Dict, Any, Optional
import logging

# Add parent directory to path to import app modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy.exc import IntegrityError

from app.database import SessionLocal, engine
from app.models import Base, SpatialPaper
from app.schemas import SpatialPaperCreate

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def clean_value(value: str) -> Optional[str]:
    """Clean and validate CSV values"""
    if value is None or value.strip() == '':
        return None
    return value.strip()

def parse_integer(value: str) -> Optional[int]:
    """Parse integer fields"""
    if not value or value.strip() == '':
        return None
    try:
        # Remove any non-numeric characters and convert to int
        cleaned = ''.join(filter(str.isdigit, value.strip()))
        return int(cleaned) if cleaned else None
    except (ValueError, TypeError):
        return None

def create_paper_from_row(row: Dict[str, Any]) -> SpatialPaperCreate:
    """Create a SpatialPaperCreate object from a CSV row"""
    return SpatialPaperCreate(
        data_id=clean_value(row.get('data_id')),
        public_data_id=clean_value(row.get('public_data_id')),
        disease=clean_value(row.get('Disease')),
        species=clean_value(row.get('species')),
        n_samples=parse_integer(row.get('n_samples')),
        pmid=clean_value(row.get('pmid')),
        library_prep_protocol=clean_value(row.get('library_prep_protocol')),
        methodology=clean_value(row.get('methodology')),
        title=clean_value(row.get('title')),
        journal=clean_value(row.get('journal')),
        author=clean_value(row.get('author')),
        citation=clean_value(row.get('citation')),
        abstract=clean_value(row.get('abstract')),
        doi=clean_value(row.get('doi')),
        publish_date=clean_value(row.get('Publish date')),
        brain_region=clean_value(row.get('Brain Region')),
        drug=clean_value(row.get('drug?')),
        single_cell=clean_value(row.get('single cell?')),
        scrna_index=clean_value(row.get('scRNA index')),
        n_sample_of_scrna=parse_integer(row.get('n_sample of scRNA')),
    )

def prepare_table(replace: bool) -> None:
    """Ensure table exists; if replace, drop and recreate the paper table."""
    if replace:
        logger.warning("--replace specified: dropping and recreating 'spatial_papers' table")
        # Drop only the target table if it exists
        SpatialPaper.__table__.drop(bind=engine, checkfirst=True)
    # Create tables as needed
    Base.metadata.create_all(bind=engine)

def commit_batch_with_fallback(db, batch, duplicate_rows_counter_ref):
    """Try to commit a batch, and on IntegrityError fall back to per-row commits, skipping duplicates."""
    try:
        db.add_all(batch)
        db.commit()
        return len(batch)
    except IntegrityError as e:
        db.rollback()
        inserted = 0
        for item in batch:
            try:
                db.add(item)
                db.commit()
                inserted += 1
            except IntegrityError:
                db.rollback()
                duplicate_rows_counter_ref[0] += 1
        return inserted

def import_csv_to_database(csv_path: str, batch_size: int = 100, replace: bool = False) -> None:
    """Import CSV data into SQLite database"""
    
    # Check if CSV file exists
    if not os.path.exists(csv_path):
        logger.error(f"CSV file not found: {csv_path}")
        return
    
    # Prepare table (optionally drop + recreate)
    prepare_table(replace)
    
    db = SessionLocal()
    try:
        # If not replacing, warn if data already exists
        if not replace:
            existing_count = db.query(SpatialPaper).count()
            if existing_count > 0:
                logger.warning(f"Database already contains {existing_count} spatial papers")
                response = input("Do you want to continue and add more data? (y/N): ")
                if response.lower() != 'y':
                    logger.info("Import cancelled by user")
                    return
        
        # Read and import CSV
        logger.info(f"Starting import from {csv_path}")
        
        with open(csv_path, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            total_rows = 0
            imported_rows = 0
            skipped_rows = 0
            duplicate_rows = 0
            batch = []
            seen_ids = set()
            
            for row in reader:
                total_rows += 1
                
                try:
                    # Create paper object
                    paper_data = create_paper_from_row(row)
                    
                    # Skip if data_id is missing
                    if not paper_data.data_id:
                        logger.warning(f"Row {total_rows}: Missing data_id, skipping")
                        skipped_rows += 1
                        continue
                    
                    # De-duplicate within the same CSV
                    if paper_data.data_id in seen_ids:
                        duplicate_rows += 1
                        continue
                    seen_ids.add(paper_data.data_id)
                    
                    # If not replacing, check if paper already exists in DB
                    if not replace:
                        existing = db.query(SpatialPaper).filter(
                            SpatialPaper.data_id == paper_data.data_id
                        ).first()
                        if existing:
                            duplicate_rows += 1
                            continue
                    
                    # Create new paper using model_dump() for Pydantic v2
                    db_paper = SpatialPaper(**paper_data.model_dump())
                    batch.append(db_paper)
                    
                    # Commit batch
                    if len(batch) >= batch_size:
                        inserted = commit_batch_with_fallback(db, batch, [duplicate_rows])
                        imported_rows += inserted
                        logger.info(f"Imported batch: {imported_rows} spatial papers so far")
                        batch = []
                
                except Exception as e:
                    logger.error(f"Error processing row {total_rows}: {e}")
                    skipped_rows += 1
                    continue
            
            # Commit remaining batch
            if batch:
                inserted = commit_batch_with_fallback(db, batch, [duplicate_rows])
                imported_rows += inserted
        
        logger.info(f"Import completed!")
        logger.info(f"Total rows processed: {total_rows}")
        logger.info(f"Successfully imported: {imported_rows}")
        logger.info(f"Skipped rows (missing/invalid): {skipped_rows}")
        logger.info(f"Duplicate rows skipped: {duplicate_rows}")
        
    except Exception as e:
        logger.error(f"Error during import: {e}")
        db.rollback()
        raise
    finally:
        db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Import spatial papers CSV into SQLite")
    parser.add_argument("--csv", dest="csv_path", default="../data/ssKIND - Spatial paper for ssKIND.csv", help="Path to CSV file")
    parser.add_argument("--replace", action="store_true", help="Drop and recreate paper table before import")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for inserts")
    args = parser.parse_args()

    logger.info("Starting spatial papers import")
    logger.info(f"CSV file: {args.csv_path}")
    if args.replace:
        logger.warning("Replace mode enabled: table will be rebuilt before import")
    
    try:
        import_csv_to_database(args.csv_path, batch_size=args.batch_size, replace=args.replace)
        logger.info("Import completed successfully!")
    except Exception as e:
        logger.error(f"Import failed: {e}")
        exit(1)
