#!/usr/bin/env python3
"""
Script to import scRNA-seq papers from CSV into SQLite database
"""
import csv
import sys
import os
import argparse
from pathlib import Path
from typing import Dict, Any, Optional
import logging

from sqlalchemy.exc import IntegrityError

# Add the parent directory to the path so we can import from app
sys.path.append(str(Path(__file__).parent.parent))

from app.database import SessionLocal, engine
from app.models import Base, ScRNAPaper
from app.schemas import ScRNAPaperCreate

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def clean(value: str) -> Optional[str]:
    if value is None or value.strip() == '':
        return None
    return value.strip()

def row_to_paper(row: Dict[str, Any]) -> ScRNAPaperCreate:
    return ScRNAPaperCreate(
        paper_id=clean(row.get('Paper_ID')),
        public_data_id=clean(row.get('Public_data_id')),
        disease=clean(row.get('Disease')),
        pubmed_id=clean(row.get('Pubmed_id')),
        protocol=clean(row.get('protocol')),
        methodology=clean(row.get('methodology')),
        title=clean(row.get('title')),
        author=clean(row.get('author')),
        citation=clean(row.get('citation')),
        abstract=clean(row.get('abstract')),
        doi=clean(row.get('doi')),
        date_published=clean(row.get('Date_Published')),
        brain_region=clean(row.get('Brain_Region')),
        species=clean(row.get('Species')),
        cell_numbers=clean(row.get('Cell numbers')),
    )

def prepare_table(replace: bool) -> None:
    if replace:
        logger.warning("--replace specified: dropping and recreating 'scrna_papers' table")
        ScRNAPaper.__table__.drop(bind=engine, checkfirst=True)
    Base.metadata.create_all(bind=engine)

def commit_batch_with_fallback(db, batch, duplicate_counter_ref):
    try:
        db.add_all(batch)
        db.commit()
        return len(batch)
    except IntegrityError:
        db.rollback()
        inserted = 0
        for item in batch:
            try:
                db.add(item)
                db.commit()
                inserted += 1
            except IntegrityError:
                db.rollback()
                duplicate_counter_ref[0] += 1
        return inserted

def import_csv(csv_path: str, batch_size: int = 100, replace: bool = False):
    if not os.path.exists(csv_path):
        logger.error(f"CSV file not found: {csv_path}")
        return

    prepare_table(replace)

    db = SessionLocal()
    try:
        if not replace:
            existing = db.query(ScRNAPaper).count()
            if existing > 0:
                logger.warning(f"Database already contains {existing} papers")
                resp = input("Continue and add more? (y/N): ")
                if resp.lower() != 'y':
                    logger.info("Import cancelled")
                    return

        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            total = 0
            imported = 0
            skipped = 0
            dup = 0
            seen = set()
            batch = []

            for row in reader:
                total += 1
                try:
                    paper = row_to_paper(row)
                    if not paper.paper_id:
                        skipped += 1
                        continue
                    if paper.paper_id in seen:
                        dup += 1
                        continue
                    seen.add(paper.paper_id)

                    if not replace:
                        if db.query(ScRNAPaper).filter(ScRNAPaper.paper_id == paper.paper_id).first():
                            dup += 1
                            continue

                    db_item = ScRNAPaper(**paper.model_dump())
                    batch.append(db_item)

                    if len(batch) >= batch_size:
                        inserted = commit_batch_with_fallback(db, batch, [dup])
                        imported += inserted
                        logger.info(f"Imported batch: {imported} papers so far")
                        batch = []
                except Exception as e:
                    skipped += 1
                    logger.error(f"Error row {total}: {e}")
                    continue

            if batch:
                inserted = commit_batch_with_fallback(db, batch, [dup])
                imported += inserted

        logger.info("Import completed!")
        logger.info(f"Total rows processed: {total}")
        logger.info(f"Successfully imported: {imported}")
        logger.info(f"Skipped rows (missing/invalid): {skipped}")
        logger.info(f"Duplicate rows skipped: {dup}")

    except Exception as e:
        db.rollback()
        logger.error(f"Error during import: {e}")
        raise
    finally:
        db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Import papers CSV into SQLite")
    parser.add_argument("--csv", dest="csv_path", default=str(Path(__file__).parent.parent / "data/ssKIND - scRNAseq_paper_for_ssKIND.csv"), help="Path to CSV file")
    parser.add_argument("--replace", action="store_true", help="Drop and recreate papers table before import")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for inserts")
    args = parser.parse_args()

    logger.info(f"CSV file: {args.csv_path}")
    if args.replace:
        logger.warning("Replace mode enabled: table will be rebuilt before import")

    import_csv(args.csv_path, batch_size=args.batch_size, replace=args.replace)
