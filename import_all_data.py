#!/usr/bin/env python3
"""
Combined script to import all data types (scRNA-seq, spatial data, spatial papers) into SQLite database
Run this script from the project root directory
"""
import argparse
import logging
import subprocess
import sys
import os

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_import_script(script_name: str, csv_path: str, replace: bool = False) -> bool:
    """Run an import script and return success status"""
    # Scripts are in scripts/ directory
    script_path = os.path.join("scripts", script_name)
    cmd = [sys.executable, script_path, "--csv", csv_path]
    if replace:
        cmd.append("--replace")
    
    logger.info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"✓ {script_name} completed successfully")
        if result.stdout:
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"✗ {script_name} failed with return code {e.returncode}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="Import all data types into SQLite database")
    parser.add_argument("--replace", action="store_true", help="Drop and recreate all tables before import")
    parser.add_argument("--scrna-data", default="data/ssKIND - scRNAseq_data_for_ssKIND.csv", help="Path to scRNA-seq data CSV")
    parser.add_argument("--scrna-papers", default="data/ssKIND - scRNAseq_paper_for_ssKIND.csv", help="Path to scRNA-seq papers CSV")
    parser.add_argument("--spatial-data", default="data/ssKIND - Spatial data for ssKIND.csv", help="Path to spatial data CSV")
    parser.add_argument("--spatial-papers", default="data/ssKIND - Spatial paper for ssKIND.csv", help="Path to spatial papers CSV")
    args = parser.parse_args()

    logger.info("Starting combined data import")
    if args.replace:
        logger.warning("Replace mode enabled: all tables will be rebuilt before import")

    success_count = 0
    total_imports = 4

    # Import scRNA-seq data
    if run_import_script("import_csv.py", args.scrna_data, args.replace):
        success_count += 1

    # Import scRNA-seq papers
    if run_import_script("import_papers.py", args.scrna_papers, args.replace):
        success_count += 1

    # Import spatial data
    if run_import_script("import_spatial_data.py", args.spatial_data, args.replace):
        success_count += 1

    # Import spatial papers
    if run_import_script("import_spatial_papers.py", args.spatial_papers, args.replace):
        success_count += 1

    logger.info(f"Import completed: {success_count}/{total_imports} scripts succeeded")
    
    if success_count == total_imports:
        logger.info("🎉 All imports completed successfully!")
        return 0
    else:
        logger.error(f"❌ {total_imports - success_count} imports failed")
        return 1

if __name__ == "__main__":
    exit(main())
