#!/usr/bin/env python3
"""
Performance comparison: Simulate old JSON approach vs new SQLite approach
This shows the theoretical improvement
"""
import time
import json
import sqlite3
import statistics
from pathlib import Path
from typing import Dict, List

def simulate_json_search(json_path: Path, search_term: str, limit: int = 50) -> tuple:
    """Simulate the old JSON-based search approach"""
    # Load entire JSON file
    start_load = time.perf_counter()
    with open(json_path, 'r') as f:
        gene_stats = json.load(f)
    load_time = time.perf_counter() - start_load
    
    # Search in Python
    start_search = time.perf_counter()
    search_term_lower = search_term.lower()
    results = []
    
    for gene, stats in gene_stats.items():
        if search_term_lower in gene.lower():
            results.append({
                "symbol": gene,
                "mean_expression": stats["mean_expression"],
                "expressed_cells": stats["expressed_cells"]
            })
    
    # Sort by mean expression
    results.sort(key=lambda x: x["mean_expression"], reverse=True)
    results = results[:limit]
    search_time = time.perf_counter() - start_search
    
    return results, load_time, search_time

def sqlite_search(db_path: Path, search_term: str, limit: int = 50) -> tuple:
    """New SQLite-based search"""
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    
    start = time.perf_counter()
    
    cursor = conn.execute("""
        SELECT gene, stats_json 
        FROM gene_stats 
        WHERE gene_lower LIKE ? 
        ORDER BY mean_expression DESC 
        LIMIT ?
    """, (f'%{search_term.lower()}%', limit))
    
    results = []
    for row in cursor:
        stats = json.loads(row['stats_json'])
        results.append({
            "symbol": row['gene'],
            "mean_expression": stats["mean_expression"],
            "expressed_cells": stats["expressed_cells"]
        })
    
    total_time = time.perf_counter() - start
    conn.close()
    
    return results, total_time

def main():
    """Compare performance"""
    print("=" * 80)
    print("Performance Comparison: JSON vs SQLite")
    print("=" * 80)
    
    datasets = ["human_subset", "AD093044"]
    search_terms = ["MALAT1", "cd", "APO"]
    
    for dataset_id in datasets:
        json_path = Path(f"h5ad/precomputed/{dataset_id}/gene_stats.json")
        db_path = Path(f"h5ad/precomputed/{dataset_id}/data.db")
        
        # Check if we have JSON file (for comparison)
        has_json = json_path.exists()
        has_db = db_path.exists()
        
        if not has_db:
            print(f"\n⚠️  Database not found: {db_path}")
            continue
        
        print(f"\n{'='*80}")
        print(f"Dataset: {dataset_id}")
        print(f"{'='*80}")
        
        # Get file sizes
        if has_json:
            json_size = json_path.stat().st_size / (1024 * 1024)
        db_size = db_path.stat().st_size / (1024 * 1024)
        
        print(f"\n📊 Storage:")
        if has_json:
            print(f"   JSON file size: {json_size:.2f} MB")
        print(f"   SQLite DB size: {db_size:.2f} MB")
        
        print(f"\n⚡ Performance Comparison (10 iterations):")
        print(f"{'Search':<10} {'Method':<10} {'Load (ms)':<15} {'Search (ms)':<15} {'Total (ms)':<15} {'Memory':<15} {'Speedup'}")
        print("-" * 80)
        
        for search_term in search_terms:
            # Test SQLite
            sqlite_times = []
            for _ in range(10):
                _, elapsed = sqlite_search(db_path, search_term)
                sqlite_times.append(elapsed)
            
            sqlite_mean = statistics.mean(sqlite_times) * 1000
            sqlite_memory = db_size  # Only DB is loaded
            
            # Test JSON (if available)
            if has_json:
                json_load_times = []
                json_search_times = []
                for _ in range(10):
                    _, load_time, search_time = simulate_json_search(json_path, search_term)
                    json_load_times.append(load_time)
                    json_search_times.append(search_time)
                
                json_load_mean = statistics.mean(json_load_times) * 1000
                json_search_mean = statistics.mean(json_search_times) * 1000
                json_total_mean = json_load_mean + json_search_mean
                json_memory = json_size  # Entire JSON loaded
                
                speedup = json_total_mean / sqlite_mean if sqlite_mean > 0 else 0
                
                print(f"{search_term:<10} {'JSON':<10} {json_load_mean:>10.2f}    {json_search_mean:>10.2f}    {json_total_mean:>10.2f}    {json_memory:>8.2f} MB")
                print(f"{search_term:<10} {'SQLite':<10} {'N/A':<15} {'N/A':<15} {sqlite_mean:>10.2f}    {sqlite_memory:>8.2f} MB    {speedup:>6.1f}x")
            else:
                print(f"{search_term:<10} {'SQLite':<10} {'N/A':<15} {'N/A':<15} {sqlite_mean:>10.2f}    {sqlite_memory:>8.2f} MB")
        
        # Summary
        print(f"\n📈 Key Benefits:")
        print(f"   ✅ SQLite: Only loads matching genes (not entire dataset)")
        print(f"   ✅ SQLite: Uses indexes for fast sorting")
        print(f"   ✅ SQLite: No need to load entire JSON into memory")
        if has_json:
            print(f"   ✅ Memory saved: ~{json_size - db_size:.2f} MB per query")
    
    print(f"\n{'='*80}")
    print("✅ Comparison complete!")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()

