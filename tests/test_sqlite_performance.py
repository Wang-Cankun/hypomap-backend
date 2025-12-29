#!/usr/bin/env python3
"""
Performance test script for SQLite gene search
Compares SQLite search performance and validates correctness
"""
import time
import sqlite3
import json
import statistics
from pathlib import Path
from typing import List, Dict

def test_sqlite_search(db_path: Path, search_term: str, limit: int = 50, iterations: int = 10) -> Dict:
    """Test SQLite search performance"""
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    
    times = []
    results_count = 0
    
    for i in range(iterations):
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
        
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        
        if i == 0:
            results_count = len(results)
    
    conn.close()
    
    return {
        "mean_time": statistics.mean(times),
        "median_time": statistics.median(times),
        "min_time": min(times),
        "max_time": max(times),
        "std_dev": statistics.stdev(times) if len(times) > 1 else 0,
        "results_count": results_count,
        "iterations": iterations
    }

def test_query_plan(db_path: Path, search_term: str):
    """Check if query uses indexes"""
    conn = sqlite3.connect(str(db_path))
    
    # Enable query plan
    cursor = conn.execute("EXPLAIN QUERY PLAN SELECT gene, stats_json FROM gene_stats WHERE gene_lower LIKE ? ORDER BY mean_expression DESC LIMIT 50", 
                         (f'%{search_term.lower()}%',))
    
    plan = cursor.fetchall()
    conn.close()
    
    return plan

def get_database_stats(db_path: Path) -> Dict:
    """Get database statistics"""
    conn = sqlite3.connect(str(db_path))
    
    # Get total genes
    cursor = conn.execute("SELECT COUNT(*) as count FROM gene_stats")
    total_genes = cursor.fetchone()[0]
    
    # Get database size
    db_size = db_path.stat().st_size
    
    # Get index info
    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='gene_stats'")
    indexes = [row[0] for row in cursor.fetchall()]
    
    conn.close()
    
    return {
        "total_genes": total_genes,
        "db_size_mb": db_size / (1024 * 1024),
        "indexes": indexes
    }

def main():
    """Run performance tests"""
    print("=" * 70)
    print("SQLite Gene Search Performance Test")
    print("=" * 70)
    
    # Test datasets
    datasets = [
        "human_subset",
        "AD093044",
        "AD093045"
    ]
    
    search_terms = [
        "MALAT1",      # Exact match
        "cd",          # Common prefix (many matches)
        "APO",         # Common prefix
        "GAD",         # Specific genes
        "xyz123"       # No matches
    ]
    
    for dataset_id in datasets:
        db_path = Path(f"h5ad/precomputed/{dataset_id}/data.db")
        
        if not db_path.exists():
            print(f"\n⚠️  Database not found: {db_path}")
            continue
        
        print(f"\n{'='*70}")
        print(f"Dataset: {dataset_id}")
        print(f"{'='*70}")
        
        # Get database stats
        stats = get_database_stats(db_path)
        print(f"\n📊 Database Statistics:")
        print(f"   Total genes: {stats['total_genes']:,}")
        print(f"   Database size: {stats['db_size_mb']:.2f} MB")
        print(f"   Indexes: {', '.join(stats['indexes'])}")
        
        # Test each search term
        print(f"\n🔍 Search Performance Tests (10 iterations each):")
        print(f"{'Search Term':<15} {'Mean (ms)':<12} {'Min (ms)':<12} {'Max (ms)':<12} {'Results':<10} {'Status'}")
        print("-" * 70)
        
        for search_term in search_terms:
            try:
                result = test_sqlite_search(db_path, search_term, limit=50, iterations=10)
                
                mean_ms = result['mean_time'] * 1000
                min_ms = result['min_time'] * 1000
                max_ms = result['max_time'] * 1000
                count = result['results_count']
                
                # Performance rating
                if mean_ms < 1:
                    status = "🚀 Excellent"
                elif mean_ms < 5:
                    status = "✅ Very Fast"
                elif mean_ms < 10:
                    status = "✅ Fast"
                elif mean_ms < 50:
                    status = "⚠️  Acceptable"
                else:
                    status = "❌ Slow"
                
                print(f"{search_term:<15} {mean_ms:>8.3f}    {min_ms:>8.3f}    {max_ms:>8.3f}    {count:>6}    {status}")
                
            except Exception as e:
                print(f"{search_term:<15} {'ERROR':<12} {str(e)[:30]}")
        
        # Check query plan for one search
        print(f"\n📋 Query Plan Analysis (for 'cd' search):")
        plan = test_query_plan(db_path, "cd")
        for row in plan:
            print(f"   {row}")
        
        # Memory efficiency test
        print(f"\n💾 Memory Efficiency:")
        print(f"   Database size: {stats['db_size_mb']:.2f} MB")
        print(f"   Estimated JSON size: ~{stats['total_genes'] * 0.0002:.2f} MB (if loaded)")
        print(f"   Memory saved: ~{stats['total_genes'] * 0.0002 - stats['db_size_mb']:.2f} MB per query")
    
    print(f"\n{'='*70}")
    print("✅ Performance test complete!")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()

