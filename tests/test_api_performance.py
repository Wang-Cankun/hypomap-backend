#!/usr/bin/env python3
"""
Test actual API endpoint performance
"""
import time
import requests
import statistics
from typing import Dict

def test_api_search(base_url: str, dataset_id: str, search_term: str, iterations: int = 10) -> Dict:
    """Test API endpoint performance"""
    url = f"{base_url}/api/v1/h5ad/{dataset_id}/genes/search"
    params = {"q": search_term, "limit": 50}
    
    times = []
    results_count = 0
    errors = 0
    
    for i in range(iterations):
        try:
            start = time.perf_counter()
            response = requests.get(url, params=params, timeout=10)
            elapsed = time.perf_counter() - start
            
            if response.status_code == 200:
                data = response.json()
                times.append(elapsed)
                if i == 0:
                    results_count = len(data)
            else:
                errors += 1
                print(f"  Error: {response.status_code} - {response.text[:100]}")
        except Exception as e:
            errors += 1
            print(f"  Exception: {e}")
    
    if not times:
        return {"error": "All requests failed"}
    
    return {
        "mean_time_ms": statistics.mean(times) * 1000,
        "median_time_ms": statistics.median(times) * 1000,
        "min_time_ms": min(times) * 1000,
        "max_time_ms": max(times) * 1000,
        "results_count": results_count,
        "errors": errors,
        "iterations": iterations
    }

def main():
    """Test API performance"""
    base_url = "http://localhost:9117/sskind-backend"
    
    print("=" * 70)
    print("API Endpoint Performance Test")
    print("=" * 70)
    
    # Check if server is running
    try:
        response = requests.get(f"{base_url}/api/v1/h5ad/datasets", timeout=5)
        if response.status_code != 200:
            print(f"❌ Server not responding correctly: {response.status_code}")
            return
    except Exception as e:
        print(f"❌ Cannot connect to server at {base_url}")
        print(f"   Error: {e}")
        print(f"\n💡 Make sure the server is running:")
        print(f"   /Users/wang.13246/.local/share/mamba/envs/sskind/bin/python main.py")
        return
    
    datasets = ["human_subset", "AD093044"]
    search_terms = ["MALAT1", "cd", "APO"]
    
    print(f"\n🔍 Testing API Endpoints (10 iterations each):")
    print(f"{'Dataset':<15} {'Search':<10} {'Mean (ms)':<12} {'Min (ms)':<12} {'Max (ms)':<12} {'Results':<10} {'Status'}")
    print("-" * 90)
    
    for dataset_id in datasets:
        for search_term in search_terms:
            result = test_api_search(base_url, dataset_id, search_term, iterations=10)
            
            if "error" in result:
                print(f"{dataset_id:<15} {search_term:<10} {'ERROR':<12}")
                continue
            
            mean_ms = result['mean_time_ms']
            min_ms = result['min_time_ms']
            max_ms = result['max_time_ms']
            count = result['results_count']
            
            # Performance rating
            if mean_ms < 10:
                status = "🚀 Excellent"
            elif mean_ms < 50:
                status = "✅ Very Fast"
            elif mean_ms < 100:
                status = "✅ Fast"
            elif mean_ms < 500:
                status = "⚠️  Acceptable"
            else:
                status = "❌ Slow"
            
            print(f"{dataset_id:<15} {search_term:<10} {mean_ms:>8.2f}    {min_ms:>8.2f}    {max_ms:>8.2f}    {count:>6}    {status}")
    
    print(f"\n{'='*70}")
    print("✅ API performance test complete!")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()

