#!/usr/bin/env python3
"""
Simple test script to verify spatial API endpoints work
Run this after starting the server and importing spatial data
"""
import requests
import json

BASE_URL = "http://localhost:9117/sskind-backend/api/v1"

def test_endpoint(endpoint, description):
    """Test an API endpoint and print results"""
    print(f"\n🔍 Testing: {description}")
    print(f"URL: {BASE_URL}{endpoint}")
    
    try:
        response = requests.get(f"{BASE_URL}{endpoint}")
        print(f"Status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            if isinstance(data, list):
                print(f"Results: {len(data)} items")
                if data:
                    print(f"First item keys: {list(data[0].keys())}")
            else:
                print(f"Response: {json.dumps(data, indent=2)}")
        else:
            print(f"Error: {response.text}")
    except requests.exceptions.ConnectionError:
        print("❌ Connection failed - make sure the server is running")
    except Exception as e:
        print(f"❌ Error: {e}")

def main():
    print("🧪 Testing Spatial API Endpoints")
    print("=" * 50)
    
    # Test spatial datasets endpoints
    test_endpoint("/spatial-datasets/", "Get all spatial datasets")
    test_endpoint("/spatial-datasets/stats/", "Get spatial datasets statistics")
    test_endpoint("/spatial-datasets/search/disease/AD", "Search spatial datasets by disease (AD)")
    test_endpoint("/spatial-datasets/search/species/human", "Search spatial datasets by species (human)")
    test_endpoint("/spatial-datasets/search/methodology/Visium", "Search spatial datasets by methodology (Visium)")
    test_endpoint("/spatial-datasets/ST024001", "Get specific spatial dataset")
    
    # Test spatial papers endpoints
    test_endpoint("/spatial-papers/", "Get all spatial papers")
    test_endpoint("/spatial-papers/stats/", "Get spatial papers statistics")
    test_endpoint("/spatial-papers/search/disease/AD", "Search spatial papers by disease (AD)")
    test_endpoint("/spatial-papers/search/species/human", "Search spatial papers by species (human)")
    test_endpoint("/spatial-papers/search/methodology/Visium", "Search spatial papers by methodology (Visium)")
    test_endpoint("/spatial-papers/ST024", "Get specific spatial paper")
    
    print("\n✅ API testing completed!")

if __name__ == "__main__":
    main()
