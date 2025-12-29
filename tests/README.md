# Tests

This directory contains test scripts for the ssKIND backend.

## Performance Tests

### `test_sqlite_performance.py`
Benchmarks SQLite gene search performance across different datasets and search terms.

**Usage:**
```bash
python tests/test_sqlite_performance.py
```

**What it tests:**
- Search performance for various query types (exact match, prefix, contains)
- Database statistics and index usage
- Query execution plans
- Memory efficiency

### `test_api_performance.py`
Tests actual API endpoint performance.

**Usage:**
```bash
# Make sure server is running first
python tests/test_api_performance.py
```

**What it tests:**
- API endpoint response times
- End-to-end performance including network overhead
- Multiple datasets and search terms

### `test_performance_comparison.py`
Compares theoretical JSON-based approach vs SQLite approach.

**Usage:**
```bash
python tests/test_performance_comparison.py
```

**What it tests:**
- Performance comparison between JSON and SQLite
- Memory usage differences
- Speed improvements

## API Tests

### `test_spatial_api.py`
Tests spatial transcriptomics API endpoints.

**Usage:**
```bash
python tests/test_spatial_api.py
```

**What it tests:**
- Spatial dataset endpoints
- Spatial coordinates
- SVG (Spatially Variable Genes) endpoints
- Precomputed DEG endpoints

## Running All Tests

```bash
# Run all performance tests
python tests/test_sqlite_performance.py
python tests/test_api_performance.py
python tests/test_performance_comparison.py

# Run API tests
python tests/test_spatial_api.py
```

## Notes

- Performance tests require preprocessed datasets in `h5ad/precomputed/`
- API tests require the server to be running
- Some tests may take a few minutes to complete

