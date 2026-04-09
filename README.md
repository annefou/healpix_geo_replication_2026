# HEALPix-geo Benchmark Replication

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19488374.svg)](https://doi.org/10.5281/zenodo.19488374)
![Build Docker Image](https://github.com/annefou/healpix_geo_replication_2026/actions/workflows/docker-build.yml/badge.svg)
![Run Replication](https://github.com/annefou/healpix_geo_replication_2026/actions/workflows/run-replication.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Replication study** of the DGGS benchmarks from Law & Ardo (2024) using HEALPix indexing via [`healpix-geo`](https://github.com/EOPF-DGGS/healpix-geo), with sphere vs. WGS84 ellipsoid comparison.

This repository provides a **reproducible environment** for replicating the benchmarks from:

> Law, R.M. & Ardo, J. (2024). "Using a discrete global grid system for a scalable, interoperable, and reproducible system of land-use mapping"
> *Big Earth Data*, DOI: [10.1080/20964471.2024.2429847](https://doi.org/10.1080/20964471.2024.2429847)

**Companion to:** The [H3 replication study](https://github.com/annefou/dggs_replication_2026) which validated the original paper using H3. This repository extends the replication to HEALPix using `healpix-geo`, adding WGS84 ellipsoid support.

---

## Key Findings

### HEALPix-geo matches H3 performance

All three DGGS implementations (H3, HEALPix/sphere, HEALPix/WGS84) validate the paper's claims with consistent results:

| Method | Max speedup vs vector | Crossover point |
|--------|----------------------|-----------------|
| H3 (sphere, reference) | ~5,800x | ~5 layers |
| HEALPix / sphere | ~5,691x | ~5 layers |
| HEALPix / WGS84 | ~5,603x | ~5 layers |

> All benchmarks use **HEALPix depth 9** (~1 km^2 cells), chosen to match H3 resolution 9 for a like-for-like comparison.

### WGS84 ellipsoid matters for European EO data

A key finding is the impact of the reference surface on cell assignment:

| Region | Center latitude | Pixels in different cell | Jaccard similarity |
|--------|----------------|--------------------------|-------------------|
| Equatorial | 0 deg | 27% | 0.9951 |
| Mid-latitude (Mediterranean) | +48 deg | **98%** | 0.9843 |
| High-latitude (Scandinavia) | +62 deg | **91%** | 0.9868 |
| Arctic | +78 deg | 53% | 0.9908 |

**Key insight:** For European EO data (Copernicus/Sentinel, 45-65 deg N), sphere-based HEALPix indexing assigns almost every pixel to the wrong cell. WGS84 indexing via `healpix-geo` is strongly recommended for production workflows.

---

## Purpose

This replication study provides:

1. **HEALPix benchmarks** — same methodology as Law & Ardo (2024), using `healpix-geo` instead of H3
2. **Sphere vs WGS84 comparison** — quantifying geodetic differences across latitude bands
3. **Cross-DGGS comparison** — comparing HEALPix-geo against H3 reference results
4. **Docker container** with all dependencies pinned
5. **GitHub Actions** for automated CI/CD and continuous verification

---

## Scripts

| Script | Description |
|--------|-------------|
| `run_healpix_geo_replication.py` | HEALPix benchmarks via healpix-geo (sphere + WGS84) |
| `run_comparison.py` | Cross-DGGS comparison against H3 reference results |

---

## Quick Start

### Option 1: Docker (Recommended)

```bash
# Pull pre-built image from GitHub Container Registry
docker pull ghcr.io/annefou/healpix_geo_replication_2026:latest

# Run benchmarks + comparison
docker run -v $(pwd)/results:/app/results ghcr.io/annefou/healpix_geo_replication_2026:latest

# Or build locally
docker build -t healpix-geo-replication .
docker run -v $(pwd)/results:/app/results healpix-geo-replication
```

### Option 2: Local Python Environment

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run healpix-geo benchmarks
python run_healpix_geo_replication.py --all --output results

# Run comparison against H3 reference
python run_comparison.py --healpix-geo results --output results
```

### Option 3: Make

```bash
make local-setup
make run
```

---

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `VECTOR_LAYERS` | `5,10,20,50,100` | Comma-separated layer counts for vector benchmark |
| `RASTER_LAYERS` | `10,50,100,500` | Comma-separated layer counts for raster benchmark |
| `HEALPIX_DEPTH` | `9` | HEALPix depth for all benchmarks |
| `RANDOM_SEED` | `42` | Random seed for reproducibility |

### CLI Arguments

```bash
# Full benchmark suite
python run_healpix_geo_replication.py --all

# Ellipsoid analysis only
python run_healpix_geo_replication.py --ellipsoid-only

# Custom parameters
python run_healpix_geo_replication.py --all \
    --vector-layers 5,10,20,50 \
    --healpix-depth 9

# Comparison with custom H3 results
python run_comparison.py \
    --h3 /path/to/custom_h3_results \
    --healpix-geo results
```

---

## H3 Reference Data

Pre-computed H3 benchmark results are bundled in `reference_h3/` from the
[original H3 replication study](https://github.com/annefou/dggs_replication_2026) (v2.0.0).
The comparison script uses these by default.

---

## Output Files

```
results/
├── system_info.json                    # Hardware/software environment
├── vector_benchmark_healpix_geo.csv    # HEALPix sphere+WGS84 vector timings
├── raster_benchmark_healpix_geo.csv    # HEALPix sphere+WGS84 raster timings
├── ellipsoid_analysis.json             # Sphere vs WGS84 indexing difference
├── indexing_healpix_geo.json           # Indexing metadata
├── summary_healpix_geo.json            # HEALPix structured summary
├── benchmark_healpix_geo.png           # HEALPix benchmark plots
├── benchmark_healpix_geo.pdf           # HEALPix benchmark plots (PDF)
├── comparison_table.csv                # Cross-DGGS unified table
├── comparison_summary.json             # Cross-DGGS structured summary
├── comparison.png                      # Cross-DGGS comparison plot
└── comparison.pdf                      # Cross-DGGS comparison plot (PDF)
```

---

## Methodology

### Vector Benchmark (Figure 6)

Following the paper's Section 3.2.1 methodology:

1. **Data Generation**: Random points -> Voronoi polygons (`scipy.spatial.Voronoi`)
2. **Values**: Each polygon assigned 0 or 1 randomly
3. **Dissolve**: Polygons dissolved by value before overlay (as per paper)
4. **Traditional Method**: Unary union (spatial overlay) of all dissolved layers
5. **DGGS Method**: **Polyfill** polygons to HEALPix cells -> join on cell ID
6. **Classification**: 7 functions (prime, perfect, triangular, square, pentagonal, hexagonal, Fibonacci) -> 7-bit class

### Raster Benchmark (Figure 7)

1. **Data Generation**: Spatially-correlated rasters (Gaussian smoothing)
2. **Traditional Method**: NumPy array stacking and classification
3. **DGGS Method**: Index raster cells to HEALPix -> aggregate -> classify

### Ellipsoid Analysis

HEALPix benchmarks are run twice -- once with `ellipsoid='sphere'` and once with `ellipsoid='WGS84'` -- using `healpix-geo`. The analysis measures the percentage of pixels assigned to a different cell depending on the reference surface, across four latitude bands.

### Reproduction vs Replication

| Term | Definition | Implementation |
|------|------------|----------------|
| **Reproduction** | Same methodology, same tools | H3 (in companion repo) |
| **Replication** | Same methodology, different tools | healpix-geo for HEALPix sphere+WGS84 |

---

## GitHub Actions

### Triggering a Replication Run

#### Via GitHub UI

1. Go to **Actions** -> **Run Replication**
2. Click **Run workflow**
3. Choose options:
   - `benchmark_type`: `ci-test`, `quick-test`, or `full`
   - `random_seed`: For reproducibility (default: 42)

#### Via GitHub CLI

```bash
gh workflow run run-replication.yml -f benchmark_type=ci-test
```

### Scheduled Runs

The replication runs automatically every Sunday at 00:00 UTC.

---

## Related Work

- **Original paper**: Law & Ardo (2024), DOI: [10.1080/20964471.2024.2429847](https://doi.org/10.1080/20964471.2024.2429847)
- **Original benchmark code**: [dggsBenchmarks v1.1.1](https://github.com/manaakiwhenua/dggsBenchmarks/releases/tag/v1.1.1)
- **H3 replication study**: [annefou/dggs_replication_2026](https://github.com/annefou/dggs_replication_2026)
- **healpix-geo library**: [EOPF-DGGS/healpix-geo](https://github.com/EOPF-DGGS/healpix-geo)

---

## Citation

If you use this replication environment, please cite both:

```bibtex
@article{law2024dggs,
  title={Using a discrete global grid system for a scalable, interoperable,
         and reproducible system of land-use mapping},
  author={Law, Richard M and Ardo, James},
  journal={Big Earth Data},
  volume={9},
  number={1},
  pages={29--46},
  year={2024},
  doi={10.1080/20964471.2024.2429847}
}

@software{healpix_geo_replication_2026,
  title={HEALPix-geo Benchmark Replication},
  author={Fouilloux, Anne},
  year={2026},
  publisher={Zenodo},
  doi={10.5281/zenodo.19488374},
  url={https://doi.org/10.5281/zenodo.19488374},
  note={HEALPix replication of Law \& Ardo (2024) with WGS84 ellipsoid extension}
}
```

---

## License

This replication code is released under the MIT License.

---

## Contact

- **Author:** Anne Fouilloux (ORCID: [0000-0002-1784-2920](https://orcid.org/0000-0002-1784-2920))
- **Affiliation:** LifeWatch ERIC

---

## Acknowledgments

- Original research by Richard M. Law and James Ardo at Manaaki Whenua -- Landcare Research
- [healpix-geo](https://github.com/EOPF-DGGS/healpix-geo) library for WGS84-aware HEALPix indexing
- This replication follows the framework from the [Replication Handbook](https://forrt.org/replication_handbook/)
