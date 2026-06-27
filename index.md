# HEALPix-geo benchmark replication

A **replication study** of the discrete-global-grid benchmarks from Law & Ardo (2024),
*"Using a discrete global grid system for a scalable, interoperable, and reproducible
system of land-use mapping"* (*Big Earth Data*,
[doi:10.1080/20964471.2024.2429847](https://doi.org/10.1080/20964471.2024.2429847)),
re-run on **HEALPix** via [`healpix-geo`](https://github.com/GRID4EARTH/healpix-geo) —
with a **sphere vs WGS84 ellipsoid** comparison the original H3-based study could not make.

It is the HEALPix companion to the [H3 replication](https://github.com/annefou/dggs_replication_2026):
same paper, same benchmarks, but indexed on the GRID4EARTH ellipsoidal-HEALPix stack so
the effect of the reference surface becomes measurable.

## What the replication tests

The paper's central claim is that indexing vector geometries onto a discrete global grid
makes raster–vector land-use operations dramatically faster while staying reproducible.
We reproduce that benchmark with HEALPix **cells** (addressed by `depth`, the GRID4EARTH
term for the refinement level), on both a bare sphere and the WGS84 ellipsoid.

## Key findings

### HEALPix-geo matches H3 performance

All three DGGS implementations validate the paper's headline result with consistent numbers
(all at **HEALPix depth 9**, matching H3 resolution 9 for a like-for-like comparison):

| Method | Max speedup vs vector | Crossover point |
|---|---|---|
| H3 (sphere, reference) | ~5,800× | ~5 layers |
| HEALPix / sphere | ~5,691× | ~5 layers |
| HEALPix / WGS84 | ~5,603× | ~5 layers |

### The WGS84 ellipsoid matters for European EO data

Because HEALPix on the sphere and on the WGS84 ellipsoid place a latitude-dependent share
of points in a **different cell**, the reference surface is not a free choice for
Copernicus-resolution analysis:

| Region | Centre latitude | Points in a different cell |
|---|---|---|
| Equatorial | 0° | low |
| Mid-latitude (Mediterranean) | +48° | the large majority |

The notebook [**Why the reference surface matters**](notebooks/01_wgs84_matters.ipynb)
reproduces this pattern from scratch, in seconds, with `healpix-geo`.

## Reproduce it

- **Source & Docker:** [`annefou/healpix_geo_replication_2026`](https://github.com/annefou/healpix_geo_replication_2026)
- **Archived release:** [doi:10.5281/zenodo.19488374](https://doi.org/10.5281/zenodo.19488374)
- **Core library:** [`healpix-geo`](https://github.com/GRID4EARTH/healpix-geo) (GRID4EARTH) — HEALPix on the WGS84 ellipsoid
