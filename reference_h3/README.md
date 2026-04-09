# H3 Reference Results

Pre-computed H3 benchmark results from the original replication study:
https://github.com/annefou/dggs_replication_2026 (v2.0.0)

These results serve as a reference baseline for comparing HEALPix-geo
performance against H3. They were generated using:

- H3 resolution 9 (~1 km^2 cells)
- Random seed 42
- Vector layers: 5, 10, 20, 50, 100
- Raster layers: 10, 50, 100, 500

To regenerate these results with the original repo:

```bash
git clone https://github.com/annefou/dggs_replication_2026
cd dggs_replication_2026
pip install -r requirements.txt
python run_replication.py --all --output results_h3
```
