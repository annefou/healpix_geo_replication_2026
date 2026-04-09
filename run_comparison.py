#!/usr/bin/env python3
"""
HEALPix-geo vs H3 Cross-DGGS Comparison

Reads results produced by the healpix-geo benchmark script and compares
against pre-computed H3 reference results from the original replication
study (https://github.com/annefou/dggs_replication_2026).

Default usage (uses bundled H3 reference data):

  python run_healpix_geo_replication.py --all --output results_healpix_geo
  python run_comparison.py \
      --healpix-geo results_healpix_geo \
      --output results_comparison

Or with custom H3 results:

  python run_comparison.py \
      --h3 /path/to/h3_results \
      --healpix-geo results_healpix_geo \
      --output results_comparison

Scientific questions answered:
  Does HEALPix (via healpix-geo) validate Law & Ardo (2024) like H3 does?
  Does the sphere/WGS84 choice affect performance or crossover point?

Repository: https://github.com/annefou/healpix_geo_replication_2026
Original H3 replication: https://github.com/annefou/dggs_replication_2026

Author: Anne Fouilloux
Date: 2026-04-09
"""

import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

CODE_VERSION = "2026-04-09-comparison-v1.0.0"

# Colours consistent across all plots
COLORS = {
    "h3":             "#E07B39",   # orange
    "healpix_sphere": "#5B8DB8",   # blue
    "healpix_wgs84":  "#8B5BA8",   # purple
    "vector":         "#AAAAAA",   # grey
    "raster":         "#A0C878",   # green
}
LABELS = {
    "h3":             "H3 (sphere)",
    "healpix_sphere": "HEALPix / sphere",
    "healpix_wgs84":  "HEALPix / WGS84",
    "vector":         "Vector overlay",
    "raster":         "Raster (numpy)",
}

# Default path to bundled H3 reference results
REFERENCE_H3_DIR = Path(__file__).parent / "reference_h3"


# =============================================================================
# Loaders
# =============================================================================

def load_h3(results_dir: Path) -> dict:
    """
    Load results from H3 replication (run_replication.py).

    Vector CSV columns:
        num_layers, dggs_total_time, vector_total_time, vector_success, ...
    Raster CSV columns:
        num_layers, raster_total, dggs_preindex_total, ...
    """
    out = {}
    vec_path = results_dir / "vector_benchmark.csv"
    if vec_path.exists():
        df = pd.read_csv(vec_path)
        out["vector"] = pd.DataFrame({
            "num_layers":  df["num_layers"],
            "dggs_total":  df["dggs_total_time"],
            "vector_total": df["vector_total_time"],
            "vector_success": df["vector_success"],
        })

    ras_path = results_dir / "raster_benchmark.csv"
    if ras_path.exists():
        df = pd.read_csv(ras_path)
        out["raster"] = pd.DataFrame({
            "num_layers":  df["num_layers"],
            "dggs_total":  df["dggs_preindex_total"],
            "raster_total": df["raster_total"],
        })
    return out


def load_healpix_geo(results_dir: Path) -> dict:
    """
    Load results from run_healpix_geo_replication.py (healpix-geo).

    Vector CSV columns:
        num_layers, healpix_sphere_total, healpix_wgs84_total,
        vector_total, vector_success, ...
    Raster CSV columns:
        num_layers, raster_total, healpix_sphere_total, healpix_wgs84_total
    """
    out = {}
    vec_path = results_dir / "vector_benchmark_healpix_geo.csv"
    if vec_path.exists():
        df = pd.read_csv(vec_path)
        out["vector_sphere"] = pd.DataFrame({
            "num_layers":    df["num_layers"],
            "dggs_total":    df["healpix_sphere_total"],
            "vector_total":  df["vector_total"],
            "vector_success": df["vector_success"],
        })
        out["vector_wgs84"] = pd.DataFrame({
            "num_layers":    df["num_layers"],
            "dggs_total":    df["healpix_wgs84_total"],
            "vector_total":  df["vector_total"],
            "vector_success": df["vector_success"],
        })

    ras_path = results_dir / "raster_benchmark_healpix_geo.csv"
    if ras_path.exists():
        df = pd.read_csv(ras_path)
        out["raster_sphere"] = pd.DataFrame({
            "num_layers":   df["num_layers"],
            "dggs_total":   df["healpix_sphere_total"],
            "raster_total": df["raster_total"],
        })
        out["raster_wgs84"] = pd.DataFrame({
            "num_layers":   df["num_layers"],
            "dggs_total":   df["healpix_wgs84_total"],
            "raster_total": df["raster_total"],
        })

    # Ellipsoid analysis
    ea_path = results_dir / "ellipsoid_analysis.json"
    if ea_path.exists():
        with open(ea_path) as f:
            out["ellipsoid_analysis"] = json.load(f)

    return out


# =============================================================================
# Crossover analysis
# =============================================================================

def find_crossover(layers: np.ndarray, dggs_times: np.ndarray,
                    vector_times: np.ndarray) -> Optional[float]:
    """
    Find the number of layers at which DGGS becomes faster than vector overlay.
    Returns None if DGGS is always faster or always slower in the tested range.
    """
    speedups = vector_times / dggs_times
    for i in range(len(speedups) - 1):
        if speedups[i] < 1.0 and speedups[i + 1] >= 1.0:
            x0, x1 = layers[i], layers[i + 1]
            y0, y1 = speedups[i], speedups[i + 1]
            return float(x0 + (1.0 - y0) * (x1 - x0) / (y1 - y0))
        if speedups[i] >= 1.0:
            return float(layers[i])
    return None


# =============================================================================
# Build unified comparison table
# =============================================================================

def build_comparison_table(h3: dict, geo: dict) -> pd.DataFrame:
    """
    Merge all vector results into a single wide table indexed by num_layers.
    """
    frames = {}
    if "vector" in h3:
        frames["h3"] = h3["vector"][["num_layers", "dggs_total", "vector_total"]]
    if "vector_sphere" in geo:
        frames["healpix_sphere"] = geo["vector_sphere"][["num_layers", "dggs_total"]]
    if "vector_wgs84" in geo:
        frames["healpix_wgs84"] = geo["vector_wgs84"][["num_layers", "dggs_total"]]

    if not frames:
        return pd.DataFrame()

    base_key = "h3" if "h3" in frames else next(iter(frames))
    result = frames[base_key].rename(columns={"dggs_total": "h3_total"})

    for key, df in frames.items():
        if key == base_key:
            continue
        col = f"{key}_total"
        result = result.merge(
            df.rename(columns={"dggs_total": col}),
            on="num_layers", how="outer"
        )

    result = result.sort_values("num_layers").reset_index(drop=True)

    # Speedup columns
    for col in ["h3_total", "healpix_sphere_total", "healpix_wgs84_total"]:
        sp_col = col.replace("_total", "_speedup")
        if col in result.columns and "vector_total" in result.columns:
            result[sp_col] = result["vector_total"] / result[col]

    return result


# =============================================================================
# Plotting
# =============================================================================

def plot_comparison(h3: dict, geo: dict,
                     table: pd.DataFrame, output_dir: Path):
    """Produce 6-panel comparison figure."""
    if not HAS_MATPLOTLIB:
        print("matplotlib not available -- skipping plots")
        return

    fig = plt.figure(figsize=(18, 12))
    fig.suptitle(
        "HEALPix-geo vs H3 -- Replication of Law & Ardo (2024)\n"
        "H3 (reference) vs HEALPix/sphere vs HEALPix/WGS84 | "
        "DOI: 10.1080/20964471.2024.2429847",
        fontsize=13, fontweight="bold",
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # -- Panel 1: Vector timing (log-log) --
    ax = fig.add_subplot(gs[0, 0])
    plotted_vector = False
    for method, src, key, color in [
        ("h3",             h3,  "vector",        COLORS["h3"]),
        ("healpix_sphere", geo, "vector_sphere",  COLORS["healpix_sphere"]),
        ("healpix_wgs84",  geo, "vector_wgs84",   COLORS["healpix_wgs84"]),
    ]:
        if key in src:
            df = src[key]
            valid = df[df["vector_success"] == True] if "vector_success" in df.columns else df
            ax.loglog(valid["num_layers"], valid["dggs_total"], "o-",
                      color=color, label=LABELS[method], linewidth=2, markersize=6)
            if not plotted_vector and "vector_total" in valid.columns:
                ax.loglog(valid["num_layers"], valid["vector_total"], "s--",
                          color=COLORS["vector"], label=LABELS["vector"],
                          linewidth=1.5, markersize=5, alpha=0.7)
                plotted_vector = True
    ax.set_xlabel("Number of layers")
    ax.set_ylabel("Time (s)")
    ax.set_title("Vector Benchmark -- Timing")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # -- Panel 2: Vector speedup vs vector overlay --
    ax = fig.add_subplot(gs[0, 1])
    for col, color, label in [
        ("h3_speedup",             COLORS["h3"],             LABELS["h3"]),
        ("healpix_sphere_speedup", COLORS["healpix_sphere"], LABELS["healpix_sphere"]),
        ("healpix_wgs84_speedup",  COLORS["healpix_wgs84"],  LABELS["healpix_wgs84"]),
    ]:
        if col in table.columns:
            valid = table.dropna(subset=[col])
            ax.semilogy(valid["num_layers"], valid[col], "o-",
                        color=color, label=label, linewidth=2, markersize=6)
    ax.axhline(y=1, color="gray", linestyle="--", linewidth=1, label="Break-even")
    ax.set_xlabel("Number of layers")
    ax.set_ylabel("Speedup vs vector overlay")
    ax.set_title("Vector Speedup")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # -- Panel 3: Raster timing --
    ax = fig.add_subplot(gs[0, 2])
    plotted_raster = False
    for method, src, key, color in [
        ("h3",             h3,  "raster",        COLORS["h3"]),
        ("healpix_sphere", geo, "raster_sphere",  COLORS["healpix_sphere"]),
        ("healpix_wgs84",  geo, "raster_wgs84",   COLORS["healpix_wgs84"]),
    ]:
        if key in src:
            df = src[key]
            ax.loglog(df["num_layers"], df["dggs_total"], "o-",
                      color=color, label=LABELS[method], linewidth=2, markersize=6)
            if not plotted_raster and "raster_total" in df.columns:
                ax.loglog(df["num_layers"], df["raster_total"], "s--",
                          color=COLORS["raster"], label=LABELS["raster"],
                          linewidth=1.5, markersize=5, alpha=0.7)
                plotted_raster = True
    ax.set_xlabel("Number of layers")
    ax.set_ylabel("Time (s)")
    ax.set_title("Raster Benchmark -- Timing")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # -- Panel 4: Crossover comparison bar chart --
    ax = fig.add_subplot(gs[1, 0])
    crossovers = {}
    for method, src, key, color in [
        ("h3",             h3,  "vector",        COLORS["h3"]),
        ("healpix_sphere", geo, "vector_sphere",  COLORS["healpix_sphere"]),
        ("healpix_wgs84",  geo, "vector_wgs84",   COLORS["healpix_wgs84"]),
    ]:
        if key in src:
            df = src[key].dropna(subset=["dggs_total", "vector_total"])
            if not df.empty:
                co = find_crossover(
                    df["num_layers"].values,
                    df["dggs_total"].values,
                    df["vector_total"].values,
                )
                crossovers[method] = co
    if crossovers:
        names = [LABELS[m] for m in crossovers]
        vals = [v if v is not None else 0 for v in crossovers.values()]
        colors = [COLORS[m] for m in crossovers]
        bars = ax.bar(range(len(names)), vals, color=colors, alpha=0.85, edgecolor="white")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=15, ha="right", fontsize=8)
        for bar, val in zip(bars, crossovers.values()):
            label = f"~{val:.0f}" if val else "N/A"
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.3,
                    label, ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_ylabel("Layers at crossover point")
    ax.set_title("Crossover: DGGS becomes\nfaster than vector overlay")
    ax.grid(True, alpha=0.3, axis="y")

    # -- Panel 5: Speedup at max tested layers --
    ax = fig.add_subplot(gs[1, 1])
    max_speedups = {}
    for col, method in [
        ("h3_speedup",             "h3"),
        ("healpix_sphere_speedup", "healpix_sphere"),
        ("healpix_wgs84_speedup",  "healpix_wgs84"),
    ]:
        if col in table.columns:
            valid = table.dropna(subset=[col])
            if not valid.empty:
                max_speedups[method] = valid[col].max()
    if max_speedups:
        names = [LABELS[m] for m in max_speedups]
        vals = list(max_speedups.values())
        colors = [COLORS[m] for m in max_speedups]
        bars = ax.bar(range(len(names)), vals, color=colors, alpha=0.85, edgecolor="white")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=15, ha="right", fontsize=8)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{val:.0f}x", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
    ax.set_ylabel("Max speedup vs vector overlay")
    ax.set_title("Peak Speedup (all tested layers)")
    ax.grid(True, alpha=0.3, axis="y")

    # -- Panel 6: Ellipsoid pixel-assignment difference --
    ax = fig.add_subplot(gs[1, 2])
    ea = geo.get("ellipsoid_analysis", {})
    if ea:
        regions = list(ea.values())
        names = [r["region"].replace("_", "\n") for r in regions]
        diffs = [r["pixels_different_pct"] for r in regions]
        lats  = [r["center_lat"] for r in regions]
        bars  = ax.bar(names, diffs, color=COLORS["healpix_wgs84"],
                       alpha=0.85, edgecolor="white")
        for bar, d in zip(bars, diffs):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{d:.0f}%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
        ax.set_ylim(0, 110)
        ax.set_ylabel("Pixels in different cell (%)")
        ax.set_title("Sphere vs WGS84 Indexing\nDifference by Region")
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(range(len(lats)))
        ax2.set_xticklabels([f"{l:+.0f}" for l in lats],
                             fontsize=8, color="gray")
        ax2.set_xlabel("Center latitude", fontsize=8, color="gray")
    else:
        ax.text(0.5, 0.5, "No ellipsoid analysis\ndata available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="gray")
        ax.set_title("Sphere vs WGS84 Indexing\nDifference by Region")
    ax.grid(True, alpha=0.3, axis="y")

    plt.savefig(output_dir / "comparison.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "comparison.pdf", bbox_inches="tight")
    print(f"Plots saved: comparison.png / .pdf")
    plt.close()


# =============================================================================
# Summary
# =============================================================================

def generate_summary(table: pd.DataFrame, h3: dict,
                      geo: dict, output_dir: Path) -> dict:
    """Print summary table and save JSON."""

    summary = {
        "generated": datetime.now().isoformat(),
        "code_version": CODE_VERSION,
        "paper": {"doi": "10.1080/20964471.2024.2429847"},
        "h3_reference": "https://github.com/annefou/dggs_replication_2026",
        "methods": {},
    }

    print("\n" + "=" * 70)
    print("CROSS-DGGS COMPARISON SUMMARY")
    print("HEALPix-geo vs H3 (reference)")
    print("=" * 70)
    print(f"\n{'Method':<22} {'Max speedup':>12} {'Crossover':>12} {'Claim validated':>16}")
    print("-" * 64)

    for method, sp_col, src, key in [
        ("H3 (reference)",   "h3_speedup",             h3,  "vector"),
        ("HEALPix/sphere",   "healpix_sphere_speedup",  geo, "vector_sphere"),
        ("HEALPix/WGS84",    "healpix_wgs84_speedup",   geo, "vector_wgs84"),
    ]:
        max_sp = "--"
        crossover = "--"
        validated = "--"

        if sp_col in table.columns:
            valid_sp = table.dropna(subset=[sp_col])
            if not valid_sp.empty:
                ms = valid_sp[sp_col].max()
                max_sp = f"{ms:.0f}x"
                validated = "YES" if ms > 1.0 else "NO"
                summary["methods"][method] = {"max_speedup": round(ms, 1)}

        if key in src:
            df = src[key].dropna(subset=["dggs_total", "vector_total"])
            if not df.empty:
                co = find_crossover(
                    df["num_layers"].values,
                    df["dggs_total"].values,
                    df["vector_total"].values,
                )
                crossover = f"~{co:.0f} layers" if co else "not reached"
                if method in summary["methods"]:
                    summary["methods"][method]["crossover_layers"] = co

        print(f"{method:<22} {max_sp:>12} {crossover:>12} {validated:>16}")

    # Ellipsoid impact
    ea = geo.get("ellipsoid_analysis", {})
    if ea:
        print("\n  Sphere vs WGS84 pixel-assignment difference:")
        print(f"  {'Region':<16} {'Center lat':>12} {'Pixels differ':>14} {'Jaccard':>10}")
        print("  " + "-" * 54)
        for name, v in ea.items():
            jac = v.get("polygon_jaccard_similarity")
            jac_str = f"{jac:.4f}" if jac is not None else "--"
            print(f"  {name:<16} {v['center_lat']:>+10.0f}   "
                  f"{v['pixels_different_pct']:>11.1f}%  {jac_str:>10}")
        summary["ellipsoid_analysis"] = {
            name: {
                "center_lat": v["center_lat"],
                "pixels_different_pct": v["pixels_different_pct"],
                "jaccard": v.get("polygon_jaccard_similarity"),
            }
            for name, v in ea.items()
        }

    # Save
    with open(output_dir / "comparison_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Save comparison table
    if not table.empty:
        table.to_csv(output_dir / "comparison_table.csv", index=False)
        print(f"\n  Comparison table -> comparison_table.csv")

    print(f"\n  Full summary     -> comparison_summary.json")
    return summary


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="HEALPix-geo vs H3 comparison: cross-DGGS benchmark analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Run the healpix-geo benchmark first, then compare against H3 reference:

  python run_healpix_geo_replication.py --all --output results_healpix_geo
  python run_comparison.py --healpix-geo results_healpix_geo

H3 reference results are bundled in reference_h3/. To use custom H3 results:

  python run_comparison.py \\
      --h3 /path/to/h3_results \\
      --healpix-geo results_healpix_geo
        """,
    )
    parser.add_argument("--h3", type=Path, default=None,
                        help=f"Directory with H3 results (default: {REFERENCE_H3_DIR})")
    parser.add_argument("--healpix-geo", type=Path, default=None,
                        help="Directory with run_healpix_geo_replication.py results")
    parser.add_argument("--output", "-o", type=Path, default=Path("results_comparison"),
                        help="Output directory (default: results_comparison)")
    args = parser.parse_args()

    # Default H3 to bundled reference data
    if args.h3 is None:
        args.h3 = REFERENCE_H3_DIR

    # Require healpix-geo results
    if args.healpix_geo is None:
        parser.error("--healpix-geo is required (run run_healpix_geo_replication.py first)")

    args.output.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("HEALPIX-GEO vs H3 COMPARISON")
    print("=" * 70)
    print(f"Version: {CODE_VERSION}")
    print(f"Output:  {args.output}")

    # Load results
    h3_data, geo_data = {}, {}

    if args.h3.exists():
        h3_data = load_h3(args.h3)
        source = "bundled reference" if args.h3 == REFERENCE_H3_DIR else str(args.h3)
        print(f"H3            : loaded from {source}  "
              f"(vector={'vector' in h3_data}, raster={'raster' in h3_data})")
    else:
        print(f"WARNING: H3 path not found: {args.h3}")

    if args.healpix_geo.exists():
        geo_data = load_healpix_geo(args.healpix_geo)
        print(f"HEALPix-geo   : loaded from {args.healpix_geo}  "
              f"(sphere={'vector_sphere' in geo_data}, "
              f"wgs84={'vector_wgs84' in geo_data}, "
              f"ellipsoid_analysis={'ellipsoid_analysis' in geo_data})")
    else:
        print(f"WARNING: healpix-geo path not found: {args.healpix_geo}")

    # Build unified table
    table = build_comparison_table(h3_data, geo_data)

    # Plot + summary
    plot_comparison(h3_data, geo_data, table, args.output)
    generate_summary(table, h3_data, geo_data, args.output)

    print(f"\nResults saved to: {args.output}/")
    for f in sorted(args.output.iterdir()):
        print(f"   {f.name}")


if __name__ == "__main__":
    main()
