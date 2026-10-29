#!/usr/bin/env python3
"""
Median OrthoFinder score-normalisation parameters over the species pairs fitted with enough points.

  python fit_norm_params.py normfit.tsv [--min-hits 30]

normfit.tsv is written by the patched OrthoFinder when OF_NORM_LOG is set (columns: species i, species j, points used
in the fit, b, a). Prints "b,a", the value to give to OF_NORM_PARAMS.
"""
import argparse
import statistics


# Print the median parameters of the pairs with >= --min-hits points.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log")
    ap.add_argument("--min-hits", type=int, default=30)
    a = ap.parse_args()
    rows = [line.rstrip("\n").split("\t") for line in open(a.log) if line.strip()]
    sel = [r for r in rows if int(r[2]) >= a.min_hits]
    if not sel:
        raise SystemExit(f"no species pair with >= {a.min_hits} points in {a.log}")
    print(f"{statistics.median(float(r[3]) for r in sel):.3f},{statistics.median(float(r[4]) for r in sel):.3f}")


if __name__ == "__main__":
    main()
