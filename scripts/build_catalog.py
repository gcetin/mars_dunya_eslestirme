"""Build the Earth terrain catalog.

Examples:
    python scripts/build_catalog.py --region 36 42 26 45 --output data/catalog/turkiye.parquet
    python scripts/build_catalog.py --global --output data/catalog/earth_global.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analog_finder.catalog import GLOBAL_REGIONS, build_catalog


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--region", nargs=4, type=float, action="append", metavar=("SOUTH", "NORTH", "WEST", "EAST"))
    parser.add_argument("--global", dest="whole_earth", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--window-km", type=float, default=60.0)
    parser.add_argument("--stride-km", type=float, default=20.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--threads", type=int, default=12)
    args = parser.parse_args()

    if args.whole_earth == bool(args.region):
        parser.error("Give either --global or at least one --region.")

    build_catalog(
        GLOBAL_REGIONS if args.whole_earth else [tuple(r) for r in args.region],
        args.output,
        window_m=args.window_km * 1000,
        stride_m=args.stride_km * 1000,
        workers=args.workers,
        threads=args.threads,
        log=lambda message: print(message, flush=True),
    )


if __name__ == "__main__":
    main()
