import argparse
from pathlib import Path
from .DNA_clustering import DNA_clustering

def main():
    parser = argparse.ArgumentParser(description="DNA Fragment Grouper")
    parser.add_argument("csv", help="Input CSV file")
    parser.add_argument("--aggressive", action="store_true",
                        help="Combine singleton groups")
    args = parser.parse_args()

    out = DNA_clustering(Path(args.csv), aggressive=args.aggressive)
    print(f"Output written to {out}")