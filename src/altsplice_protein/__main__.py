from __future__ import annotations
import argparse
from .pipeline import run


def main(argv: list[str] | None = None) -> dict:
    p = argparse.ArgumentParser(
        prog="altsplice_protein",
        description="Resolve before/after protein sequences for significant "
                    "SpliceSeq events.",
    )
    p.add_argument("csv", help="Path to spliceseq_info_*.csv")
    p.add_argument("--data-dir", default="data",
                   help="Where Ensembl files are downloaded/cached")
    p.add_argument("--results-dir", default="results",
                   help="Where output files are written")
    p.add_argument("--fdr-max", type=float, default=0.05)
    p.add_argument("--dpsi-min", type=float, default=0.1)
    p.add_argument("--no-rest", action="store_true",
                   help="Disable Ensembl REST fallback for unresolved names")
    a = p.parse_args(argv)
    stats = run(
        a.csv, a.data_dir, a.results_dir, a.fdr_max, a.dpsi_min, not a.no_rest
    )
    for k, v in stats.items():
        print(f"{k}: {v}")
    return stats


if __name__ == "__main__":
    main()
