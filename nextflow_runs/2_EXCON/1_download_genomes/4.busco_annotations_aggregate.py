#!/usr/bin/env python3
import argparse
import os
import re
import glob
import csv


def parse_short_summary(path):
    """Parse BUSCO short_summary.txt and return a dict of metrics."""
    with open(path, "r") as fh:
        lines = fh.readlines()

    busco_version = None
    lineage = None
    mode = None
    input_file = None
    C_pct = S_pct = D_pct = F_pct = M_pct = None
    n_total = None
    n_C = n_S = n_D = n_F = n_M = None

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("# BUSCO version is"):
            busco_version = line.split(":", 1)[1].strip()

        elif line.startswith("# The lineage dataset is"):
            tmp = line.split(":", 1)[1].strip()
            lineage = tmp.split()[0]

        elif line.startswith("# Summarized benchmarking in BUSCO notation for file"):
            parts = line.split("for file", 1)
            if len(parts) == 2:
                input_file = parts[1].strip()

        elif line.startswith("# BUSCO was run in mode"):
            mode = line.split(":", 1)[1].strip()

        elif line.startswith("C:") and "n:" in line:
            m = re.match(
                r"C:(\d+(?:\.\d+)?)%\[S:(\d+(?:\.\d+)?)%,D:(\d+(?:\.\d+)?)%\],"
                r"F:(\d+(?:\.\d+)?)%,M:(\d+(?:\.\d+)?)%,n:(\d+)",
                line,
            )
            if m:
                C_pct, S_pct, D_pct, F_pct, M_pct, n_total = m.groups()

        elif "Complete BUSCOs (C)" in line:
            n_C = int(line.split()[0])

        elif "Complete and single-copy BUSCOs (S)" in line:
            n_S = int(line.split()[0])

        elif "Complete and duplicated BUSCOs (D)" in line:
            n_D = int(line.split()[0])

        elif "Fragmented BUSCOs (F)" in line:
            n_F = int(line.split()[0])

        elif "Missing BUSCOs (M)" in line:
            n_M = int(line.split()[0])

        elif "Total BUSCO groups searched" in line and n_total is None:
            n_total = line.split()[0]

    return {
        "busco_version": busco_version,
        "lineage": lineage,
        "mode": mode,
        "input_file": input_file,
        "C_pct": C_pct,
        "S_pct": S_pct,
        "D_pct": D_pct,
        "F_pct": F_pct,
        "M_pct": M_pct,
        "n_total": n_total,
        "n_C": n_C,
        "n_S": n_S,
        "n_D": n_D,
        "n_F": n_F,
        "n_M": n_M,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Aggregate BUSCO annotation runs into a single TSV."
    )
    ap.add_argument(
        "--busco_dir",
        default="busco_annotation_all",
        help="Directory containing BUSCO runs (each run in its own subfolder).",
    )
    ap.add_argument(
        "--out",
        default="000_BUSCO_annotations_summary.tsv",
        help="Output TSV file.",
    )
    args = ap.parse_args()

    rows = []

    for run_name in sorted(os.listdir(args.busco_dir)):
        run_path = os.path.join(args.busco_dir, run_name)
        if not os.path.isdir(run_path):
            continue

        # Infer species + source from run_name
        species = run_name
        source = "unknown"
        if run_name.endswith("_ENSEMBL_prot_long"):
            species = run_name[: -len("_ENSEMBL_prot_long")]
            source = "ENSEMBL"
        elif run_name.endswith("_NCBI_prot_long"):
            species = run_name[: -len("_NCBI_prot_long")]
            source = "NCBI"

        # Find short_summary + (optionally) full_table
        short_candidates = glob.glob(
            os.path.join(run_path, "run_*", "short_summary*.txt")
        )
        if not short_candidates:
            # Nothing useful here
            continue
        short_path = short_candidates[0]

        full_candidates = glob.glob(
            os.path.join(run_path, "run_*", "full_table*.tsv")
        )
        full_path = full_candidates[0] if full_candidates else ""

        metrics = parse_short_summary(short_path)

        row = {
            "run_name": run_name,
            "species": species,
            "source": source,
            "short_summary": short_path,
            "full_table": full_path,
            "busco_version": metrics["busco_version"],
            "lineage": metrics["lineage"],
            "mode": metrics["mode"],
            "input_file": metrics["input_file"],
            "busco_C_pct": metrics["C_pct"],
            "busco_S_pct": metrics["S_pct"],
            "busco_D_pct": metrics["D_pct"],
            "busco_F_pct": metrics["F_pct"],
            "busco_M_pct": metrics["M_pct"],
            "busco_n": metrics["n_total"],
            "n_complete": metrics["n_C"],
            "n_single": metrics["n_S"],
            "n_duplicated": metrics["n_D"],
            "n_fragmented": metrics["n_F"],
            "n_missing": metrics["n_M"],
        }
        rows.append(row)

    # Write TSV
    fieldnames = [
        "run_name",
        "species",
        "source",
        "busco_version",
        "lineage",
        "mode",
        "input_file",
        "busco_C_pct",
        "busco_S_pct",
        "busco_D_pct",
        "busco_F_pct",
        "busco_M_pct",
        "busco_n",
        "n_complete",
        "n_single",
        "n_duplicated",
        "n_fragmented",
        "n_missing",
        "short_summary",
        "full_table",
    ]

    with open(args.out, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()

