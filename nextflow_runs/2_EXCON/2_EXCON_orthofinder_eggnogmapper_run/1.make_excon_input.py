#!/usr/bin/env python3
"""Build the three-column EXCON input CSV from genome and BUSCO summaries."""

import argparse
import csv
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DOWNLOAD_DIR = SCRIPT_DIR.parent / "1_download_genomes"


def resolve_data_path(value: str, summary_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = summary_dir / path
    return path.resolve()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create EXCON Name,genome_path,gff_path input rows."
    )
    parser.add_argument(
        "--busco_tsv",
        type=Path,
        default=DOWNLOAD_DIR / "000_BUSCO_annotations_summary.tsv",
        help="BUSCO summary defining the retained species.",
    )
    parser.add_argument(
        "--summary_tsv",
        type=Path,
        default=DOWNLOAD_DIR / "genomes_out" / "download_summary_all.tsv",
        help="Genome-download summary containing FASTA and GFF paths.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=SCRIPT_DIR / "input_excon.csv",
        help="Output EXCON CSV (no header).",
    )
    args = parser.parse_args()

    busco_tsv = args.busco_tsv.expanduser().resolve()
    summary_tsv = args.summary_tsv.expanduser().resolve()
    out_csv = args.out.expanduser().resolve()

    keep = set()
    with busco_tsv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species = (row.get("species") or "").strip()
            if species:
                keep.add(species)
    print(f"[INFO] species from BUSCO summary: {len(keep)}")

    mapping = {}
    with summary_tsv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species_raw = (row.get("Species") or "").strip()
            source = (row.get("Source") or "").strip().lower()
            fasta_raw = (row.get("FASTA") or "").strip()
            gff_raw = (row.get("GFF") or "").strip()
            if not species_raw or not fasta_raw or not gff_raw:
                continue

            species = species_raw.replace(" ", "_")
            fasta = resolve_data_path(fasta_raw, summary_tsv.parent)
            gff = resolve_data_path(gff_raw, summary_tsv.parent)
            if not fasta.is_file():
                print(f"[WARN] FASTA path does not exist, skipping: {fasta}")
                continue
            if not gff.is_file():
                print(f"[WARN] GFF path does not exist, skipping: {gff}")
                continue

            record = (str(fasta), str(gff), source)
            if source == "ncbi" or species not in mapping:
                mapping[species] = record

    print(f"[INFO] species with usable genome paths: {len(mapping)}")

    rows = []
    for species in sorted(keep):
        if species not in mapping:
            print(f"[WARN] retained species has no valid genome/GFF pair: {species}")
            continue
        fasta, gff, _source = mapping[species]
        rows.append((species, fasta, gff))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(rows)
    print(f"[OK] wrote {len(rows)} rows to {out_csv}")


if __name__ == "__main__":
    main()
