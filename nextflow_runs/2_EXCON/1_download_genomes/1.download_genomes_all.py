#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""download_genomes_all.py

Unified genome downloader for:
  - Ensembl-only species
  - NCBI-only species
  - Species where you want BOTH Ensembl + NCBI

Works with three tables:
  - Genomes_selection_ensembl.csv
  - Genomes_selection_ncbi.csv
  - Genomes_selection_both.csv

Folder layout (for --out genomes_out):

    genomes_out/
      ensembl/
        Acromyrmex_echinatior/
          <Ensembl FASTA/GFF>
      ncbi/
        Acromyrmex_echinatior/
          <NCBI FASTA/GFF>
      both/
        Bombus_pascuorum/
          ensembl_softmasked.fa or ensembl_genome.fa
          ensembl_annotation.gff3 / .gtf
          ncbi_GCF_905332965.1_iyBomPasc1.1_genomic.fna
          ncbi_GCF_905332965.1_iyBomPasc1.1_genomic.gff

And a unified summary:
    genomes_out/download_summary_all.tsv

Usage examples (from Download_genomes/):

  python3 1.download_genomes_all.py \
      --ensembl_csv Genomes_selection_ensembl.csv \
      --out genomes_out

  python3 1.download_genomes_all.py \
      --ncbi_csv Genomes_selection_ncbi.csv \
      --out genomes_out

  python3 1.download_genomes_all.py \
      --both_csv Genomes_selection_both.csv \
      --out genomes_out

  # Or everything in one go:
  python3 1.download_genomes_all.py \
      --ensembl_csv Genomes_selection_ensembl.csv \
      --ncbi_csv Genomes_selection_ncbi.csv \
      --both_csv Genomes_selection_both.csv \
      --out genomes_out
"""

import argparse
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import unicodedata
import zipfile
from io import StringIO
from typing import List, Optional, Tuple

import requests
from bs4 import BeautifulSoup

# ---------------------------------------------------------------------------
# Global config / helpers
# ---------------------------------------------------------------------------

UA = {"User-Agent": "genome-downloader/unified/1.0"}
TIMEOUT = 45
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DEFAULT = os.path.join(SCRIPT_DIR, "genomes_out")
ACC_RE = re.compile(r"^(GC[AF])_(\d{9})\.(\d+)$", re.I)


def norm(s: Optional[str]) -> str:
    """Normalise whitespace and weird unicode."""
    if s is None:
        return ""
    s = unicodedata.normalize("NFKC", str(s))
    s = s.replace("\u00A0", " ").replace("\u202F", " ").replace("\ufeff", "")
    return re.sub(r"\s+", " ", s).strip()


def sanitize_species(name: str) -> str:
    """Turn free-text species name into safe folder name (underscores)."""
    s = re.sub(r"[^\w\s-]", "", norm(name))
    s = s.replace(" ", "_").replace("-", "_")
    return re.sub(r"_+", "_", s)


def detect_delim(first_line: str) -> str:
    """Guess delimiter from first line of CSV/TSV."""
    if not first_line:
        return ","
    scores = {
        ",": first_line.count(","),
        "\t": first_line.count("\t"),
        ";": first_line.count(";"),
    }
    return max(scores, key=scores.get)


def read_table(path: str) -> Tuple[csv.DictReader, dict]:
    """Load CSV/TSV into DictReader with automatic delimiter detection.

    Returns (reader, header_map) where header_map maps lowercase header -> original header.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    text = raw.decode("utf-8-sig", errors="ignore")
    first = text.splitlines()[0] if text else ""
    delim = detect_delim(first)
    rdr = csv.DictReader(StringIO(text), delimiter=delim)
    hdr = {(h or "").strip().lower(): h for h in (rdr.fieldnames or [])}
    return rdr, hdr


def choose_column(hdr: dict, *candidates: str) -> Optional[str]:
    """Find a column in hdr that matches any of the candidate names."""
    # exact
    for c in candidates:
        if c.lower() in hdr:
            return hdr[c.lower()]
    # substring / fuzzy
    for k, orig in hdr.items():
        for c in candidates:
            if c.lower() in k:
                return orig
    return None


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


# ---------------------------------------------------------------------------
# Ensembl helpers
# ---------------------------------------------------------------------------

def ens_slug(species: str) -> str:
    parts = norm(species).split()
    if not parts:
        return ""
    parts[0] = parts[0].capitalize()
    for i in range(1, len(parts)):
        parts[i] = parts[i].lower()
    return "_".join(parts)


def hget(url: str) -> requests.Response:
    r = requests.get(url, headers=UA, timeout=TIMEOUT)
    r.raise_for_status()
    return r


def list_links(url: str) -> List[str]:
    try:
        r = hget(url)
    except Exception:
        return []
    soup = BeautifulSoup(r.text, "html.parser")
    out: List[str] = []
    for a in soup.find_all("a"):
        href = a.get("href")
        if not href or href in ("../",):
            continue
        href = href.split("?")[0].split("#")[0]
        out.append(href)
    return out


def choose_first(items: List[str], patterns: Tuple[re.Pattern, ...]) -> Optional[str]:
    for pat in patterns:
        for it in items:
            if pat.search(it):
                return it
    return None


def crawl(base: str, patterns: Tuple[re.Pattern, ...], max_depth: int = 2) -> Optional[str]:
    base = base if base.endswith("/") else base + "/"
    seen = set()
    q = [(base, 0)]
    while q:
        url, d = q.pop(0)
        if url in seen or d > max_depth:
            continue
        seen.add(url)
        hrefs = list_links(url)
        files = [h for h in hrefs if not h.endswith("/")]
        m = choose_first(files, patterns)
        if m:
            return url + m
        for sd in [h for h in hrefs if h.endswith("/") and h != "../"]:
            q.append((url + sd, d + 1))
    return None


ENS_SOFTMASK = (
    re.compile(r".*soft[_-]?masked.*\.fa(\.gz)?$", re.I),
    re.compile(r".*dna_sm.*toplevel\.fa\.gz$", re.I),
)
ENS_FASTA_FALLBACKS = (
    re.compile(r".*dna\.toplevel\.fa\.gz$", re.I),
    re.compile(r".*toplevel\.fa\.gz$", re.I),
    re.compile(r".*primary_assembly\.fa\.gz$", re.I),
    re.compile(r".*\.fa\.gz$", re.I),
)
GFF_PATS = (
    re.compile(r".*genes\.gff3\.gz$", re.I),
    re.compile(r".*\.gff3?\.gz$", re.I),
    re.compile(r".*\.gtf\.gz$", re.I),
    re.compile(r".*\.gff3?$", re.I),
    re.compile(r".*\.gtf$", re.I),
)


def ensembl_asm_base(species: str) -> Optional[str]:
    root = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ens_slug(species)}/"
    dirs = [h for h in list_links(root) if h.endswith("/") and re.match(r"GC[AF]_\d+\.\d+/", h)]
    if not dirs:
        return None

    def key(d: str):
        m = re.match(r"GC[AF]_(\d+)\.(\d+)/", d)
        return (int(m.group(1)) if m else 0, int(m.group(2)) if m else 0)

    dirs.sort(key=key, reverse=True)
    return root + dirs[0]


def find_ensembl(species: str) -> Tuple[Optional[str], Optional[str]]:
    base = ensembl_asm_base(species)
    if not base:
        return None, None
    fasta = crawl(base + "genome/", ENS_SOFTMASK, 2) or crawl(base + "genome/", ENS_FASTA_FALLBACKS, 2)
    gff = None
    for root in ("ensembl/geneset/", "braker/geneset/", "geneset/"):
        gff = crawl(base + root, GFF_PATS, 2)
        if gff:
            break
    return fasta, gff


def save_unzipped(url: str, dest: str) -> None:
    ensure_dir(os.path.dirname(dest))
    with requests.get(url, stream=True, headers=UA, timeout=max(TIMEOUT, 120)) as r:
        r.raise_for_status()
        if url.endswith(".gz"):
            buf = io.BytesIO(r.content)
            with gzip.GzipFile(fileobj=buf) as gzf, open(dest, "wb") as out:
                shutil.copyfileobj(gzf, out)
        else:
            with open(dest, "wb") as out:
                shutil.copyfileobj(r.raw, out)


# ---------------------------------------------------------------------------
# NCBI helpers
# ---------------------------------------------------------------------------

def looks_like_fasta(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            chunk = fh.read(4096).lstrip()
        return chunk.startswith(b">")
    except Exception:
        return False


def choose_accession(gcf: str, gca: str) -> Optional[str]:
    gcf = norm(gcf)
    gca = norm(gca)
    m = ACC_RE.match(gcf)
    if m:
        return m.group(0).upper()
    m = ACC_RE.match(gca)
    if m:
        return m.group(0).upper()
    return None


def have_datasets_cli() -> bool:
    return shutil.which("datasets") is not None


def cli_download_zip(acc: str, dest_zip: str) -> None:
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        acc,
        "--include",
        "genome,gff3",
        "--filename",
        dest_zip,
        "--no-progressbar",
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def api_download_zip(acc: str, dest_zip: str) -> None:
    headers = {"Accept": "application/zip"}
    urls = [
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF&filename={acc}.zip",
        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF&filename={acc}.zip",
    ]
    last_err: Optional[Exception] = None
    for url in urls:
        try:
            with requests.get(url, headers=headers, stream=True, timeout=180) as r:
                r.raise_for_status()
                with open(dest_zip, "wb") as out:
                    for chunk in r.iter_content(1024 * 256):
                        if chunk:
                            out.write(chunk)
            return
        except Exception as e:
            last_err = e
    raise RuntimeError(f"REST download failed for {acc}: {last_err}")


def _extract_file_member(zf: zipfile.ZipFile, member_name: str, dest_path: str) -> None:
    """Extract a member from the zip. Gunzip if needed."""
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with zf.open(member_name) as src:
        if member_name.endswith(".gz"):
            with gzip.GzipFile(fileobj=src) as gzf, open(dest_path, "wb") as dst:
                shutil.copyfileobj(gzf, dst)
        else:
            with open(dest_path, "wb") as dst:
                shutil.copyfileobj(src, dst)


def extract_and_rename(zip_path: str, species_dir: str, acc: str, prefix: str = ""):
    """Extract FASTA+GFF from an NCBI datasets zip.

    If prefix is provided (e.g. 'ncbi_'), file names become:
      prefix + ACC_ASM_genomic.fna / .gff / .gtf
    """
    os.makedirs(species_dir, exist_ok=True)
    asm_name = None
    got_fna = got_gff = False
    fna_dest = gff_dest = ""

    with zipfile.ZipFile(zip_path) as z:
        members = z.namelist()

        fna_member = next((m for m in members if re.search(r"/[^/]*_genomic\.fna(\.gz)?$", m, re.I)), None)
        gff_member = next((m for m in members if re.search(r"/genomic\.(gff3?|gtf)(\.gz)?$", m, re.I)), None)

        if fna_member:
            m = re.search(r"/(GC[AF]_\d+\.\d+)_([^/]+)_genomic\.fna(\.gz)?$", fna_member, re.I)
            if m:
                acc_in_zip, asm_name = m.group(1), m.group(2)
            else:
                # fallback: dataset_catalog.json
                try:
                    cat_member = next(mm for mm in members if mm.endswith("dataset_catalog.json"))
                    cat = json.loads(z.read(cat_member).decode("utf-8"))
                    for df in cat.get("dataFiles", []):
                        for f in df.get("files", []):
                            fp = f.get("filePath", "")
                            m2 = re.search(r"/(GC[AF]_\d+\.\d+)_([^/]+)_genomic\.fna(\.gz)?$", fp, re.I)
                            if m2:
                                asm_name = m2.group(2)
                                break
                except Exception:
                    pass

            base_name = f"{acc}_{asm_name}_genomic.fna" if asm_name else f"{acc}_genomic.fna"
            out_name = f"{prefix}{base_name}"
            fna_dest = os.path.join(species_dir, out_name)
            _extract_file_member(z, fna_member, fna_dest)
            got_fna = looks_like_fasta(fna_dest)

        if gff_member:
            is_gtf = bool(re.search(r"\.gtf(\.gz)?$", gff_member, re.I))
            if asm_name:
                base_name = f"{acc}_{asm_name}_genomic.{'gtf' if is_gtf else 'gff'}"
            else:
                base_name = f"{acc}_genomic.{'gtf' if is_gtf else 'gff'}"
            out_name = f"{prefix}{base_name}"
            gff_dest = os.path.join(species_dir, out_name)
            _extract_file_member(z, gff_member, gff_dest)
            got_gff = True

    return got_fna, got_gff, fna_dest, gff_dest


def process_accession(acc: str, species: str, out_root: str, use_cli: bool, prefix: str = ""):
    """Download + extract one NCBI accession into out_root/species_dir.

    prefix is added in front of file names (used for 'both' species to get 'ncbi_*').
    """
    sp_dir = os.path.join(out_root, sanitize_species(species))
    os.makedirs(sp_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        zip_path = os.path.join(td, f"{acc}.zip")
        if use_cli:
            cli_download_zip(acc, zip_path)
        else:
            api_download_zip(acc, zip_path)

        ok_fa, ok_gff, fa_path, gff_path = extract_and_rename(zip_path, sp_dir, acc, prefix=prefix)

    msgs = []
    if ok_fa:
        msgs.append(f"FASTA -> {fa_path}")
    else:
        msgs.append("FASTA missing or invalid (no leading '>').")
    if ok_gff:
        msgs.append(f"GFF   -> {gff_path}")
    else:
        msgs.append("GFF not available for this assembly.")
    return "; ".join(msgs), (fa_path if ok_fa else ""), (gff_path if ok_gff else "")


# ---------------------------------------------------------------------------
# Download routines
# ---------------------------------------------------------------------------

def run_ensembl(csv_path: str, out_root: str, summary_writer: csv.writer) -> None:
    """Download genomes from Ensembl for all rows in csv_path into out_root."""
    rdr, hdr = read_table(csv_path)
    c_species = choose_column(hdr, "species", "organism", "taxon")
    if not c_species:
        raise SystemExit(f"[Ensembl] Missing 'Species' column in {csv_path}")

    for row in rdr:
        sp = norm(row.get(c_species, ""))
        if not sp:
            continue

        print(f"[ENSEMBL] {sp}")
        fa_url, gff_url = find_ensembl(sp)
        if not fa_url and not gff_url:
            print("  .. not found on Ensembl")
            summary_writer.writerow([sp, "ensembl", "", "", ""])
            continue

        sp_dir = os.path.join(out_root, sanitize_species(sp))
        ensure_dir(sp_dir)

        fa_out = ""
        gff_out = ""

        if fa_url:
            # Keep original basename in Ensembl-only case
            fa_out = os.path.join(sp_dir, os.path.basename(re.sub(r"\.gz$", "", fa_url)))
            print(f"  + FASTA -> {fa_out}")
            save_unzipped(fa_url, fa_out)
        if gff_url:
            gff_out = os.path.join(sp_dir, os.path.basename(re.sub(r"\.gz$", "", gff_url)))
            print(f"  + GFF   -> {gff_out}")
            save_unzipped(gff_url, gff_out)

        summary_writer.writerow([
            sp,
            "ensembl",
            "",  # accession not applicable
            os.path.abspath(fa_out) if fa_out else "",
            os.path.abspath(gff_out) if gff_out else "",
        ])


def run_ncbi(csv_path: str, out_root: str, summary_writer: csv.writer) -> None:
    """Download genomes from NCBI for all rows in csv_path into out_root."""
    rdr, hdr = read_table(csv_path)
    c_species = choose_column(hdr, "species", "organism", "taxon")
    c_gcf = choose_column(hdr, "gcf", "refseq", "ncbi_gcf")
    c_gca = choose_column(hdr, "gca", "genbank", "ncbi_gca")

    if not c_species:
        raise SystemExit(f"[NCBI] Missing 'Species' column in {csv_path}")
    if not (c_gcf or c_gca):
        raise SystemExit(f"[NCBI] Need at least one of 'GCF' or 'GCA' columns in {csv_path}")

    use_cli = have_datasets_cli()
    if use_cli:
        print("[NCBI] Using 'datasets' CLI.")
    else:
        print("[NCBI] 'datasets' CLI not found; using REST API.")

    for row in rdr:
        species = norm(row.get(c_species, ""))
        if not species:
            continue
        acc = choose_accession(row.get(c_gcf, ""), row.get(c_gca, ""))
        if not acc:
            print(f"[skip] {species}: no valid GCF/GCA")
            summary_writer.writerow([species, "ncbi", "", "", ""])
            continue

        try:
            # No prefix => classic NCBI names in the 'ncbi' subtree
            msg, fa_path, gff_path = process_accession(acc, species, out_root, use_cli, prefix="")
            print(f"[OK] {species} ({acc}): {msg}")
            summary_writer.writerow([
                species,
                "ncbi",
                acc,
                os.path.abspath(fa_path) if fa_path else "",
                os.path.abspath(gff_path) if gff_path else "",
            ])
        except subprocess.CalledProcessError as e:
            print(f"[fail] {species}: {acc} failed ({e})")
            summary_writer.writerow([species, "ncbi", acc, "", ""])
        except Exception as e:
            print(f"[fail] {species} ({acc}): {e}")
            summary_writer.writerow([species, "ncbi", acc, "", ""])


def run_both(csv_path: str, out_root: str, summary_writer: csv.writer) -> None:
    """For each row in csv_path, download from BOTH Ensembl and NCBI into out_root."""
    rdr, hdr = read_table(csv_path)
    c_species = choose_column(hdr, "species", "organism", "taxon")
    c_gcf = choose_column(hdr, "gcf", "refseq", "ncbi_gcf")
    c_gca = choose_column(hdr, "gca", "genbank", "ncbi_gca")

    if not c_species:
        raise SystemExit(f"[BOTH] Missing 'Species' column in {csv_path}")
    if not (c_gcf or c_gca):
        raise SystemExit(f"[BOTH] Need at least one of 'GCF' or 'GCA' columns in {csv_path}")

    use_cli = have_datasets_cli()
    if use_cli:
        print("[NCBI] Using 'datasets' CLI.")
    else:
        print("[NCBI] 'datasets' CLI not found; using REST API.")

    for row in rdr:
        species = norm(row.get(c_species, ""))
        if not species:
            continue

        sp_dir = os.path.join(out_root, sanitize_species(species))
        ensure_dir(sp_dir)

        # ---- Ensembl part (into out_root/both/<species>) ----
        print(f"[BOTH/ENSEMBL] {species}")
        fa_url, gff_url = find_ensembl(species)
        fa_out_ens = ""
        gff_out_ens = ""
        if not fa_url and not gff_url:
            print("  .. not found on Ensembl")
        else:
            if fa_url:
                # Decide whether we call it softmasked or generic
                if re.search(r"(soft[_-]?masked|dna_sm)", fa_url, re.I):
                    fa_name = "ensembl_softmasked.fa"
                else:
                    fa_name = "ensembl_genome.fa"
                fa_out_ens = os.path.join(sp_dir, fa_name)
                print(f"  + FASTA (Ensembl) -> {fa_out_ens}")
                save_unzipped(fa_url, fa_out_ens)
            if gff_url:
                # Decide extension
                if re.search(r"\.gtf(\.gz)?$", gff_url, re.I):
                    ext = "gtf"
                else:
                    ext = "gff3"
                gff_out_ens = os.path.join(sp_dir, f"ensembl_annotation.{ext}")
                print(f"  + GFF   (Ensembl) -> {gff_out_ens}")
                save_unzipped(gff_url, gff_out_ens)

        summary_writer.writerow([
            species,
            "ensembl",
            "",
            os.path.abspath(fa_out_ens) if fa_out_ens else "",
            os.path.abspath(gff_out_ens) if gff_out_ens else "",
        ])

        # ---- NCBI part (same <species> dir under out_root/both) ----
        acc = choose_accession(row.get(c_gcf, ""), row.get(c_gca, ""))
        if not acc:
            print(f"[BOTH/skip] {species}: no valid GCF/GCA")
            summary_writer.writerow([species, "ncbi", "", "", ""])
            continue

        try:
            # prefix="ncbi_" so names become ncbi_GCF_..._genomic.fna / .gff
            msg, fa_path, gff_path = process_accession(acc, species, out_root, use_cli, prefix="ncbi_")
            print(f"[BOTH/OK] {species} ({acc}): {msg}")
            summary_writer.writerow([
                species,
                "ncbi",
                acc,
                os.path.abspath(fa_path) if fa_path else "",
                os.path.abspath(gff_path) if gff_path else "",
            ])
        except Exception as e:
            print(f"[BOTH/fail] {species} ({acc}): {e}")
            summary_writer.writerow([species, "ncbi", acc, "", ""])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Unified downloader for Ensembl + NCBI genomes.")
    ap.add_argument("--ensembl_csv", help="CSV/TSV with species to fetch from Ensembl.")
    ap.add_argument("--ncbi_csv", help="CSV/TSV with species to fetch from NCBI.")
    ap.add_argument("--both_csv", help="CSV/TSV with species to fetch from BOTH Ensembl + NCBI.")
    ap.add_argument("--out", default=OUT_DEFAULT, help=f"Output root (default: {OUT_DEFAULT})")
    args = ap.parse_args()

    if not (args.ensembl_csv or args.ncbi_csv or args.both_csv):
        ap.error("Please provide at least one of --ensembl_csv, --ncbi_csv, or --both_csv")

    # Top-level out root (e.g. genomes_out)
    ensure_dir(args.out)
    # Sub-roots
    ens_root = os.path.join(args.out, "ensembl")
    ncbi_root = os.path.join(args.out, "ncbi")
    both_root = os.path.join(args.out, "both")
    ensure_dir(ens_root)
    ensure_dir(ncbi_root)
    ensure_dir(both_root)

    summary_path = os.path.join(args.out, "download_summary_all.tsv")

    with open(summary_path, "w", newline="") as sf:
        sw = csv.writer(sf, delimiter="\t")
        sw.writerow(["Species", "Source", "Accession", "FASTA", "GFF"])

        if args.ensembl_csv:
            print(f"=== ENSEMBL from {args.ensembl_csv} ===")
            run_ensembl(args.ensembl_csv, ens_root, sw)

        if args.ncbi_csv:
            print(f"\n=== NCBI from {args.ncbi_csv} ===")
            run_ncbi(args.ncbi_csv, ncbi_root, sw)

        if args.both_csv:
            print(f"\n=== BOTH from {args.both_csv} ===")
            run_both(args.both_csv, both_root, sw)

    print(f"\n[summary] wrote {summary_path}")


if __name__ == "__main__":
    main()
