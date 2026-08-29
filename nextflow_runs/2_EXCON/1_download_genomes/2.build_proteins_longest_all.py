#!/usr/bin/env python3
import argparse
import os
import re
import subprocess

def log(msg):
    print(msg, flush=True)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def list_subdirs(path):
    if not os.path.isdir(path):
        return []
    return sorted(
        d for d in os.listdir(path)
        if os.path.isdir(os.path.join(path, d)) and not d.startswith(".")
    )

# ---------- Pick genome + GFF for each species / source ----------

def choose_genome_and_gff(species_dir, source, root_label):
    """
    Heuristic to pick genome FASTA and GFF/GTF for a given species/source.

    source: 'ensembl' or 'ncbi'
    root_label: 'ensembl', 'ncbi', or 'both'
    """
    files = [
        f for f in os.listdir(species_dir)
        if os.path.isfile(os.path.join(species_dir, f))
    ]
    if not files:
        return None, None

    def candidates(exts, prefix_filter=None, exclude_prefix=None):
        out = []
        for f in files:
            fl = f.lower()
            if not any(fl.endswith(ext) for ext in exts):
                continue
            if prefix_filter is not None and not f.startswith(prefix_filter):
                continue
            if exclude_prefix is not None and f.startswith(exclude_prefix):
                continue
            out.append(f)
        return out

    genome = None
    gff = None

    if root_label == "both":
        if source == "ensembl":
            # Genome: prefer standard names, else any non-ncbi .fa/.fna
            for name in ("ensembl_softmasked.fa", "ensembl_genome.fa"):
                if name in files:
                    genome = name
                    break
            if genome is None:
                g_cands = candidates((".fa", ".fna", ".fasta"), exclude_prefix="ncbi_")
                if g_cands:
                    g_cands_sorted = sorted(
                        g_cands,
                        key=lambda x: (not re.search(r"(softmasked|dna_sm)", x.lower()), x),
                    )
                    genome = g_cands_sorted[0]

            # GFF: prefer ensembl_annotation.*, else any non-ncbi GFF/GTF
            for name in ("ensembl_annotation.gff3", "ensembl_annotation.gtf"):
                if name in files:
                    gff = name
                    break
            if gff is None:
                gff_cands = candidates((".gff3", ".gff", ".gtf"), exclude_prefix="ncbi_")
                if gff_cands:
                    gff = sorted(gff_cands)[0]

        elif source == "ncbi":
            # Genome: prefer ncbi_*.fna/fa
            g_cands = candidates((".fna", ".fa", ".fasta"), prefix_filter="ncbi_")
            if g_cands:
                genome = sorted(g_cands)[0]
            # GFF: prefer ncbi_*.gff/gtf
            gff_cands = candidates((".gff3", ".gff", ".gtf"), prefix_filter="ncbi_")
            if gff_cands:
                gff = sorted(gff_cands)[0]

    else:
        # root_label is 'ensembl' or 'ncbi' with only one source
        g_cands = candidates((".fna", ".fa", ".fasta"))
        if g_cands:
            g_cands_sorted = sorted(
                g_cands,
                key=lambda x: (not re.search(r"(softmasked|dna_sm)", x.lower()), x),
            )
            genome = g_cands_sorted[0]

        gff_cands = candidates((".gff3", ".gff", ".gtf"))
        if gff_cands:
            gff = sorted(gff_cands)[0]

    if genome:
        genome = os.path.join(species_dir, genome)
    if gff:
        gff = os.path.join(species_dir, gff)
    return genome, gff

# ---------- GFF strand fix ----------

def fix_gff_strand(gff_path):
    """Replace '?' strand in GFF column 7 with '.' (in-place, with .bak backup)."""
    if not os.path.exists(gff_path):
        return

    with open(gff_path, "r") as fh:
        txt = fh.read()

    if "\t?\t" not in txt:
        return  # nothing to fix

    bak = gff_path + ".bak"
    if not os.path.exists(bak):
        with open(bak, "w") as out:
            out.write(txt)

    out_lines = []
    for line in txt.splitlines(keepends=False):
        if not line or line.startswith("#"):
            out_lines.append(line)
            continue
        cols = line.split("\t")
        if len(cols) >= 8 and cols[6] == "?":
            cols[6] = "."
            line = "\t".join(cols)
        out_lines.append(line)

    with open(gff_path, "w") as out:
        out.write("\n".join(out_lines) + "\n")

# ---------- Longest-isoform helper functions ----------

def parse_attrs(field: str) -> dict:
    """Parse GFF3 or GTF attribute column into a dict."""
    s = field.strip().strip(";")
    attrs = {}
    if not s:
        return attrs

    # GTF-like?
    if 'transcript_id "' in s or 'gene_id "' in s:
        parts = s.split(";")
        for part in parts:
            part = part.strip()
            if not part:
                continue
            if " " not in part:
                continue
            key, val = part.split(" ", 1)
            val = val.strip().strip('"')
            attrs[key.strip()] = val
    else:
        # GFF3-like
        parts = s.split(";")
        for part in parts:
            part = part.strip()
            if not part or "=" not in part:
                continue
            key, val = part.split("=", 1)
            attrs[key.strip()] = val.strip()
    return attrs


def build_tx_to_gene(gff_path: str):
    tx_to_gene = {}
    with open(gff_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2].lower()
            if ftype not in ("mrna", "transcript", "ncrna", "lncrna", "lnc_rna", "rrna", "trna"):
                continue
            attrs = parse_attrs(cols[8])
            tid = attrs.get("ID") or attrs.get("transcript_id")
            if not tid:
                continue
            gid = (
                attrs.get("gene_id")
                or attrs.get("locus_tag")
                or attrs.get("gene")
                or attrs.get("Parent")
                or attrs.get("parent")
            )
            if not gid:
                continue
            tx_to_gene[tid] = gid
    return tx_to_gene


def read_fasta(path: str):
    seqs = {}
    cur_id = None
    chunks = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(chunks)
                cur_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if cur_id is not None:
            seqs[cur_id] = "".join(chunks)
    return seqs


def select_longest_isoforms(tx_to_gene, prot_seqs):
    best = {}
    for tx, gene in tx_to_gene.items():
        seq = prot_seqs.get(tx)
        if not seq:
            continue
        L = len(seq)
        if gene not in best or L > best[gene][1]:
            best[gene] = (tx, L)
    keep = {tx for tx, _ in best.values()}
    return keep


def write_filtered_fasta(in_fa, out_fa, keep_ids):
    with open(in_fa, "r") as inp, open(out_fa, "w") as out:
        cur_id = None
        cur_header = None
        buf = []

        def flush():
            nonlocal cur_id, cur_header, buf
            if cur_id is not None and cur_id in keep_ids:
                out.write(cur_header + "\n")
                out.write("".join(buf))

        for line in inp:
            if line.startswith(">"):
                flush()
                cur_header = line.rstrip("\n")
                cur_id = cur_header[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        flush()

# ---------- gffread wrapper + per-species processing ----------

def run_gffread(gff, genome, out_prot, gffread_bin="gffread"):
    cmd = [gffread_bin, gff, "-g", genome, "-y", out_prot]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        log(f"[ERROR] gffread failed on {gff} with genome {genome}: {res.stderr.strip()}")
        return False
    return True


def process_one(species_dir, proteins_out_dir, species_name, source, root_label, gffread_bin="gffread"):
    genome, gff = choose_genome_and_gff(species_dir, source, root_label)
    if not genome or not gff:
        log(f"[WARN] {root_label}/{species_name} ({source}): missing genome or GFF, skipping")
        return

    log(f"[INFO] {root_label}/{species_name} ({source}): genome={os.path.basename(genome)}, gff={os.path.basename(gff)}")

    fix_gff_strand(gff)

    # Decide base name for proteins
    if root_label == "both":
        if source == "ensembl":
            base = f"{species_name}_ENSEMBL_proteins"
        else:
            base = f"{species_name}_NCBI_proteins"
    elif root_label == "ensembl":
        base = f"{species_name}_ENSEMBL_proteins"
    else:  # root_label == 'ncbi'
        base = f"{species_name}_NCBI_proteins"

    out_prot = os.path.join(proteins_out_dir, base + ".fa")
    out_long = os.path.join(proteins_out_dir, base + "_longest.fa")

    ensure_dir(proteins_out_dir)

    log(f"[INFO]  - running gffread -> {out_prot}")
    ok = run_gffread(gff, genome, out_prot, gffread_bin=gffread_bin)
    if not ok:
        return

    log(f"[INFO]  - building transcript→gene map from {gff}")
    tx_to_gene = build_tx_to_gene(gff)
    log(f"[INFO]    transcripts with gene mapping: {len(tx_to_gene)}")

    log(f"[INFO]  - reading proteins from {out_prot}")
    prot_seqs = read_fasta(out_prot)
    log(f"[INFO]    proteins: {len(prot_seqs)}")

    keep_ids = select_longest_isoforms(tx_to_gene, prot_seqs)
    log(f"[INFO]  - keeping {len(keep_ids)} transcripts (longest isoform per gene)")

    log(f"[INFO]  - writing longest-isoform proteins to {out_long}")
    write_filtered_fasta(out_prot, out_long, keep_ids)

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Extract proteins with gffread and keep longest isoform per gene, "
                    "for all species in genomes_out/ensembl, genomes_out/ncbi and genomes_out/both."
    )
    ap.add_argument(
        "--base",
        default="genomes_out",
        help="Base genomes_out directory (default: genomes_out)",
    )
    ap.add_argument(
        "--gffread",
        default="gffread",
        help="Path to gffread binary (default: gffread)",
    )
    args = ap.parse_args()

    base = os.path.abspath(args.base)
    gffread_bin = args.gffread

    log(f"[INFO] Base directory: {base}")
    for root_label in ("ensembl", "ncbi", "both"):
        root_dir = os.path.join(base, root_label)
        if not os.path.isdir(root_dir):
            continue
        proteins_out_dir = os.path.join(root_dir, "00_proteins_out")
        ensure_dir(proteins_out_dir)

        log(f"[INFO] Processing root: {root_label} ({root_dir})")
        for species_name in list_subdirs(root_dir):
            if species_name == "00_proteins_out":
                continue
            species_dir = os.path.join(root_dir, species_name)
            if not os.path.isdir(species_dir):
                continue

            if root_label == "both":
                # For species with both, do Ensembl and NCBI
                process_one(species_dir, proteins_out_dir, species_name, "ensembl", root_label, gffread_bin=gffread_bin)
                process_one(species_dir, proteins_out_dir, species_name, "ncbi", root_label, gffread_bin=gffread_bin)
            elif root_label == "ensembl":
                process_one(species_dir, proteins_out_dir, species_name, "ensembl", root_label, gffread_bin=gffread_bin)
            else:  # ncbi
                process_one(species_dir, proteins_out_dir, species_name, "ncbi", root_label, gffread_bin=gffread_bin)

    log("[OK] Finished processing all species.")

if __name__ == "__main__":
    main()

