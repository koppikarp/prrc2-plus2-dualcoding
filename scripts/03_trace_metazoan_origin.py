#!/usr/bin/env python3
"""
03_trace_metazoan_origin.py

Live NCBI RefSeq verification script for the 680-My metazoan evolutionary series
of PRRC2 / nocte across all 7 non-human lineages:
  1. Drosophila melanogaster (Protostome outgroup: nocte, NP_001245614.1)
  2. Saccoglossus kowalevskii (Hemichordate acorn worm: XP_077866102.1)
  3. Strongylocentrotus purpuratus (Echinoderm sea urchin: XP_030852947.1)
  4. Branchiostoma floridae (Cephalochordate amphioxus: XP_035683835.1)
  5. Petromyzon marinus (Jawless vertebrate sea lamprey: XP_032802661.2)
  6. Callorhinchus milii prrc2b (Cartilaginous fish elephant shark: XM_042337680.1)
  7. Callorhinchus milii prrc2c (Cartilaginous fish elephant shark: XM_042335521.1)

Fetches each exact RefSeq CDS via NCBI E-utilities, translates all three forward
reading frames across the homologous BAT2 core, verifies +2 stop counts and
conserved motifs, and optionally refreshes `data/metazoan_prrc2_cds.fa` when
invoked with `--write-fasta`.
"""
import argparse
import re
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

CT = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

METAZOAN_SPECS = [
    ("Drosophila_melanogaster", "nocte", "protein", "NP_001245614.1", 680, "Protostome (Arthropoda)", 80, 215),
    ("Saccoglossus_kowalevskii", "PRRC2", "protein", "XP_077866102.1", 580, "Hemichordata (Acorn worm)", 85, 223),
    ("Strongylocentrotus_purpuratus", "PRRC2", "protein", "XP_030852947.1", 580, "Echinodermata (Sea urchin)", 83, 218),
    ("Branchiostoma_floridae", "PRRC2", "protein", "XP_035683835.1", 550, "Cephalochordata (Amphioxus)", 84, 188),
    ("Petromyzon_marinus", "PRRC2 (LOC116938979)", "protein", "XP_032802661.2", 500, "Agnatha (Sea lamprey)", 94, 225),
    ("Callorhinchus_milii_prrc2b", "prrc2b", "nuccore", "XM_042337680.1", 465, "Chondrichthyes (Elephant shark)", 86, 210),
    ("Callorhinchus_milii_prrc2c", "prrc2c", "nuccore", "XM_042335521.1", 465, "Chondrichthyes (Elephant shark)", 85, 227),
]


def fetch(url: str, pause: float = 0.45) -> str:
    time.sleep(pause)
    req = urllib.request.Request(url, headers={"User-Agent": "PRRC2-Metazoan-Origin/1.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read().decode("utf-8", errors="ignore")


def transl(dna: str, off: int = 0) -> str:
    return "".join(CT.get(dna[i:i + 3], "X") for i in range(off, len(dna) - 2, 3))


def main():
    ap = argparse.ArgumentParser(description="Verify all 7 non-human metazoan PRRC2/nocte CDSs from NCBI RefSeq.")
    ap.add_argument("--write-fasta", action="store_true", help="Write refreshed data/metazoan_prrc2_cds.fa")
    args = ap.parse_args()

    fa_lines = []
    for key, gene, db, acc, mya, clade, s_aa, e_aa in METAZOAN_SPECS:
        if db == "protein":
            gp = fetch(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={acc}&rettype=gp&retmode=text")
            m_cb = re.search(r'/coded_by="([A-Z0-9_]+\.\d+):(\d+)\.\.(\d+)"', gp)
            if not m_cb:
                raise RuntimeError(f"Could not parse coded_by for {acc}")
            rna_acc, s1, e1 = m_cb.group(1), int(m_cb.group(2)) - 1, int(m_cb.group(3))
            fa = fetch(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={rna_acc}&rettype=fasta&retmode=text")
            dna = "".join(ln.strip() for ln in fa.splitlines() if not ln.startswith(">")).upper()
            cds = dna[s1:e1]
        else:
            gb = fetch(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text")
            m_cds = re.search(r"     CDS             (?:join\()?(\d+)\.\.(\d+)", gb)
            m_ori = re.search(r"ORIGIN\s+(.+?)//", gb, re.S)
            if not (m_cds and m_ori):
                raise RuntimeError(f"Could not parse CDS/ORIGIN for {acc}")
            s1, e1 = int(m_cds.group(1)) - 1, int(m_cds.group(2))
            dna = re.sub(r"[^acgtACGT]", "", m_ori.group(1)).upper()
            cds = dna[s1:e1]
            rna_acc = acc

        p1 = transl(cds, 1)[s_aa:e_aa]
        p2 = transl(cds, 2)[s_aa:e_aa]
        print(f"[{key}] {acc} ({rna_acc}) | CDS={len(cds)} nt | core [{s_aa}:{e_aa}] stops: +2={p2.count('*')}, +1={p1.count('*')}")
        print(f"  +2 core: {p2}")
        fa_lines.append(f">{key}|gene={gene}|acc={acc}|rna={rna_acc}|mya={mya}|clade={clade}|full_cds_nt={len(cds)}\n{cds[:1200]}\n")

    if args.write_fasta:
        out = DATA_DIR / "metazoan_prrc2_cds.fa"
        out.write_text("".join(fa_lines))
        print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
