#!/usr/bin/env python3
"""
1. Fetch Drosophila melanogaster nocte, Octopus, Ciona, and Trichoplax PRRC2/BAT2 CDSs
   via UniProt -> EMBL/RefSeq to pinpoint the EXACT evolutionary origin of alt-PRRC2
   (Deuterostome birth ~600 My ago vs Bilaterian/Metazoan birth!).
2. Print the exact nucleotide codon-by-codon correspondence between Frame-0 (BAT2)
   and Frame +2 (alt-PRRC2) across Sea Urchin, Amphioxus, Lamprey, Elephant Shark,
   Human PRRC2A, Human PRRC2B, and Human PRRC2C!
"""
import json
import re
import time
import urllib.parse
import urllib.request

CT = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def fetch(url, pause=0.45):
    time.sleep(pause)
    req = urllib.request.Request(url, headers={"User-Agent": "AntigravityResearch/1.0"})
    with urllib.request.urlopen(req, timeout=25) as r:
        return r.read().decode("utf-8", errors="ignore")

def transl(dna, off=0):
    return "".join(CT.get(dna[i:i+3], "X") for i in range(off, len(dna) - 2, 3))

# Query UniProt for InterPro IPR009738 (BAT2 N-terminal domain!) in key invertebrate phyla:
# Drosophila (7227), Octopus/mollusc (6447), Nematode/C. elegans (6239), Cnidaria/Nematostella (45351), Trichoplax (10228), Ciona (7719), Acorn worm/Saccoglossus (10224)
taxa = [
    ("Drosophila_melanogaster", "7227"),
    ("Ciona_intestinalis", "7719"),
    ("Saccoglossus_kowalevskii", "10224"),
    ("Octopus_bimaculoides", "37653"),
    ("Nematostella_vectensis", "45351"),
    ("Trichoplax_adhaerens", "10228"),
]
for name, taxid in taxa:
    u = f"https://rest.uniprot.org/uniprotkb/search?query=xref:interpro-IPR009738+AND+organism_id:{taxid}&format=json&size=2"
    try:
        d = json.loads(fetch(u, 0.3))
        res = d.get("results", [])
        if not res:
            print(f"[{name}] No IPR009738 hit in UniProt")
            continue
        r = res[0]
        acc = r["primaryAccession"]
        # Find EMBL or RefSeq nucleotide xref
        embl_ids = [x["id"] for x in r.get("uniProtKBCrossReferences", []) if x["database"] == "RefSeq"]
        embl_cds = []
        for x in r.get("uniProtKBCrossReferences", []):
            if x["database"] == "EMBL":
                props = {p["key"]: p["value"] for p in x.get("properties", [])}
                if props.get("ProteinId") and props.get("ProteinId") != "-":
                    embl_cds.append((x["id"], props["ProteinId"], props.get("MoleculeType")))
        print(f"\n[{name}] UniProt={acc} len={r['sequence']['length']} | RefSeq={embl_ids[:2]} | EMBL={embl_cds[:2]}")
        # Fetch CDS via ENA REST if EMBL ProteinId exists, or NCBI RefSeq
        dna = None
        if embl_cds:
            pid = embl_cds[0][1]
            try:
                fa = fetch(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{pid}?download=false", 0.3)
                if fa.startswith(">"):
                    dna = "".join(ln.strip() for ln in fa.splitlines() if not ln.startswith(">")).upper()
            except Exception:
                pass
        if not dna and embl_ids:
            try:
                gp = fetch(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={embl_ids[0]}&rettype=gp&retmode=text", 0.45)
                m_cb = re.search(r'/coded_by="([A-Z0-9_]+\.\d+):(\d+)\.\.(\d+)"', gp)
                if m_cb:
                    racc, s1, e1 = m_cb.group(1), int(m_cb.group(2))-1, int(m_cb.group(3))
                    fa = fetch(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={racc}&rettype=fasta&retmode=text", 0.45)
                    seq = "".join(ln.strip() for ln in fa.splitlines() if not ln.startswith(">")).upper()
                    dna = seq[s1:e1]
            except Exception as e:
                print("  RefSeq fetch err:", e)
        if dna and len(dna) >= 600:
            p0, p1, p2 = transl(dna, 0), transl(dna, 1), transl(dna, 2)
            print(f"  CDS={len(dna)}nt ({len(p0)}aa)")
            print(f"  F0[0:95]  : {p0[:95]}")
            print(f"  F0[75:225]: {p0[75:225]}")
            print(f"  +2[75:225]: {p2[75:225]}")
            print(f"  +1[75:225]: {p1[75:225]}")
            print(f"  Stops in [80:215]: +1={p1[80:215].count('*')}, +2={p2[80:215].count('*')}")
        else:
            print("  Could not fetch CDS DNA; F0 head:", r["sequence"]["value"][:180])
    except Exception as e:
        print(f"[{name}] ERR: {e}")
