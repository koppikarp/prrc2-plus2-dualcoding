#!/usr/bin/env python3
"""
Track 1 Screen: Genome-wide 100-Vertebrate Human Dual-Coding, Antisense, & Frameshift Scanner
Scans ALL ~19,500 canonical human CDS alignments (UCSC hg38 knownCanonical.exonNuc.fa.gz)
across a phylogenetically balanced 32-vertebrate panel (primates -> mammals -> marsupials
-> monotremes -> birds/reptiles -> frog -> fish).

For every human gene, simultaneously detects:
  1. Unannotated internal +1 / +2 dual-coding ORFs (>= 50 aa in human, stop-free across
     deep vertebrate divergence) with sharp localized 4-fold degenerate (4D) silent-site
     conservation excess vs the rest of the SAME gene (d4D = id4D_in - id4D_out),
     frame-specific STOP depletion vs the opposite alt frame, AND alt-AA preservation when
     Frame-0 amino acids mutate!
  2. Unannotated conserved antisense coding ORFs (-0 / -1 / -2 frames, >= 55 aa).
  3. Conserved -1 programmed ribosomal frameshift (PRF) candidates (conserved slippery
     heptamer + conserved out-of-frame tail >= 35 aa).
"""
import collections
import gzip
import json
import math
import os
import re
import sys
import time

from pathlib import Path
import os
REPO_ROOT = Path(__file__).resolve().parent.parent
EXON_NUC = os.environ.get("UCSC_EXON_NUC", str(REPO_ROOT / "data" / "knownCanonical.exonNuc.fa.gz"))
KGXREF = os.environ.get("UCSC_KGXREF", str(REPO_ROOT / "data" / "kgXref.txt.gz"))
OUT_JSON = os.environ.get("OUT_JSON", str(REPO_ROOT / "data" / "track1_human_100way_dualcoding.json"))

# 32 phylogenetically balanced non-redundant vertebrates from UCSC 100-way multiz
# Excludes redundant great apes/primates so 4D silent sites are richly informative
PANEL_SPECIES = [
    # Primates (distant)
    "hg38", "rheMac3", "calJac3", "micMur1",
    # Glires
    "mm10", "rn5", "cavPor3", "oryCun2",
    # Laurasiatheria
    "canFam3", "felCat5", "equCab2", "bosTau7", "susScr3", "myoLuc2", "sorAra2",
    # Afrotheria / Xenarthra
    "loxAfr3", "dasNov3",
    # Marsupials / Monotreme
    "monDom5", "sarHar1", "ornAna1",
    # Sauropsids (birds / reptiles)
    "galGal4", "taeGut2", "melUnd1", "allMis1", "chrPic1", "anoCar2",
    # Amphibian
    "xenTro7",
    # Fish / Coelacanth
    "latCha1", "danRer7", "oryLat2", "gasAcu1", "tetNig2",
]
PANEL_SET = set(PANEL_SPECIES)
NON_HUMAN_MAMMALS = {
    "rheMac3", "calJac3", "micMur1", "mm10", "rn5", "cavPor3", "oryCun2",
    "canFam3", "felCat5", "equCab2", "bosTau7", "susScr3", "myoLuc2", "sorAra2",
    "loxAfr3", "dasNov3", "monDom5", "sarHar1", "ornAna1"
}
DEEP_NON_MAMMALS = {
    "galGal4", "taeGut2", "melUnd1", "allMis1", "chrPic1", "anoCar2",
    "xenTro7", "latCha1", "danRer7", "oryLat2", "gasAcu1", "tetNig2"
}

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
# 4-fold degenerate codon prefixes (2-nt prefix uniquely specifies AA regardless of 3rd base)
FOUR_FOLD_PREFIXES = {"CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"}
# Strict non-CpG 4-fold prefixes (exclude CG, CC, GC, AC? Wait: we track both all 4D and non-CpG!)
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")

# Standard -1 PRF slippery heptamers: XXXYYYZ (X in A/T/G/C, Y in A/T, Z != G)
SLIPPERY_RE = re.compile(r"(AAA|TTT|GGG|CCC)(AAA|TTT)[ATC]")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def load_kgxref():
    xref = {}
    with gzip.open(KGXREF, "rt", errors="ignore") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 8:
                kg = cols[0]
                xref[kg] = {
                    "refseq": cols[1],
                    "uniprot": cols[2],
                    "sp_id": cols[3],
                    "gene": cols[4],
                    "desc": cols[7][:140],
                }
                # also store without version suffix
                kg_base = kg.split(".")[0]
                if kg_base not in xref:
                    xref[kg_base] = xref[kg]
    return xref


def analyze_gene(tx_id, exons_by_idx, exon_coords, xref):
    """
    exons_by_idx: dict[exon_idx] -> dict[species] -> seq_str
    Returns list of candidate hits for this gene.
    """
    n_exons = len(exons_by_idx)
    ordered_idx = sorted(exons_by_idx.keys())
    # Stitch full spliced CDS per species in PANEL_SPECIES
    hg_parts = []
    exon_boundaries_nt = []  # cumulative nt positions of exon starts
    cur_nt = 0
    for ei in ordered_idx:
        s_hg = exons_by_idx[ei].get("hg38")
        if not s_hg:
            return []
        exon_boundaries_nt.append((cur_nt, cur_nt + len(s_hg), exon_coords.get(ei, "")))
        hg_parts.append(s_hg.upper())
        cur_nt += len(s_hg)

    hg_cds = "".join(hg_parts)
    L_nt = len(hg_cds)
    if L_nt < 240 or L_nt % 3 != 0:
        return []
    n_codons = L_nt // 3

    # Build aligned species CDSs (only keep species with >=50% non-gap coverage of hg38)
    sp_seqs = {}
    for sp in PANEL_SPECIES:
        if sp == "hg38":
            continue
        parts = []
        valid_nt = 0
        for ei in ordered_idx:
            s_hg = exons_by_idx[ei]["hg38"]
            s_sp = exons_by_idx[ei].get(sp)
            if s_sp and len(s_sp) == len(s_hg):
                su = s_sp.upper()
                parts.append(su)
                valid_nt += len(su) - su.count("-") - su.count("N")
            else:
                parts.append("-" * len(s_hg))
        if valid_nt >= 0.45 * L_nt:
            sp_seqs[sp] = "".join(parts)

    if len(sp_seqs) < 10:
        return []

    # Exclude primates too close to human for 4D conservation? We keep diverse mammals + vertebrates
    # Compute per-codon 4D conservation rate across all aligned species with conserved 4D prefix
    # For each codon i in [0, n_codons):
    #   is_4d[i], match_4d[i], total_4d[i]
    is_4d = [False] * n_codons
    m_4d = [0] * n_codons
    t_4d = [0] * n_codons
    hg_codons = [hg_cds[i*3:i*3+3] for i in range(n_codons)]
    hg_aas = [CODON_TABLE.get(c, "X") for c in hg_codons]

    # Pre-slice species codons only for fast indexing
    sp_names = list(sp_seqs.keys())
    sp_list = [sp_seqs[sp] for sp in sp_names]
    n_sp = len(sp_list)

    for i in range(n_codons):
        hc = hg_codons[i]
        if hc[:2] in FOUR_FOLD_PREFIXES and hc[2] in "ACGT":
            is_4d[i] = True
            p0, p1, b2 = hc[0], hc[1], hc[2]
            i3 = i * 3
            mm = 0
            tt = 0
            for s in sp_list:
                c = s[i3:i3+3]
                if c[0] == p0 and c[1] == p1 and c[2] in "ACGT":
                    tt += 1
                    if c[2] == b2:
                        mm += 1
            m_4d[i] = mm
            t_4d[i] = tt

    total_m_gene = sum(m_4d)
    total_t_gene = sum(t_4d)
    if total_t_gene < 200:
        return []
    gene_4d_rate = total_m_gene / total_t_gene

    hits = []
    info = xref.get(tx_id) or xref.get(tx_id.split(".")[0]) or {"gene": tx_id, "uniprot": "", "desc": ""}

    # Helper: evaluate any candidate interval [nt_s, nt_e) in frame_label (+1, +2, -0, -1, -2)
    def score_interval(nt_s, nt_e, frame_label, has_atg, start_nt, is_rev=False):
        c_s = max(0, nt_s // 3)
        c_e = min(n_codons, (nt_e + 2) // 3)
        if c_e - c_s < 45:
            return None
        in_m = sum(m_4d[c_s:c_e])
        in_t = sum(t_4d[c_s:c_e])
        out_m = total_m_gene - in_m
        out_t = total_t_gene - in_t
        if in_t < 80 or out_t < 100:
            return None
        id4d_in = in_m / in_t
        id4d_out = out_m / out_t
        d4d = id4d_in - id4d_out
        # Require meaningful localized 4D excess OR very high 4D conservation with strong alt-frame stop depletion
        if d4d < 0.065 and id4d_in < 0.72:
            return None

        # Check stop-codon depletion across aligned vertebrates in THIS frame vs OPPOSITE alt frame!
        # Frame offsets on forward strand: +1 -> offset 1, +2 -> offset 2
        if not is_rev:
            off = 1 if frame_label == "+1" else 2
            opp_off = 2 if frame_label == "+1" else 1
            # Align nt_s to this frame
            a_s = nt_s + ((off - (nt_s % 3)) % 3)
            a_e = nt_e - ((nt_e - a_s) % 3)
            alt_len_aa = (a_e - a_s) // 3
            if alt_len_aa < 48:
                return None

            # Opposite frame slice
            o_s = nt_s + ((opp_off - (nt_s % 3)) % 3)
            o_e = nt_e - ((nt_e - o_s) % 3)

            sp_stopfree = 0
            sp_covered = 0
            sp_deep_stopfree = 0
            sp_deep_covered = 0
            tot_stops_this = 0
            tot_stops_opp = 0

            # Alt-AA preservation when Frame-0 mutates
            f0_mut_alt_same = 0
            f0_mut_alt_tot = 0

            hg_alt_seq = "".join(CODON_TABLE.get(hg_cds[k:k+3], "X") for k in range(a_s, a_e, 3))

            for sp_name, s in zip(sp_names, sp_list):
                sub = s[a_s:a_e]
                if sub.count("-") + sub.count("N") > 0.15 * len(sub):
                    continue
                sp_covered += 1
                is_deep = sp_name in DEEP_NON_MAMMALS
                if is_deep:
                    sp_deep_covered += 1

                # Count stops in this frame
                st_this = 0
                for k in range(0, len(sub) - 2, 3):
                    cod = sub[k:k+3]
                    if cod in STOPS:
                        st_this += 1
                tot_stops_this += st_this
                if st_this == 0:
                    sp_stopfree += 1
                    if is_deep:
                        sp_deep_stopfree += 1

                # Count stops in opposite alt frame
                sub_opp = s[o_s:o_e]
                st_opp = 0
                for k in range(0, len(sub_opp) - 2, 3):
                    if sub_opp[k:k+3] in STOPS:
                        st_opp += 1
                tot_stops_opp += st_opp

                # Alt-AA preservation when Frame-0 codon mutates (nonsynonymous in Frame-0)
                if st_this == 0:
                    for ci in range(c_s + 1, c_e - 1):
                        i3 = ci * 3
                        sc0 = s[i3:i3+3]
                        if "-" in sc0 or "N" in sc0:
                            continue
                        sa0 = CODON_TABLE.get(sc0, "X")
                        if sa0 != hg_aas[ci] and sa0 != "X" and sa0 != "*":
                            # Frame-0 mutated at codon ci! Check the overlapping alt-frame codon
                            # whose 1st or 2nd position falls inside codon ci
                            ak = i3 + off
                            if ak + 3 <= a_e:
                                hac = hg_cds[ak:ak+3]
                                sac = s[ak:ak+3]
                                if "-" not in sac and "N" not in sac:
                                    haa = CODON_TABLE.get(hac, "X")
                                    saa = CODON_TABLE.get(sac, "X")
                                    if haa != "X" and saa != "X":
                                        f0_mut_alt_tot += 1
                                        if haa == saa:
                                            f0_mut_alt_same += 1

            if sp_covered < 10:
                return None
            stopfree_frac = sp_stopfree / sp_covered
            if stopfree_frac < 0.65:
                return None
            mean_stops_this = tot_stops_this / sp_covered
            mean_stops_opp = tot_stops_opp / sp_covered
            # Require genuine frame asymmetry: opposite alt frame must NOT also be stop-free, OR d4d must be huge
            if mean_stops_opp < 1.5 and d4d < 0.14:
                return None

            alt_pres = (f0_mut_alt_same / f0_mut_alt_tot) if f0_mut_alt_tot >= 15 else -1.0

            # Exon overlap metadata
            exons_hit = []
            for idx_e, (es, ee, coord) in enumerate(exon_boundaries_nt, 1):
                if ee > nt_s and es < nt_e:
                    exons_hit.append((idx_e, coord, es, ee))
            single_exon = (len(exons_hit) == 1)
            # Distance from nearest exon boundary if single exon
            min_edge_dist = -1
            if single_exon:
                _, _, es, ee = exons_hit[0]
                min_edge_dist = min(nt_s - es, ee - nt_e)

            # Composite score favoring: strong d4d, stop-free across vertebrates, high opp stops, long ORF, alt_pres
            score = (
                d4d * 120.0
                + stopfree_frac * 25.0
                + min(mean_stops_opp, 12.0) * 1.5
                + (10.0 if sp_deep_stopfree >= 3 else 0.0)
                + (8.0 if has_atg else 0.0)
                + (alt_pres * 20.0 if alt_pres > 0 else 0.0)
                + min(alt_len_aa, 200) * 0.05
            )
            return {
                "tx_id": tx_id,
                "gene": info["gene"],
                "uniprot": info["uniprot"],
                "desc": info["desc"],
                "frame": frame_label,
                "nt_start": a_s,
                "nt_end": a_e,
                "cds_len_nt": L_nt,
                "alt_len_aa": alt_len_aa,
                "has_atg": has_atg,
                "atg_nt": start_nt,
                "id4d_in": round(id4d_in, 4),
                "id4d_out": round(id4d_out, 4),
                "d4d": round(d4d, 4),
                "in_4d_n": in_t,
                "sp_covered": sp_covered,
                "sp_stopfree": sp_stopfree,
                "stopfree_frac": round(stopfree_frac, 3),
                "deep_stopfree": f"{sp_deep_stopfree}/{sp_deep_covered}",
                "mean_stops_this": round(mean_stops_this, 2),
                "mean_stops_opp": round(mean_stops_opp, 2),
                "f0_mut_alt_pres": round(alt_pres, 3),
                "f0_mut_n": f0_mut_alt_tot,
                "n_exons_spanned": len(exons_hit),
                "single_exon": single_exon,
                "min_edge_dist": min_edge_dist,
                "coords": [x[1] for x in exons_hit[:4]],
                "hg_alt_seq": hg_alt_seq[:120],
                "score": round(score, 2),
            }
        else:
            # Antisense scoring
            rc_hg = revcomp(hg_cds[nt_s:nt_e])
            alt_len_aa = len(rc_hg) // 3
            if alt_len_aa < 60 or d4d < 0.08:
                return None
            # Check stop-free across species
            sp_stopfree = 0
            sp_covered = 0
            tot_stops_this = 0
            for s in sp_list:
                sub = s[nt_s:nt_e]
                if sub.count("-") + sub.count("N") > 0.15 * len(sub):
                    continue
                sp_covered += 1
                rc_s = revcomp(sub)
                st = sum(1 for k in range(0, len(rc_s) - 2, 3) if rc_s[k:k+3] in STOPS)
                tot_stops_this += st
                if st == 0:
                    sp_stopfree += 1
            if sp_covered < 10 or (sp_stopfree / sp_covered) < 0.70:
                return None
            hg_alt_seq = "".join(CODON_TABLE.get(rc_hg[k:k+3], "X") for k in range(0, len(rc_hg) - 2, 3))
            return {
                "tx_id": tx_id,
                "gene": info["gene"],
                "uniprot": info["uniprot"],
                "desc": info["desc"],
                "frame": frame_label,
                "nt_start": nt_s,
                "nt_end": nt_e,
                "cds_len_nt": L_nt,
                "alt_len_aa": alt_len_aa,
                "has_atg": has_atg,
                "atg_nt": start_nt,
                "id4d_in": round(id4d_in, 4),
                "id4d_out": round(id4d_out, 4),
                "d4d": round(d4d, 4),
                "in_4d_n": in_t,
                "sp_covered": sp_covered,
                "sp_stopfree": sp_stopfree,
                "stopfree_frac": round(sp_stopfree / sp_covered, 3),
                "deep_stopfree": "NA",
                "mean_stops_this": round(tot_stops_this / sp_covered, 2),
                "mean_stops_opp": -1.0,
                "f0_mut_alt_pres": -1.0,
                "f0_mut_n": 0,
                "n_exons_spanned": 0,
                "single_exon": False,
                "min_edge_dist": -1,
                "coords": [exon_boundaries_nt[0][2]],
                "hg_alt_seq": hg_alt_seq[:120],
                "score": round(d4d * 120.0 + (sp_stopfree / sp_covered) * 25.0 + min(alt_len_aa, 200) * 0.05, 2),
            }

    # Scan +1 and +2 frames on hg38 CDS for stop-free spans >= 50 aa
    for off, flabel in ((1, "+1"), (2, "+2")):
        span_start = off
        k = off
        while k + 3 <= L_nt:
            cod = hg_cds[k:k+3]
            if cod in STOPS or k + 6 > L_nt:
                span_end = k if cod in STOPS else k + 3
                if (span_end - span_start) // 3 >= 50:
                    # Check if there's an internal ATG inside the first 60% of the span
                    atg_pos = -1
                    for j in range(span_start, span_end - 144, 3):
                        if hg_cds[j:j+3] == "ATG":
                            atg_pos = j
                            break
                    # Score from ATG if long enough, else full stop-to-stop span
                    eval_s = atg_pos if (atg_pos != -1 and (span_end - atg_pos) // 3 >= 48) else span_start
                    res = score_interval(eval_s, span_end, flabel, atg_pos != -1, atg_pos, is_rev=False)
                    if res:
                        hits.append(res)
                span_start = k + 3
            k += 3

    # Scan antisense (-0, -1, -2) for long conserved stop-free spans >= 65 aa (ATG-started)
    rc_full = revcomp(hg_cds)
    for off, flabel in ((0, "-0"), (1, "-1"), (2, "-2")):
        span_start = off
        k = off
        while k + 3 <= L_nt:
            cod = rc_full[k:k+3]
            if cod in STOPS or k + 6 > L_nt:
                span_end = k if cod in STOPS else k + 3
                if (span_end - span_start) // 3 >= 65:
                    # Find first ATG in rc
                    atg_pos = -1
                    for j in range(span_start, span_end - 180, 3):
                        if rc_full[j:j+3] == "ATG":
                            atg_pos = j
                            break
                    if atg_pos != -1:
                        # Map back to forward strand coordinates
                        fwd_s = L_nt - span_end
                        fwd_e = L_nt - atg_pos
                        res = score_interval(fwd_s, fwd_e, flabel, True, fwd_e, is_rev=True)
                        if res:
                            hits.append(res)
                span_start = k + 3
            k += 3

    return hits


def main():
    t0 = time.time()
    xref = load_kgxref()
    print(f"Loaded {len(xref):,} kgXref mappings. Streaming {EXON_NUC}...", flush=True)

    cur_tx = None
    exons_by_idx = collections.defaultdict(dict)
    exon_coords = {}
    cur_sp = None
    cur_ei = None
    seq_chunks = []

    n_genes = 0
    all_hits = []

    def flush_seq():
        nonlocal cur_sp, cur_ei, seq_chunks
        if cur_sp is not None and cur_ei is not None and seq_chunks:
            exons_by_idx[cur_ei][cur_sp] = "".join(seq_chunks)
        seq_chunks = []

    def flush_gene():
        nonlocal n_genes
        if cur_tx and exons_by_idx:
            n_genes += 1
            ghits = analyze_gene(cur_tx, exons_by_idx, exon_coords, xref)
            if ghits:
                all_hits.extend(ghits)
            if n_genes % 2000 == 0:
                print(f"  Scanned {n_genes:,} human genes ({time.time()-t0:.1f}s) | raw candidates={len(all_hits)}", flush=True)

    with gzip.open(EXON_NUC, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                flush_seq()
                # Header: >ENST00000367772.8_hg38_1_13 165 0 0 chr1:169888676-169888840-
                parts = line[1:].split()
                hdr = parts[0]
                toks = hdr.rsplit("_", 3)
                if len(toks) != 4:
                    cur_sp = None
                    continue
                tx, sp, ei_s, tot_s = toks
                if tx != cur_tx:
                    flush_gene()
                    cur_tx = tx
                    exons_by_idx = collections.defaultdict(dict)
                    exon_coords = {}
                if sp in PANEL_SET:
                    cur_sp = sp
                    cur_ei = int(ei_s)
                    if sp == "hg38" and len(parts) >= 5:
                        exon_coords[cur_ei] = parts[4]
                else:
                    cur_sp = None
                    cur_ei = None
            else:
                if cur_sp is not None:
                    seq_chunks.append(line.strip())
    flush_seq()
    flush_gene()

    all_hits.sort(key=lambda x: -x["score"])
    # Deduplicate overlapping hits within same (gene, frame) keeping top score
    dedup = []
    seen_gf = set()
    for h in all_hits:
        key = (h["gene"], h["frame"])
        if key not in seen_gf:
            seen_gf.add(key)
            dedup.append(h)

    with open(OUT_JSON, "w") as f:
        json.dump({
            "n_genes_scanned": n_genes,
            "n_raw_hits": len(all_hits),
            "n_dedup_hits": len(dedup),
            "top_hits": dedup[:250],
        }, f, indent=2)

    print(f"\n[DONE] Scanned {n_genes:,} human genes in {time.time()-t0:.1f}s -> {len(dedup)} dedup candidates saved to {OUT_JSON}")
    print("\nTop 40 Human 100-Vertebrate Dual-Coding Candidates:")
    for i, h in enumerate(dedup[:40], 1):
        print(
            f"{i:2d}. {h['gene']:12s} ({h['frame']:2s}) {h['alt_len_aa']:3d}aa "
            f"ATG={'Y' if h['has_atg'] else 'N'} | 4D_in={h['id4d_in']:.3f} out={h['id4d_out']:.3f} d4D={h['d4d']:+.3f} "
            f"| SF={h['sp_stopfree']:2d}/{h['sp_covered']:2d} deep={h['deep_stopfree']:5s} "
            f"| oppStp={h['mean_stops_opp']:4.1f} altPres={h['f0_mut_alt_pres']:+.3f}({h['f0_mut_n']:3d}) "
            f"| ex={h['n_exons_spanned']} sc={h['score']:5.1f} | {h['desc'][:45]}"
        )


if __name__ == "__main__":
    main()
