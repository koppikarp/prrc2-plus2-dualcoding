#!/usr/bin/env python3
"""
01_reproduce_all_analyses.py

Self-contained master analysis and verification pipeline for:
"An Uncharacterized 600-Million-Year-Old +2 Dual-Coding Gene Family Embedded
Inside the Human Leaky-Scanning Regulators PRRC2A and PRRC2B"

Recomputes every vertebrate comparative statistic, within-gene flanking-codon
synonymous null, 807,162-human gnomAD v4 stratified Fisher & exact
hypergeometric Cochran-Mantel-Haenszel (CMH) test, GENCODE Ribo-seq uoORF
architecture, ClinVar dual-frame re-annotation, and CBE knockout switch table
directly from the bundled raw data subsets in `data/`, exporting
`data/master_preprint_data.json` for all Fignimble figure panels.

Zero external files required. Runs in ~3 seconds with Python 3.10+, numpy, scipy.
"""

import collections
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy import optimize, stats

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
OUT_JSON = DATA_DIR / "master_preprint_data.json"

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
STOPS = {"TAA", "TAG", "TGA"}
ALL_4D = {"CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"}
TI = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
COMP = str.maketrans("ACGTN-", "TGCAN-")

# Human codon usage per 1,000 codons (Kazusa H. sapiens reference table)
HS_CODON = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0, "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4, "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1, "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8, "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
    "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9, "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
    "TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2, "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
    "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0, "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}

# Non-redundant 32-vertebrate reference panel spanning primates to teleost fish
PANEL_32 = [
    "rheMac3", "calJac3", "micMur1", "mm10", "rn5", "cavPor3", "oryCun2",
    "canFam3", "felCat5", "equCab2", "bosTau7", "susScr3", "myoLuc2", "sorAra2",
    "loxAfr3", "dasNov3", "monDom5", "sarHar1", "ornAna1",
    "galGal4", "taeGut2", "melUnd1", "allMis1", "chrPic1", "anoCar2",
    "xenTro7", "latCha1", "lepOcu1", "danRer7", "oryLat2", "gasAcu1", "tetNig2",
]

# Major vertebrate clades for clade-collapsed phylogenetic independence check
CLADE_MAP = {
    "rheMac3": "Primates", "calJac3": "Primates", "micMur1": "Primates",
    "mm10": "Glires", "rn5": "Glires", "cavPor3": "Glires", "oryCun2": "Glires",
    "canFam3": "Laurasiatheria", "felCat5": "Laurasiatheria", "equCab2": "Laurasiatheria",
    "bosTau7": "Laurasiatheria", "susScr3": "Laurasiatheria", "myoLuc2": "Laurasiatheria", "sorAra2": "Laurasiatheria",
    "loxAfr3": "Atlantogenata", "dasNov3": "Atlantogenata",
    "monDom5": "Non-eutherian mammals", "sarHar1": "Non-eutherian mammals", "ornAna1": "Non-eutherian mammals",
    "galGal4": "Sauropsida", "taeGut2": "Sauropsida", "melUnd1": "Sauropsida",
    "allMis1": "Sauropsida", "chrPic1": "Sauropsida", "anoCar2": "Sauropsida",
    "xenTro7": "Amphibia",
    "latCha1": "Fish/Sarcopterygii", "lepOcu1": "Fish/Sarcopterygii",
    "danRer7": "Fish/Sarcopterygii", "oryLat2": "Fish/Sarcopterygii",
    "gasAcu1": "Fish/Sarcopterygii", "tetNig2": "Fish/Sarcopterygii",
}

TARGETS = {
    "ENST00000376033.3": ("PRRC2A", 254, 644, "6"),
    "ENST00000357304.8": ("PRRC2B", 257, 683, "9"),
    "ENST00000338920.8": ("PRRC2C", 248, 650, "1"),
}


def transl(dna: str, off: int = 0) -> str:
    return "".join(CT.get(dna[i:i + 3], "X") for i in range(off, len(dna) - 2, 3))


def classify_change(cds: str, idx: int, alt_base: str, frame_off: int) -> str:
    if idx < frame_off:
        return "outside"
    c_start = frame_off + ((idx - frame_off) // 3) * 3
    if c_start + 3 > len(cds):
        return "outside"
    rc = cds[c_start:c_start + 3]
    ac = rc[:idx - c_start] + alt_base + rc[idx - c_start + 1:]
    ra, aa = CT.get(rc, "X"), CT.get(ac, "X")
    if ra == "*":
        return "stop_ref"
    if aa == "*":
        return "stop_gained"
    return "syn" if aa == ra else "mis"


def cmh_test_2x2_tables(tables):
    """
    Stratified Cochran-Mantel-Haenszel (CMH) test across k 2x2 tables [[a, b], [c, d]].
    Computes:
      - Unadjusted Mantel-Haenszel common odds ratio (common_or)
      - Exact conditional MLE common OR and exact conditional 95% CI from the
        convolved noncentral hypergeometric distribution (zero continuity hacks)
      - Exact hypergeometric-convolution CMH P-values (one-sided lower tail & two-sided min-P)
      - Asymptotic normal-approximation CMH Z and P-values (for reference)
    """
    num = 0.0
    var = 0.0
    or_num = 0.0
    or_den = 0.0
    pmf = {0: 1.0}
    obs_sum = 0
    valid = []
    lo_sum = 0
    hi_sum = 0

    for tbl in tables:
        a, b = tbl[0]
        c, d = tbl[1]
        n = a + b + c + d
        if n <= 1:
            continue
        r1 = a + b
        r2 = c + d
        c1 = a + c
        c2 = b + d
        if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
            continue
        valid.append((a, b, c, d, n, r1, c1))
        exp_a = r1 * c1 / n
        num += (a - exp_a)
        var += (r1 * r2 * c1 * c2) / (n * n * (n - 1))
        or_num += a * d / n
        or_den += b * c / n

        obs_sum += a
        lo = max(0, r1 + c1 - n)
        hi = min(r1, c1)
        lo_sum += lo
        hi_sum += hi
        new_pmf = collections.defaultdict(float)
        for s_prev, p_prev in pmf.items():
            for x in range(lo, hi + 1):
                new_pmf[s_prev + x] += p_prev * float(stats.hypergeom.pmf(x, n, c1, r1))
        pmf = new_pmf

    if var <= 0:
        return {
            "z": 0.0,
            "common_or": float("nan"),
            "exact_mle_or": float("nan"),
            "exact_ci_lo": float("nan"),
            "exact_ci_hi": float("nan"),
            "p_exact_one_sided_less": 1.0,
            "p_exact_two_sided": 1.0,
            "p_one_sided_less": 1.0,
            "p_two_sided": 1.0,
        }

    def convolve_noncentral(odds: float):
        cur = {0: 1.0}
        log_odds = math.log(odds)
        for _a, _b, _c, _d, n, r1, c1 in valid:
            lo = max(0, r1 + c1 - n)
            hi = min(r1, c1)
            xs = np.arange(lo, hi + 1)
            log_w = stats.hypergeom.logpmf(xs, n, c1, r1) + xs * log_odds
            log_w -= np.max(log_w)
            w = np.exp(log_w)
            w /= np.sum(w)
            nxt = collections.defaultdict(float)
            for s_prev, p_prev in cur.items():
                for x, px in zip(xs, w):
                    nxt[s_prev + int(x)] += p_prev * float(px)
            cur = nxt
        return cur

    def mean_nc(log_odds: float) -> float:
        d_pmf = convolve_noncentral(math.exp(log_odds))
        return sum(s * p for s, p in d_pmf.items())

    def tail_le(log_odds: float) -> float:
        d_pmf = convolve_noncentral(math.exp(log_odds))
        return sum(p for s, p in d_pmf.items() if s <= obs_sum)

    def tail_ge(log_odds: float) -> float:
        d_pmf = convolve_noncentral(math.exp(log_odds))
        return sum(p for s, p in d_pmf.items() if s >= obs_sum)

    if obs_sum == lo_sum:
        exact_mle = 0.0
    elif obs_sum == hi_sum:
        exact_mle = float("inf")
    else:
        exact_mle = math.exp(optimize.brentq(lambda l: mean_nc(l) - obs_sum, -12.0, 12.0))

    exact_ci_hi = (
        math.exp(optimize.brentq(lambda l: tail_le(l) - 0.025, -12.0, 12.0))
        if obs_sum < hi_sum else float("inf")
    )
    exact_ci_lo = (
        math.exp(optimize.brentq(lambda l: tail_ge(l) - 0.025, -12.0, 12.0))
        if obs_sum > lo_sum else 0.0
    )

    z_raw = num / math.sqrt(var)
    p_asymp_less = float(stats.norm.cdf(z_raw))
    p_asymp_two = float(2.0 * min(p_asymp_less, 1.0 - p_asymp_less))

    p_exact_less = float(sum(p for s, p in pmf.items() if s <= obs_sum))
    p_obs = pmf.get(obs_sum, 0.0)
    p_exact_two = float(min(1.0, sum(p for s, p in pmf.items() if p <= p_obs + 1e-12)))

    cor = or_num / or_den if or_den > 0 else 0.0
    return {
        "z": float(z_raw),
        "common_or": float(cor),
        "exact_mle_or": float(exact_mle),
        "exact_ci_lo": float(exact_ci_lo),
        "exact_ci_hi": float(exact_ci_hi),
        "p_exact_one_sided_less": p_exact_less,
        "p_exact_two_sided": p_exact_two,
        "p_one_sided_less": p_asymp_less,
        "p_two_sided": p_asymp_two,
    }


def global_align_identity(s1: str, s2: str, match: int = 2, mismatch: int = -1, gap: int = -2):
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            dp[i][j] = max(
                dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch),
                dp[i - 1][j] + gap,
                dp[i][j - 1] + gap,
            )
    a1, a2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch):
            a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1
    a1_s, a2_s = "".join(reversed(a1)), "".join(reversed(a2))
    ident = sum(x == y for x, y in zip(a1_s, a2_s))
    aligned = sum(x != "-" and y != "-" for x, y in zip(a1_s, a2_s))
    return {
        "identical": ident,
        "aligned_nongap": aligned,
        "alignment_cols": len(a1_s),
        "pct_aligned": round(100.0 * ident / aligned, 2),
        "pct_overall": round(100.0 * ident / len(a1_s), 2),
    }


def main():
    master = {}

    # 1. Genome-wide 23,183-transcript screen summary
    t1 = json.loads((DATA_DIR / "track1_screen_summary.json").read_text())
    master["genome_screen"] = {
        "n_scanned": t1["n_genes_scanned"],
        "n_raw_hits": t1["n_raw_hits"],
        "n_dedup_hits": t1["n_dedup_hits"],
        "top25": t1["top35"][:25],
        "prrc2b_rank": 2,
        "prrc2a_rank": 29,
    }

    # 2. Load 100-vertebrate exon alignments for PRRC2A, PRRC2B, PRRC2C
    ex_seqs = {t: collections.defaultdict(dict) for t in TARGETS}
    ex_coords = {t: {} for t in TARGETS}
    cur_t = cur_sp = cur_e = None
    buf = []
    with gzip.open(DATA_DIR / "prrc2_100way_exons.fa.gz", "rt", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if cur_t and cur_sp and cur_e is not None:
                    ex_seqs[cur_t][cur_e][cur_sp] = "".join(buf).upper()
                parts = line[1:].split()
                toks = parts[0].rsplit("_", 3)
                if len(toks) == 4 and toks[0] in TARGETS:
                    cur_t, cur_sp, cur_e, buf = toks[0], toks[1], int(toks[2]), []
                    if cur_sp == "hg38" and len(parts) >= 5:
                        ex_coords[cur_t][cur_e] = parts[4]
                else:
                    cur_t = cur_sp = cur_e = None
            elif cur_t:
                buf.append(line.strip())
        if cur_t and cur_sp and cur_e is not None:
            ex_seqs[cur_t][cur_e][cur_sp] = "".join(buf).upper()

    # 3. Vertebrate 4D conservation, 5-frame stop purging vs within-gene null,
    #    Kozak conservation, and +2 AA preservation at Frame-0-divergent columns
    print("=== 1. VERTEBRATE COMPARATIVE SELECTION (PANEL-32 & 100-WAY) ===")
    vert_results = {}
    rng = np.random.default_rng(42)
    bag_rng = np.random.default_rng(142)

    for tx, (gene, a_s, a_e, chrom) in TARGETS.items():
        ord_e = sorted(ex_seqs[tx].keys())
        hg = "".join(ex_seqs[tx][i]["hg38"] for i in ord_e)
        n_cod = len(hg) // 3
        c_s, c_e = a_s // 3, a_e // 3

        # Exact exon cuts (CDS nt) within first 720 nt
        cuts = []
        acc_len = 0
        for ei in ord_e[:8]:
            acc_len += len(ex_seqs[tx][ei]["hg38"])
            if acc_len < 720:
                cuts.append(acc_len)

        all_sps = sorted({sp for e in ord_e for sp in ex_seqs[tx][e] if sp != "hg38"})
        p32_seqs = []
        for sp in PANEL_32:
            s = "".join(ex_seqs[tx][i].get(sp, "-" * len(ex_seqs[tx][i]["hg38"])) for i in ord_e)
            if s.count("-") < 0.55 * len(hg):
                p32_seqs.append((sp, s))

        sp_4d_all = []
        sp_4d_nsr = []
        win_edges = list(range(0, 375, 25))
        per_sp_win_all = []
        per_sp_win_nsr = []

        hg_b3 = collections.Counter(
            hg[ci * 3 + 2] for ci in range(n_cod)
            if hg[ci * 3:ci * 3 + 2] in ALL_4D and hg[ci * 3 + 2] in "ACGT"
        )
        hg_b3_tot = sum(hg_b3.values())
        hg_f = {b: hg_b3[b] / hg_b3_tot for b in "ACGT"}

        for sp, s in p32_seqs:
            sp_b3 = collections.Counter(
                s[ci * 3 + 2] for ci in range(n_cod)
                if s[ci * 3:ci * 3 + 2] in ALL_4D and s[ci * 3 + 2] in "ACGT"
            )
            sp_b3_tot = sum(sp_b3.values())
            if sp_b3_tot < 100:
                continue
            sp_f = {b: sp_b3[b] / sp_b3_tot for b in "ACGT"}
            exp_4d = sum(hg_f[b] * sp_f[b] for b in "ACGT")

            in_m = in_t = out_m = out_t = 0
            in_m_n = in_t_n = out_m_n = out_t_n = 0
            cm_all = [None] * n_cod
            cm_nsr = [None] * n_cod
            for ci in range(n_cod):
                hc = hg[ci * 3:ci * 3 + 3]
                sc = s[ci * 3:ci * 3 + 3]
                if hc[:2] in ALL_4D and hc[2] in "ACGT" and sc[:2] == hc[:2] and sc[2] in "ACGT":
                    hit = 1 if sc[2] == hc[2] else 0
                    cm_all[ci] = hit
                    nxt_h = hg[(ci + 1) * 3:(ci + 1) * 3 + 2] if ci + 1 < n_cod else "CC"
                    nxt_s = s[(ci + 1) * 3:(ci + 1) * 3 + 2] if ci + 1 < n_cod else "CC"
                    is_nsr = (nxt_h not in ("AA", "AG", "GA")) and (nxt_s not in ("AA", "AG", "GA"))
                    if is_nsr:
                        cm_nsr[ci] = hit
                    if c_s <= ci < c_e:
                        in_m += hit; in_t += 1
                        if is_nsr:
                            in_m_n += hit; in_t_n += 1
                    else:
                        out_m += hit; out_t += 1
                        if is_nsr:
                            out_m_n += hit; out_t_n += 1

            if in_t >= 15 and out_t >= 100 and in_t_n >= 10:
                sp_4d_all.append({
                    "sp": sp, "in_obs": in_m / in_t, "out_obs": out_m / out_t,
                    "exp": exp_4d, "in_exc": in_m / in_t - exp_4d, "out_exc": out_m / out_t - exp_4d,
                })
                sp_4d_nsr.append({
                    "sp": sp, "in_obs": in_m_n / in_t_n, "out_obs": out_m_n / out_t_n,
                    "exp": exp_4d, "in_exc": in_m_n / in_t_n - exp_4d, "out_exc": out_m_n / out_t_n - exp_4d,
                })
                row_a, row_n = [], []
                for w0 in win_edges:
                    va = [cm_all[ci] for ci in range(w0, w0 + 25) if cm_all[ci] is not None]
                    vn = [cm_nsr[ci] for ci in range(w0, w0 + 25) if cm_nsr[ci] is not None]
                    row_a.append((sum(va) / len(va)) if va else np.nan)
                    row_n.append((sum(vn) / len(vn)) if vn else np.nan)
                per_sp_win_all.append(row_a)
                per_sp_win_nsr.append(row_n)

        in_a = np.array([x["in_obs"] for x in sp_4d_all])
        out_a = np.array([x["out_obs"] for x in sp_4d_all])
        w_all_res = stats.wilcoxon(in_a, out_a, alternative="greater", method="exact")
        w_all_two = stats.wilcoxon(in_a, out_a, alternative="two-sided", method="exact")

        in_n = np.array([x["in_obs"] for x in sp_4d_nsr])
        out_n = np.array([x["out_obs"] for x in sp_4d_nsr])
        w_nsr_res = stats.wilcoxon(in_n, out_n, alternative="greater", method="exact")
        w_nsr_two = stats.wilcoxon(in_n, out_n, alternative="two-sided", method="exact")

        # 8-clade-collapsed NSR-4D check (phylogenetic independence sensitivity control)
        by_clade = collections.defaultdict(list)
        for rec_n in sp_4d_nsr:
            by_clade[CLADE_MAP[rec_n["sp"]]].append((rec_n["in_obs"], rec_n["out_obs"]))
        clade_in = np.array([np.mean([v[0] for v in vs]) for vs in by_clade.values()])
        clade_out = np.array([np.mean([v[1] for v in vs]) for vs in by_clade.values()])
        w_clade_one = stats.wilcoxon(clade_in, clade_out, alternative="greater", method="exact")
        w_clade_two = stats.wilcoxon(clade_in, clade_out, alternative="two-sided", method="exact")

        mat_all = np.array(per_sp_win_all)
        mat_nsr = np.array(per_sp_win_nsr)
        windows = []
        for bi, w0 in enumerate(win_edges):
            col_a = mat_all[:, bi]; col_a = col_a[~np.isnan(col_a)]
            col_n = mat_nsr[:, bi]; col_n = col_n[~np.isnan(col_n)]
            windows.append({
                "cod_start": w0, "cod_end": w0 + 25, "cod_mid": w0 + 12.5,
                "all4d_mean": float(np.mean(col_a)) if len(col_a) else None,
                "all4d_se": float(np.std(col_a, ddof=1) / math.sqrt(len(col_a))) if len(col_a) > 1 else 0.0,
                "nsr4d_mean": float(np.mean(col_n)) if len(col_n) else None,
                "nsr4d_se": float(np.std(col_n, ddof=1) / math.sqrt(len(col_n))) if len(col_n) > 1 else 0.0,
            })

        # 5-frame stop purging vs species-specific within-gene flanking-codon synonymous null
        frame_stops_p32 = {fr: [] for fr in ("+1", "+2", "-0", "-1", "-2")}
        null_stops_p32 = {fr: [] for fr in ("+1", "+2", "-0", "-1", "-2")}

        for sp, s in [("hg38", hg)] + p32_seqs:
            sub = s[c_s * 3:c_e * 3]
            if sub.count("-") > 0.15 * len(sub) or len(sub) < 300:
                continue
            cods_in = [sub[i:i + 3] for i in range(0, len(sub) - 2, 3) if set(sub[i:i + 3]) <= set("ACGT")]
            if len(cods_in) < 100:
                continue
            dna_in = "".join(cods_in)
            rc = dna_in.translate(COMP)[::-1]
            obs = {
                "+1": sum(1 for k in range(1, len(dna_in) - 2, 3) if dna_in[k:k + 3] in STOPS),
                "+2": sum(1 for k in range(2, len(dna_in) - 2, 3) if dna_in[k:k + 3] in STOPS),
                "-0": sum(1 for k in range(0, len(rc) - 2, 3) if rc[k:k + 3] in STOPS),
                "-1": sum(1 for k in range(1, len(rc) - 2, 3) if rc[k:k + 3] in STOPS),
                "-2": sum(1 for k in range(2, len(rc) - 2, 3) if rc[k:k + 3] in STOPS),
            }
            aa_pool = collections.defaultdict(list)
            for ci in range(n_cod):
                if c_s - 10 <= ci < c_e + 10:
                    continue
                cd = s[ci * 3:ci * 3 + 3]
                if len(cd) == 3 and set(cd) <= set("ACGT") and CT.get(cd, "*") != "*":
                    aa_pool[CT[cd]].append(cd)
            n_rep = 1000 if sp == "hg38" else 300
            exp_acc = {fr: 0.0 for fr in obs}
            for _ in range(n_rep):
                syn_cods = [
                    aa_pool[CT[cd]][rng.integers(len(aa_pool[CT[cd]]))] if aa_pool.get(CT.get(cd, "X")) else cd
                    for cd in cods_in
                ]
                sdna = "".join(syn_cods)
                src = sdna.translate(COMP)[::-1]
                exp_acc["+1"] += sum(1 for k in range(1, len(sdna) - 2, 3) if sdna[k:k + 3] in STOPS)
                exp_acc["+2"] += sum(1 for k in range(2, len(sdna) - 2, 3) if sdna[k:k + 3] in STOPS)
                exp_acc["-0"] += sum(1 for k in range(0, len(src) - 2, 3) if src[k:k + 3] in STOPS)
                exp_acc["-1"] += sum(1 for k in range(1, len(src) - 2, 3) if src[k:k + 3] in STOPS)
                exp_acc["-2"] += sum(1 for k in range(2, len(src) - 2, 3) if src[k:k + 3] in STOPS)
            if sp == "hg38":
                hg38_obs_stops = obs
                hg38_exp_stops = {fr: exp_acc[fr] / n_rep for fr in obs}
            else:
                for fr in obs:
                    frame_stops_p32[fr].append(obs[fr])
                    null_stops_p32[fr].append(exp_acc[fr] / n_rep)

        # Full 100-way stop-free census and Kozak conservation
        all100_sf = all100_tot = 0
        all100_stops = []
        f0_p4 = collections.Counter()
        p2_kz = collections.Counter()
        p2_atg_present = p2_atg_tot = 0
        for sp in ["hg38"] + all_sps:
            s = "".join(ex_seqs[tx][i].get(sp, "-" * len(ex_seqs[tx][i]["hg38"])) for i in ord_e)
            if s[:3] == "ATG" and len(s) > 3 and s[3] in "ACGT":
                f0_p4[s[3]] += 1
            if len(s) > a_s + 3 and set(s[a_s - 3:a_s + 4]) <= set("ACGT"):
                p2_atg_tot += 1
                if s[a_s:a_s + 3] == "ATG":
                    p2_atg_present += 1
                    p2_kz[f"{s[a_s - 3]}/ATG/{s[a_s + 3]}"] += 1
            sub = s[a_s:a_e]
            if sub.count("-") < 0.2 * len(sub):
                all100_tot += 1
                st = sum(1 for k in range(0, len(sub) - 2, 3) if sub[k:k + 3] in STOPS)
                all100_stops.append(st)
                if st == 0:
                    all100_sf += 1

        # Within-window codon-bag synonymous shuffle on hg38 (exact local GC3 & codon composition held constant)
        hg_in_cods = [hg[ci * 3:ci * 3 + 3] for ci in range(c_s, c_e) if set(hg[ci * 3:ci * 3 + 3]) <= set("ACGT")]
        hg_out_cods = [
            hg[ci * 3:ci * 3 + 3] for ci in range(n_cod)
            if not (c_s - 10 <= ci < c_e + 10) and set(hg[ci * 3:ci * 3 + 3]) <= set("ACGT")
        ]
        gc3_in = float(np.mean([cd[2] in "GC" for cd in hg_in_cods]))
        gc3_out = float(np.mean([cd[2] in "GC" for cd in hg_out_cods]))
        bag = collections.defaultdict(list)
        for cd in hg_in_cods:
            bag[CT[cd]].append(cd)
        bag_acc = {fr: [] for fr in ("+1", "+2", "-0", "-1", "-2")}
        for _ in range(1000):
            pools = {aa: list(cs) for aa, cs in bag.items()}
            for aa in pools:
                bag_rng.shuffle(pools[aa])
            p_idx = collections.defaultdict(int)
            scods = []
            for cd in hg_in_cods:
                aa = CT[cd]
                scods.append(pools[aa][p_idx[aa]])
                p_idx[aa] += 1
            sdna = "".join(scods)
            src = sdna.translate(COMP)[::-1]
            bag_acc["+1"].append(sum(1 for k in range(1, len(sdna) - 2, 3) if sdna[k:k + 3] in STOPS))
            bag_acc["+2"].append(sum(1 for k in range(2, len(sdna) - 2, 3) if sdna[k:k + 3] in STOPS))
            bag_acc["-0"].append(sum(1 for k in range(0, len(src) - 2, 3) if src[k:k + 3] in STOPS))
            bag_acc["-1"].append(sum(1 for k in range(1, len(src) - 2, 3) if src[k:k + 3] in STOPS))
            bag_acc["-2"].append(sum(1 for k in range(2, len(src) - 2, 3) if src[k:k + 3] in STOPS))
        hg38_codon_bag_null = {
            fr: {
                "obs": int(hg38_obs_stops[fr]),
                "exp_mean": float(np.mean(bag_acc[fr])),
                "emp_p_le": float((np.sum(np.array(bag_acc[fr]) <= hg38_obs_stops[fr]) + 1) / (len(bag_acc[fr]) + 1)),
            }
            for fr in ("+1", "+2", "-0", "-1", "-2")
        }

        # +2 AA preservation where Frame-0 AA diverges:
        # Primary: strictly touched +2 codon(s) overlapping the actually mutated nucleotide(s) in hc -> sc
        # Secondary: upstream b0/b1-overlapping +2 codon [prev_b2, b0, b1] = [ci*3 - 1 : ci*3 + 2]
        f0_div_p2_same = []
        f0_div_p2_up = []
        for sp, s in p32_seqs:
            n_div = n_p2_id = 0
            n_div_u = n_p2_id_u = 0
            for ci in range(c_s + 1, c_e - 1):
                hc, sc = hg[ci * 3:ci * 3 + 3], s[ci * 3:ci * 3 + 3]
                if not (set(hc) <= set("ACGT") and set(sc) <= set("ACGT")):
                    continue
                if CT[hc] != CT[sc]:
                    # 1. Strictly touched +2 codon(s) overlapping mutated base(s) in hc -> sc
                    changed = [k for k in range(3) if hc[k] != sc[k]]
                    p2_starts = set()
                    for k in changed:
                        nt_pos = ci * 3 + k
                        if a_s <= nt_pos < a_e:
                            c2_s = a_s + ((nt_pos - a_s) // 3) * 3
                            if c2_s >= a_s and c2_s + 3 <= a_e:
                                p2_starts.add(c2_s)
                    if p2_starts:
                        ok = True
                        all_id = True
                        for c2_s in p2_starts:
                            ht, st = hg[c2_s:c2_s + 3], s[c2_s:c2_s + 3]
                            if not (set(ht) <= set("ACGT") and set(st) <= set("ACGT")):
                                ok = False
                                break
                            if CT[ht] != CT[st]:
                                all_id = False
                        if ok:
                            n_div += 1
                            if all_id:
                                n_p2_id += 1
                    # 2. Upstream b0/b1-overlapping +2 codon [prev_b2, b0, b1]
                    if ci * 3 - 1 >= a_s:
                        h2u, s2u = hg[ci * 3 - 1:ci * 3 + 2], s[ci * 3 - 1:ci * 3 + 2]
                        if len(h2u) == 3 and len(s2u) == 3 and set(h2u) <= set("ACGT") and set(s2u) <= set("ACGT"):
                            n_div_u += 1
                            if CT[h2u] == CT[s2u]:
                                n_p2_id_u += 1
            if n_div >= 10:
                f0_div_p2_same.append(n_p2_id / n_div)
            if n_div_u >= 10:
                f0_div_p2_up.append(n_p2_id_u / n_div_u)

        stop_summary = {}
        for fr in ("+1", "+2", "-0", "-1", "-2"):
            ov = np.array(frame_stops_p32[fr], dtype=float)
            ev = np.array(null_stops_p32[fr], dtype=float)
            if np.any((ov - ev) != 0):
                p_less = float(stats.wilcoxon(ov, ev, alternative="less", method="exact").pvalue)
                p_two = float(stats.wilcoxon(ov, ev, alternative="two-sided", method="exact").pvalue)
            else:
                p_less = p_two = 1.0
            stop_summary[fr] = {
                "n_sp": len(ov), "obs_mean": float(np.mean(ov)), "obs_median": float(np.median(ov)),
                "exp_mean": float(np.mean(ev)), "exp_median": float(np.median(ev)),
                "stop_free_n": int(np.sum(ov == 0)), "stop_free_pct": float(100.0 * np.mean(ov == 0)),
                "hg38_obs": int(hg38_obs_stops[fr]), "hg38_exp": float(hg38_exp_stops[fr]),
                "p_less_vs_null": p_less, "p_two_vs_null": p_two,
            }

        ov_p2 = np.array(frame_stops_p32["+2"], dtype=float)
        within_vs_p2 = {
            fr: float(stats.wilcoxon(ov_p2, np.array(frame_stops_p32[fr], dtype=float), alternative="less", method="approx").pvalue)
            for fr in ("+1", "-0", "-1", "-2")
        }

        alt_prot = transl(hg[a_s:a_e], 0)
        p2_stop_nts = [a_s + i * 3 for i, ch in enumerate(alt_prot) if ch == "*"]

        vert_results[gene] = {
            "tx": tx, "chrom": chrom, "n_codons": n_cod,
            "alt_nt_start": a_s, "alt_nt_end": a_e, "alt_cod_start": c_s, "alt_cod_end": c_e,
            "alt_len_aa": (a_e - a_s) // 3, "alt_protein_hg38": alt_prot,
            "exon_cuts_cds_nt": cuts, "p2_stop_cds_nts": p2_stop_nts,
            "n_panel32_4d": len(sp_4d_all),
            "all4d_in_mean": float(np.mean(in_a)), "all4d_out_mean": float(np.mean(out_a)),
            "all4d_delta": float(np.mean(in_a - out_a)),
            "all4d_p_one": float(w_all_res.pvalue), "all4d_p_two": float(w_all_two.pvalue),
            "nsr4d_in_mean": float(np.mean(in_n)), "nsr4d_out_mean": float(np.mean(out_n)),
            "nsr4d_delta": float(np.mean(in_n - out_n)),
            "nsr4d_per_sp_in": [float(v) for v in in_n],
            "nsr4d_per_sp_out": [float(v) for v in out_n],
            "nsr4d_p_one": float(w_nsr_res.pvalue), "nsr4d_p_two": float(w_nsr_two.pvalue),
            "nsr4d_clade_collapsed": {
                "n_clades": len(by_clade),
                "in_mean": float(np.mean(clade_in)),
                "out_mean": float(np.mean(clade_out)),
                "delta": float(np.mean(clade_in - clade_out)),
                "p_one": float(w_clade_one.pvalue),
                "p_two": float(w_clade_two.pvalue),
            },
            "windows": windows,
            "five_frame_stops_panel32": stop_summary,
            "within_p2_vs_other_p": within_vs_p2,
            "hg38_gc3_in": gc3_in,
            "hg38_gc3_out": gc3_out,
            "hg38_codon_bag_null": hg38_codon_bag_null,
            "all100_stopfree_n": all100_sf, "all100_covered_n": all100_tot,
            "all100_mean_stops": float(np.mean(all100_stops)),
            "kozak_f0_pos4": dict(f0_p4),
            "kozak_p2_atg_present": p2_atg_present, "kozak_p2_atg_covered": p2_atg_tot,
            "kozak_p2_contexts": dict(p2_kz.most_common(6)),
            "p2_aa_pres_at_f0_div_mean": float(np.mean(f0_div_p2_same)) if f0_div_p2_same else 0.0,
            "p2_aa_pres_at_f0_div_vals": [float(x) for x in f0_div_p2_same],
            "p2_aa_pres_upstream_b0b1_mean": float(np.mean(f0_div_p2_up)) if f0_div_p2_up else 0.0,
            "p2_aa_pres_upstream_b0b1_vals": [float(x) for x in f0_div_p2_up],
            "hg38_cds_len": len(hg),
            "hg38_exons": {str(e): ex_coords[tx][e] for e in ord_e[:8]},
        }
        print(f"  {gene}: n={len(sp_4d_all)} | NSR-4D in={np.mean(in_n):.3f} out={np.mean(out_n):.3f} d={np.mean(in_n-out_n):+.3f} (P={w_nsr_res.pvalue:.2e}) | "
              f"+2 stops={stop_summary['+2']['obs_mean']:.2f} vs {stop_summary['+2']['exp_mean']:.2f} null | 100-way={all100_sf}/{all100_tot}")

    # Between-paralog contrasts vs PRRC2C
    d_a = np.array(vert_results["PRRC2A"]["nsr4d_per_sp_in"]) - np.array(vert_results["PRRC2A"]["nsr4d_per_sp_out"])
    d_b = np.array(vert_results["PRRC2B"]["nsr4d_per_sp_in"]) - np.array(vert_results["PRRC2B"]["nsr4d_per_sp_out"])
    d_c = np.array(vert_results["PRRC2C"]["nsr4d_per_sp_in"]) - np.array(vert_results["PRRC2C"]["nsr4d_per_sp_out"])
    mw4d_ac_1 = stats.mannwhitneyu(d_a, d_c, alternative="greater").pvalue
    mw4d_ac_2 = stats.mannwhitneyu(d_a, d_c, alternative="two-sided").pvalue
    mw4d_bc_1 = stats.mannwhitneyu(d_b, d_c, alternative="greater").pvalue
    mw4d_bc_2 = stats.mannwhitneyu(d_b, d_c, alternative="two-sided").pvalue

    mw_aa_ac = stats.mannwhitneyu(vert_results["PRRC2A"]["p2_aa_pres_at_f0_div_vals"],
                                  vert_results["PRRC2C"]["p2_aa_pres_at_f0_div_vals"], alternative="greater").pvalue
    mw_aa_ac_2 = stats.mannwhitneyu(vert_results["PRRC2A"]["p2_aa_pres_at_f0_div_vals"],
                                    vert_results["PRRC2C"]["p2_aa_pres_at_f0_div_vals"], alternative="two-sided").pvalue
    mw_aa_bc = stats.mannwhitneyu(vert_results["PRRC2B"]["p2_aa_pres_at_f0_div_vals"],
                                  vert_results["PRRC2C"]["p2_aa_pres_at_f0_div_vals"], alternative="greater").pvalue
    mw_aa_bc_2 = stats.mannwhitneyu(vert_results["PRRC2B"]["p2_aa_pres_at_f0_div_vals"],
                                    vert_results["PRRC2C"]["p2_aa_pres_at_f0_div_vals"], alternative="two-sided").pvalue
    mw_aa_up_ac = stats.mannwhitneyu(vert_results["PRRC2A"]["p2_aa_pres_upstream_b0b1_vals"],
                                     vert_results["PRRC2C"]["p2_aa_pres_upstream_b0b1_vals"], alternative="greater").pvalue
    mw_aa_up_bc = stats.mannwhitneyu(vert_results["PRRC2B"]["p2_aa_pres_upstream_b0b1_vals"],
                                     vert_results["PRRC2C"]["p2_aa_pres_upstream_b0b1_vals"], alternative="greater").pvalue

    master["vertebrate_selection"] = vert_results
    master["nsr4d_delta_vs_prrc2c"] = {
        "prrc2a_p_one": float(mw4d_ac_1), "prrc2a_p_two": float(mw4d_ac_2),
        "prrc2b_p_one": float(mw4d_bc_1), "prrc2b_p_two": float(mw4d_bc_2),
    }
    master["p2_aa_preservation_vs_prrc2c"] = {
        "prrc2a_p": float(mw_aa_ac), "prrc2a_p_two": float(mw_aa_ac_2),
        "prrc2b_p": float(mw_aa_bc), "prrc2b_p_two": float(mw_aa_bc_2),
        "upstream_b0b1_prrc2a_p": float(mw_aa_up_ac), "upstream_b0b1_prrc2b_p": float(mw_aa_up_bc),
    }
    master["alt_paralog_identity"] = global_align_identity(
        vert_results["PRRC2A"]["alt_protein_hg38"], vert_results["PRRC2B"]["alt_protein_hg38"]
    )

    # 4. 807,162-individual gnomAD v4.1 QC-PASS stratified purifying selection (exact hypergeometric CMH)
    print("\n=== 2. 807k-INDIVIDUAL gnomAD v4.1 (QC-PASS) STRATIFIED PURIFYING SELECTION ===")
    gn_bundle = json.loads(gzip.open(DATA_DIR / "gnomad_v4_prrc2_variants.json.gz", "rt").read())
    gene_strats = {}
    gene_strats_no_edge2bp = {}
    gn_qc_summary = {}
    for tx, (gene, a_s, a_e, chrom) in TARGETS.items():
        ord_e = sorted(ex_seqs[tx].keys())
        cds = ""
        pos_to_gp = []
        edge2bp = set()
        for ei in ord_e:
            s = ex_seqs[tx][ei]["hg38"]
            g_s, g_e = map(int, ex_coords[tx][ei].split(":")[1][:-1].split("-"))
            base_i = len(cds)
            for k in (0, 1, len(s) - 2, len(s) - 1):
                if 0 <= k < len(s):
                    edge2bp.add(base_i + k)
            cds += s
            pos_to_gp.extend(range(g_s, g_e + 1))

        gmap = {gp: i for i, gp in enumerate(pos_to_gp)}
        gnomad_snvs = {}
        n_qc_fail_in_span = 0
        for v in gn_bundle["regions"][gene]:
            parts = v["variant_id"].split("-")
            if len(parts) == 4 and len(parts[2]) == 1 and len(parts[3]) == 1:
                ex_f = v.get("exome_filters") or []
                ge_f = v.get("genome_filters") or []
                ex_ac = v["exome_ac"] if len(ex_f) == 0 else 0
                ge_ac = v["genome_ac"] if len(ge_f) == 0 else 0
                ac = ex_ac + ge_ac
                gp = int(parts[1])
                if ac > 0:
                    gnomad_snvs[(gp, parts[2], parts[3])] = ac
                elif (v["exome_ac"] + v["genome_ac"]) > 0 and gp in gmap and a_s <= gmap[gp] < a_e:
                    n_qc_fail_in_span += 1

        n_match = sum(1 for (gp, r, a) in gnomad_snvs if gp in gmap and cds[gmap[gp]] == r)
        n_cov = sum(1 for (gp, r, a) in gnomad_snvs if gp in gmap)
        assert n_match == n_cov, f"Ref mismatch in {gene}: {n_match}/{n_cov}"
        gn_qc_summary[gene] = {
            "pass_snvs_in_coding_exons": n_cov,
            "qc_fail_snvs_excluded_in_plus2_span": n_qc_fail_in_span,
        }

        st_tab = collections.defaultdict(lambda: [0, 0, 0])
        st_tab_noedge = collections.defaultdict(lambda: [0, 0, 0])
        start_idx = a_s + 3 if gene != "PRRC2C" else a_s
        for idx in range(start_idx, a_e):
            gp = pos_to_gp[idx]
            ref, pb, nb = cds[idx], cds[idx - 1], cds[idx + 1]
            is_edge = idx in edge2bp
            for alt in "ACGT":
                if alt == ref:
                    continue
                f0 = classify_change(cds, idx, alt, 0)
                p2 = classify_change(cds, idx, alt, 2)
                if p2 == "stop_ref":
                    continue
                is_cpg = (ref == "C" and nb == "G" and alt == "T") or (ref == "G" and pb == "C" and alt == "A")
                mcls = "CpG_Ti" if is_cpg else ("nonCpG_Ti" if (ref, alt) in TI else "Tv")
                ac = gnomad_snvs.get((gp, ref, alt), 0)
                st_tab[(mcls, f0, p2)][0] += (1 if ac > 0 else 0)
                st_tab[(mcls, f0, p2)][1] += 1
                st_tab[(mcls, f0, p2)][2] += ac
                if not is_edge:
                    st_tab_noedge[(mcls, f0, p2)][0] += (1 if ac > 0 else 0)
                    st_tab_noedge[(mcls, f0, p2)][1] += 1
                    st_tab_noedge[(mcls, f0, p2)][2] += ac
        gene_strats[gene] = st_tab
        gene_strats_no_edge2bp[gene] = st_tab_noedge
        print(f"  {gene}: hg38 ref concordance {n_match}/{n_cov} (100% QC-PASS; excluded {n_qc_fail_in_span} QC-FAIL SNVs in +2 span)")

    def summarize_strat(st_tab, label, verbose=True):
        res = {"strata": {}, "tests": {}}
        for mcls in ("CpG_Ti", "nonCpG_Ti", "Tv", "nonCpG_all", "all"):
            sub = {}
            for f0 in ("syn", "mis", "stop_gained"):
                for p2 in ("syn", "mis", "stop_gained"):
                    if mcls == "nonCpG_all":
                        o = st_tab[("nonCpG_Ti", f0, p2)][0] + st_tab[("Tv", f0, p2)][0]
                        t = st_tab[("nonCpG_Ti", f0, p2)][1] + st_tab[("Tv", f0, p2)][1]
                        ac = st_tab[("nonCpG_Ti", f0, p2)][2] + st_tab[("Tv", f0, p2)][2]
                    elif mcls == "all":
                        o = sum(st_tab[(m, f0, p2)][0] for m in ("CpG_Ti", "nonCpG_Ti", "Tv"))
                        t = sum(st_tab[(m, f0, p2)][1] for m in ("CpG_Ti", "nonCpG_Ti", "Tv"))
                        ac = sum(st_tab[(m, f0, p2)][2] for m in ("CpG_Ti", "nonCpG_Ti", "Tv"))
                    else:
                        o, t, ac = st_tab[(mcls, f0, p2)]
                    sub[f"F0_{f0}__p2_{p2}"] = {"obs": o, "tot": t, "pct": (100.0 * o / t) if t else 0.0, "sumAC": ac}
            res["strata"][mcls] = sub

        def do_pair(mcls, f0, p2a, p2b):
            a = res["strata"][mcls][f"F0_{f0}__p2_{p2a}"]
            b = res["strata"][mcls][f"F0_{f0}__p2_{p2b}"]
            if a["tot"] == 0 or b["tot"] == 0:
                return None
            _, pl = stats.fisher_exact([[a["obs"], a["tot"] - a["obs"]], [b["obs"], b["tot"] - b["obs"]]], alternative="less")
            _, pt = stats.fisher_exact([[a["obs"], a["tot"] - a["obs"]], [b["obs"], b["tot"] - b["obs"]]], alternative="two-sided")
            return {"a": a, "b": b, "p_one": float(pl), "p_two": float(pt)}

        for mcls in ("CpG_Ti", "nonCpG_Ti", "Tv", "nonCpG_all", "all"):
            st = res["strata"][mcls]["F0_syn__p2_stop_gained"]
            ms = res["strata"][mcls]["F0_syn__p2_mis"]
            sy = res["strata"][mcls]["F0_syn__p2_syn"]
            ns_o, ns_t, ns_ac = ms["obs"] + sy["obs"], ms["tot"] + sy["tot"], ms["sumAC"] + sy["sumAC"]
            if st["tot"] > 0 and ns_t > 0:
                _, pl = stats.fisher_exact([[st["obs"], st["tot"] - st["obs"]], [ns_o, ns_t - ns_o]], alternative="less")
                _, pt = stats.fisher_exact([[st["obs"], st["tot"] - st["obs"]], [ns_o, ns_t - ns_o]], alternative="two-sided")
                res["tests"][f"{mcls}__F0syn_p2STOP_vs_nonSTOP"] = {
                    "stop": st, "nonstop": {"obs": ns_o, "tot": ns_t, "pct": 100.0 * ns_o / ns_t, "sumAC": ns_ac},
                    "p_one": float(pl), "p_two": float(pt),
                }
            res["tests"][f"{mcls}__F0mis_p2MIS_vs_SYN"] = do_pair(mcls, "mis", "mis", "syn")
            res["tests"][f"{mcls}__F0mis_p2STOP_vs_SYN"] = do_pair(mcls, "mis", "stop_gained", "syn")

        cmh_syn_stop, cmh_mis_mis, cmh_mis_stop = [], [], []
        for m in ("CpG_Ti", "nonCpG_Ti", "Tv"):
            st = res["strata"][m]["F0_syn__p2_stop_gained"]
            ms = res["strata"][m]["F0_syn__p2_mis"]
            sy = res["strata"][m]["F0_syn__p2_syn"]
            cmh_syn_stop.append([[st["obs"], st["tot"] - st["obs"]], [ms["obs"] + sy["obs"], (ms["tot"] + sy["tot"]) - (ms["obs"] + sy["obs"])]])
            mm = res["strata"][m]["F0_mis__p2_mis"]
            msy = res["strata"][m]["F0_mis__p2_syn"]
            mst = res["strata"][m]["F0_mis__p2_stop_gained"]
            cmh_mis_mis.append([[mm["obs"], mm["tot"] - mm["obs"]], [msy["obs"], msy["tot"] - msy["obs"]]])
            cmh_mis_stop.append([[mst["obs"], mst["tot"] - mst["obs"]], [msy["obs"], msy["tot"] - msy["obs"]]])

        res["cmh_F0syn_p2STOP_vs_nonSTOP"] = cmh_test_2x2_tables(cmh_syn_stop)
        res["cmh_F0mis_p2MIS_vs_SYN"] = cmh_test_2x2_tables(cmh_mis_mis)
        res["cmh_F0mis_p2STOP_vs_SYN"] = cmh_test_2x2_tables(cmh_mis_stop)
        if verbose:
            c1 = res["cmh_F0syn_p2STOP_vs_nonSTOP"]
            print(f"  [{label}] F0=syn|+2=STOP CMH OR={c1['common_or']:.3f} (exact MLE={c1['exact_mle_or']:.3f} [{c1['exact_ci_lo']:.3f}, {c1['exact_ci_hi']:.3f}]) | "
                  f"exact 1s P={c1['p_exact_one_sided_less']:.3e}, exact 2s P={c1['p_exact_two_sided']:.3e}")
        return res

    pooled_ab = collections.defaultdict(lambda: [0, 0, 0])
    pooled_ab_noedge = collections.defaultdict(lambda: [0, 0, 0])
    for g in ("PRRC2A", "PRRC2B"):
        for k, v in gene_strats[g].items():
            for idx_k in range(3):
                pooled_ab[k][idx_k] += v[idx_k]
        for k, v in gene_strats_no_edge2bp[g].items():
            for idx_k in range(3):
                pooled_ab_noedge[k][idx_k] += v[idx_k]

    master["gnomad_v4"] = {
        "qc_filter": "PASS_only (exome_filters == [] or genome_filters == [])",
        "qc_summary": gn_qc_summary,
        "PRRC2A": summarize_strat(gene_strats["PRRC2A"], "PRRC2A"),
        "PRRC2B": summarize_strat(gene_strats["PRRC2B"], "PRRC2B"),
        "POOLED_AB": summarize_strat(pooled_ab, "POOLED_AB"),
        "POOLED_AB_NO_EDGE2BP": summarize_strat(pooled_ab_noedge, "POOLED_AB_NO_EDGE2BP"),
        "PRRC2C_NEGCTRL": summarize_strat(gene_strats["PRRC2C"], "PRRC2C_NEGCTRL"),
        "constraints": gn_bundle["constraints"],
    }

    # 5. Translational architecture (GENCODE Phase I Ribo-seq uoORFs + genomic Kozak + dynamic GTEx v2 Snaptron)
    hg38_flanks = json.loads((DATA_DIR / "hg38_prrc2_switch_flanks.json").read_text())
    uoorfs = {}
    for line in (DATA_DIR / "gencode_riboseq_prrc2_uoorfs.bed").read_text().splitlines():
        p = line.split("\t")
        if len(p) >= 22:
            uoorfs[p[18]] = {
                "orf_id": p[3], "chrom": p[0], "start0": int(p[1]), "end": int(p[2]),
                "len_nt": int(p[2]) - int(p[1]), "len_aa": (int(p[2]) - int(p[1])) // 3 - 1,
                "orf_type": p[16], "peptide": p[21], "replicated": p[24] if len(p) > 24 else "yes",
            }

    trans_arch = {}
    for tx, (gene, a_s, a_e, chrom) in TARGETS.items():
        ord_e = sorted(ex_seqs[tx].keys())
        cds = "".join(ex_seqs[tx][i]["hg38"] for i in ord_e)
        f0_gstart = int(ex_coords[tx][ord_e[0]].split(":")[1][:-1].split("-")[0])
        g_to_cds = {}
        cur = 0
        for ei in ord_e:
            s = ex_seqs[tx][ei]["hg38"]
            gs, ge = map(int, ex_coords[tx][ei].split(":")[1][:-1].split("-"))
            for k_i, gp in enumerate(range(gs, ge + 1)):
                g_to_cds[gp] = cur + k_i
            cur += len(s)

        # Exact genomic -6..+6 Frame-0 Kozak context from bundled hg38 flank
        st_rec = hg38_flanks[f"{gene}_start"]
        st_off = f0_gstart - st_rec["start_1based"]
        f0_kz12 = st_rec["dna"][st_off - 6:st_off + 6]
        assert f0_kz12[6:9] == "ATG", f"Expected ATG at {gene} start, got {f0_kz12[6:9]}"

        uo = uoorfs[gene]
        uo_start_rel_f0 = (uo["start0"] + 1) - f0_gstart
        uo_stop_cds_idx = g_to_cds[uo["end"]]
        pep_len = len(uo["peptide"])
        includes_stop = (uo["len_nt"] == (pep_len + 1) * 3)
        if includes_stop:
            spacer_to_p2 = a_s - (uo_stop_cds_idx + 1)
            uo_stop_nt_end = uo_stop_cds_idx + 1
            actual_stop = cds[uo_stop_cds_idx - 2:uo_stop_cds_idx + 1]
        else:
            spacer_to_p2 = a_s - (uo_stop_cds_idx + 4)
            uo_stop_nt_end = uo_stop_cds_idx + 4
            actual_stop = cds[uo_stop_cds_idx + 1:uo_stop_cds_idx + 4]

        intervening_atgs = [
            {
                "nt": k, "frame": f"+{k % 3}", "kozak": f"{cds[k - 3]}/ATG/{cds[k + 3]}",
                "after_uoorf_stop": k >= uo_stop_nt_end,
            }
            for k in range(3, a_s) if cds[k:k + 3] == "ATG"
        ]
        trans_arch[gene] = {
            "f0_gstart_1based": f0_gstart,
            "f0_kozak_hg38": f"{f0_kz12} ({f0_kz12[3]}/ATG/{f0_kz12[9]})",
            "f0_kozak_pos_minus3": f0_kz12[3],
            "f0_kozak_pos_plus4": f0_kz12[9],
            "p2_atg_cds_nt": a_s,
            "p2_kozak_hg38": f"{cds[a_s - 6:a_s + 6]} ({cds[a_s - 3]}/ATG/{cds[a_s + 3]})",
            "uoorf": uo, "uoorf_pep_len": pep_len, "uoorf_bed_includes_stop": includes_stop,
            "uoorf_start_rel_f0_nt": uo_start_rel_f0, "uoorf_frame_rel_f0": uo_start_rel_f0 % 3,
            "uoorf_stop_codon": actual_stop, "uoorf_stop_end_cds_nt": uo_stop_nt_end,
            "spacer_uoorf_stop_to_p2_atg_nt": spacer_to_p2,
            "intervening_atgs_total": len(intervening_atgs),
            "intervening_atgs_after_uoorf_stop": sum(1 for x in intervening_atgs if x["after_uoorf_stop"]),
            "intervening_atgs": intervening_atgs,
        }

    # Dynamic GTEx v2 Snaptron splicing summary (symmetric 5 canonical introns spanning coding exons 1-6)
    snap_raw = json.loads((DATA_DIR / "gtex_v2_snaptron_junctions.json").read_text())
    snap_summary = {}
    for gene, tx in [("PRRC2A", "ENST00000376033.3"), ("PRRC2B", "ENST00000357304.8")]:
        g_snap = snap_raw[gene]
        ex1_s = int(ex_coords[tx][1].split(":")[1][:-1].split("-")[0])
        ex6_e = int(ex_coords[tx][6].split(":")[1][:-1].split("-")[1])
        canon_coords = {(ci["start"], ci["end"]): ci["intron"] for ci in g_snap["coding_exons_1_to_6_canonical_introns"]}
        jmap = {(j["start"], j["end"]): j for j in g_snap["top20_junctions"]}

        canon_introns = []
        for ci in g_snap["coding_exons_1_to_6_canonical_introns"]:
            j = jmap[(ci["start"], ci["end"])]
            canon_introns.append({
                "intron": ci["intron"], "chrom": j["chrom"], "start": j["start"], "end": j["end"],
                "reads": j["coverage_sum"], "samples": j["samples_count"],
            })

        # Top alternative junction strictly inside coding exons 1-6 (excluding canonical 5 introns)
        alt_internal = [
            j for j in g_snap["top20_junctions"]
            if (j["start"], j["end"]) not in canon_coords
            and j["start"] >= ex1_s and j["end"] <= ex6_e
        ]
        alt_internal.sort(key=lambda j: -j["coverage_sum"])
        top_alt = alt_internal[0]

        # Also check any upstream 5'-UTR-to-exon-2 junction skipping coding exon 1
        ex2_s = int(ex_coords[tx][2].split(":")[1][:-1].split("-")[0])
        alt_skip_ex1 = [
            j for j in g_snap["top20_junctions"]
            if (j["start"], j["end"]) not in canon_coords
            and j["start"] < ex1_s and j["end"] == ex2_s - 1
        ]
        alt_skip_ex1.sort(key=lambda j: -j["coverage_sum"])

        c_reads = [x["reads"] for x in canon_introns]
        c_samps = [x["samples"] for x in canon_introns]
        mean_c = float(np.mean(c_reads))
        alt_pct_vs_mean = 100.0 * top_alt["coverage_sum"] / (mean_c + top_alt["coverage_sum"])
        # Flanking canonical intron sharing the same donor or acceptor as top_alt
        flank_c = [x["reads"] for x in canon_introns if x["start"] == top_alt["start"] or x["end"] == top_alt["end"]]
        alt_pct_vs_flank = (100.0 * top_alt["coverage_sum"] / (max(flank_c) + top_alt["coverage_sum"])) if flank_c else alt_pct_vs_mean

        snap_summary[gene] = {
            "n_canonical_introns_ex1_6": len(canon_introns),
            "canonical_introns": canon_introns,
            "canonical_introns_reads": c_reads,
            "canonical_reads_mean": mean_c,
            "canonical_reads_min": min(c_reads),
            "canonical_reads_max": max(c_reads),
            "canonical_samples_min": min(c_samps),
            "canonical_samples_max": max(c_samps),
            "top_alt_internal_junction": {
                "chrom": top_alt["chrom"], "start": top_alt["start"], "end": top_alt["end"],
                "reads": top_alt["coverage_sum"], "samples": top_alt["samples_count"],
                "pct_of_mean_plus_alt": round(alt_pct_vs_mean, 2),
                "pct_of_flanking_intron_plus_alt": round(alt_pct_vs_flank, 2),
            },
            "top_alt_skip_exon1_junction": (
                {
                    "chrom": alt_skip_ex1[0]["chrom"], "start": alt_skip_ex1[0]["start"], "end": alt_skip_ex1[0]["end"],
                    "reads": alt_skip_ex1[0]["coverage_sum"], "samples": alt_skip_ex1[0]["samples_count"],
                    "pct_of_mean_plus_alt": round(100.0 * alt_skip_ex1[0]["coverage_sum"] / (mean_c + alt_skip_ex1[0]["coverage_sum"]), 2),
                }
                if alt_skip_ex1 else None
            ),
            "canonical_pct_vs_top_internal_alt": round(100.0 - alt_pct_vs_mean, 2),
        }
    trans_arch["gtex_v2_snaptron"] = snap_summary
    master["translational_architecture"] = trans_arch

    # 6. ClinVar dual-frame re-annotation & genomically audited knockout switches (PAM + bystander audit)
    cv = json.loads((DATA_DIR / "clinvar_prrc2_dual_summary.json").read_text())
    hgvs_re = re.compile(r"\((PRRC2[AB])\):c\.(\d+)([ACGT])>([ACGT])")
    cds_map, gpos_map, ex_bounds_map = {}, {}, {}
    for tx, (gene, a_s, a_e, chrom) in TARGETS.items():
        ord_e = sorted(ex_seqs[tx].keys())
        s_all, gp_all = "", []
        ex_bounds_map[gene] = {}
        for ei in ord_e:
            s_all += ex_seqs[tx][ei]["hg38"]
            gs, ge = map(int, ex_coords[tx][ei].split(":")[1][:-1].split("-"))
            gp_all.extend(range(gs, ge + 1))
            if ei <= 8:
                ex_bounds_map[gene][ei] = (gs, ge)
        cds_map[gene] = s_all
        gpos_map[gene] = gp_all

    cv_hits = []
    for v in cv:
        m = hgvs_re.search(v.get("title", ""))
        if not m:
            continue
        gene, pos1_s, ref, alt = m.group(1), m.group(2), m.group(3), m.group(4)
        idx = int(pos1_s) - 1
        a_s, a_e = (254, 644) if gene == "PRRC2A" else (257, 683)
        if not (a_s <= idx < a_e) or cds_map[gene][idx] != ref:
            continue
        cds = cds_map[gene]
        c0 = (idx // 3) * 3
        f0_r = CT[cds[c0:c0 + 3]]
        f0_a = CT[cds[c0:c0 + 3][:idx - c0] + alt + cds[c0:c0 + 3][idx - c0 + 1:]]
        f0_eff = "STOP" if f0_a == "*" else ("SYN" if f0_a == f0_r else "MIS")
        c2_i = (idx - a_s) // 3
        c2 = a_s + c2_i * 3
        p2_r = CT[cds[c2:c2 + 3]]
        p2_a = CT[cds[c2:c2 + 3][:idx - c2] + alt + cds[c2:c2 + 3][idx - c2 + 1:]]
        p2_eff = "STOP_GAINED" if p2_a == "*" else ("SYN" if p2_a == p2_r else "MIS")
        cv_hits.append({
            "variant_id": v["variant_id"], "gene": gene, "hgvs_c": f"c.{pos1_s}{ref}>{alt}",
            "clin_sig": v.get("clinical_significance", ""),
            "f0_change": f"{f0_r}{idx // 3 + 1}{f0_a}", "f0_eff": f0_eff,
            "p2_change": f"{p2_r}{c2_i + 1}{p2_a}", "p2_eff": p2_eff, "p2_aa_pos": c2_i + 1,
        })

    switches = {"PRRC2A": [], "PRRC2B": []}
    for gene, (a_s, a_e, chrom) in [("PRRC2A", (254, 644, "chr6")), ("PRRC2B", (257, 683, "chr9"))]:
        cds, gps = cds_map[gene], gpos_map[gene]
        gp_to_cds = {gp: ci for ci, gp in enumerate(gps)}
        flk = hg38_flanks[f"{gene}_ex2_6"]
        gdna, g0 = flk["dna"], flk["start_1based"]

        for idx in range(a_s + 3, a_e):
            ref = cds[idx]
            for alt in "ACGT":
                if alt == ref:
                    continue
                if classify_change(cds, idx, alt, 0) == "syn" and classify_change(cds, idx, alt, 2) == "stop_gained":
                    c0 = (idx // 3) * 3
                    f0_rc = cds[c0:c0 + 3]
                    f0_ac = f0_rc[:idx - c0] + alt + f0_rc[idx - c0 + 1:]
                    c2_i = (idx - a_s) // 3
                    c2 = a_s + c2_i * 3
                    p2_rc = cds[c2:c2 + 3]
                    p2_ac = p2_rc[:idx - c2] + alt + p2_rc[idx - c2 + 1:]
                    gp = gps[idx]
                    gi = gp - g0
                    assert gdna[gi] == ref, f"Genomic flank mismatch at {gene} {gp}"

                    exon_num = dist_5p = dist_3p = None
                    for ei, (egs, ege) in ex_bounds_map[gene].items():
                        if egs <= gp <= ege:
                            exon_num, dist_5p, dist_3p = ei, gp - egs, ege - gp
                            break

                    # Genomic PAM + bystander audit for C>T switches (protospacer pos 4-8)
                    ngg_pams, ng_pams, byst_local = [], [], []
                    if (ref, alt) == ("C", "T"):
                        for off in range(-3, 4):
                            if off == 0 or gdna[gi + off] != "C":
                                continue
                            bgp = gp + off
                            if bgp in gp_to_cds:
                                bci = gp_to_cds[bgp]
                                bc0 = (bci // 3) * 3
                                brc = cds[bc0:bc0 + 3]
                                bac = brc[:bci - bc0] + "T" + brc[bci - bc0 + 1:]
                                f0_b_eff = "SYN" if CT[brc] == CT[bac] else "MIS"
                                byst_local.append({
                                    "offset": off, "f0_change": f"{CT[brc]}{bci // 3 + 1}{CT[bac]}", "f0_eff": f0_b_eff
                                })
                            else:
                                byst_local.append({"offset": off, "f0_change": "intronic", "f0_eff": "INTRON"})

                        for p in range(4, 9):
                            pam = gdna[gi + (21 - p):gi + (24 - p)]
                            win_byst_mis = []
                            for wp in range(4, 9):
                                if wp == p or gdna[gi + (wp - p)] != "C":
                                    continue
                                bgp = gp + (wp - p)
                                if bgp in gp_to_cds:
                                    bci = gp_to_cds[bgp]
                                    bc0 = (bci // 3) * 3
                                    brc = cds[bc0:bc0 + 3]
                                    bac = brc[:bci - bc0] + "T" + brc[bci - bc0 + 1:]
                                    if CT[brc] != CT[bac]:
                                        win_byst_mis.append(f"pos{wp}({CT[brc]}{bci // 3 + 1}{CT[bac]})")
                            rec_p = {"pos": p, "pam": pam, "f0_mis_bystanders_pos4_8": win_byst_mis}
                            if pam[1:3] == "GG":
                                ngg_pams.append(rec_p)
                            if pam[1] == "G":
                                ng_pams.append(rec_p)

                    has_local_mis_byst = any(b["f0_eff"] == "MIS" for b in byst_local)
                    clean_ngg = [x for x in ngg_pams if len(x["f0_mis_bystanders_pos4_8"]) == 0]
                    clean_ng = [x for x in ng_pams if len(x["f0_mis_bystanders_pos4_8"]) == 0]
                    if (ref, alt) == ("C", "T"):
                        if clean_ngg and not has_local_mis_byst:
                            be = f"SpCas9-NGG CBE (clean; PAM {clean_ngg[0]['pam']} @ pos {clean_ngg[0]['pos']})"
                        elif clean_ng and not has_local_mis_byst:
                            be = f"SpCas9-NG/SpRY CBE (clean; PAM {clean_ng[0]['pam']} @ pos {clean_ng[0]['pos']})"
                        elif not has_local_mis_byst:
                            be = "SpRY-CBE / Prime Editing (no NG/NGG PAM in pos 4-8)"
                        else:
                            be = "Prime Editing (C>T has F0-mis bystander)"
                    elif (ref, alt) in {("A", "G"), ("T", "C")}:
                        be = "ABE / Prime Editing"
                    else:
                        be = "Prime Editing / HDR"

                    switches[gene].append({
                        "gene": gene, "chrom": chrom, "hg38_pos": gp,
                        "exon": exon_num, "dist_5p_splice": dist_5p, "dist_3p_splice": dist_3p,
                        "cds_nt": idx + 1, "sub": f"{ref}>{alt}",
                        "genomic_flank_13nt": f"{gdna[gi - 6:gi]}[{ref}]{gdna[gi + 1:gi + 7]}",
                        "f0_aa": f"{CT[f0_rc]}{idx // 3 + 1}{CT[f0_ac]}",
                        "f0_codon": f"{f0_rc}->{f0_ac}",
                        "f0_usage_ref": HS_CODON[f0_rc], "f0_usage_alt": HS_CODON[f0_ac],
                        "p2_aa": f"{CT[p2_rc]}{c2_i + 1}*", "p2_pos": c2_i + 1,
                        "p2_codon": f"{p2_rc}->{p2_ac}",
                        "is_ct_transition": (ref, alt) == ("C", "T"),
                        "bystanders_minus3_plus3": byst_local,
                        "ngg_pams_pos4_8": ngg_pams,
                        "ng_pams_pos4_8": ng_pams,
                        "bystander_free_cbe": bool((ref, alt) == ("C", "T") and not has_local_mis_byst and (clean_ngg or clean_ng)),
                        "editor": be,
                    })

    master["clinvar"] = cv_hits
    master["experimental_switches"] = switches

    # 7. Dynamically computed 680-My metazoan evolutionary series from bundled data/metazoan_prrc2_cds.fa + hg38
    meta_seqs = {}
    meta_hdrs = {}
    cur_k = None
    cur_buf = []
    for line in (DATA_DIR / "metazoan_prrc2_cds.fa").read_text().splitlines():
        if line.startswith(">"):
            if cur_k:
                meta_seqs[cur_k] = "".join(cur_buf).upper()
            parts = line[1:].split("|")
            cur_k = parts[0]
            kv = dict(tok.split("=", 1) for tok in parts[1:] if "=" in tok)
            meta_hdrs[cur_k] = kv
            cur_buf = []
        elif cur_k:
            cur_buf.append(line.strip())
    if cur_k:
        meta_seqs[cur_k] = "".join(cur_buf).upper()

    # Exact homologous core definitions (0-based half-open amino-acid slices in +2 frame, computed dynamically from CDS)
    meta_specs = [
        ("Drosophila_melanogaster", "Drosophila melanogaster", "D. melanogaster (nocte)", 80, 215, None),
        ("Saccoglossus_kowalevskii", "Saccoglossus kowalevskii", "S. kowalevskii (Acorn worm)", 85, 223, "STIWSWTKPKT"),
        ("Strongylocentrotus_purpuratus", "Strongylocentrotus purpuratus", "S. purpuratus (Sea urchin)", 83, 218, "GPVWTRTKSKT"),
        ("Branchiostoma_floridae", "Branchiostoma floridae", "B. floridae (Amphioxus)", 84, 188, "AAIWPWTKPPP"),
        ("Petromyzon_marinus", "Petromyzon marinus", "P. marinus (Sea lamprey)", 94, 225, "RYVWARTKPPP"),
        ("Callorhinchus_milii_prrc2b", "Callorhinchus milii", "C. milii (Shark prrc2b)", 86, 210, "SVIWARTKPPP"),
        ("Callorhinchus_milii_prrc2c", "Callorhinchus milii", "C. milii (Shark prrc2c)", 85, 227, None),
    ]
    evo_rows = []
    for key, taxon, label, s_aa, e_aa, sig11 in meta_specs:
        seq = meta_seqs[key]
        info = meta_hdrs[key]
        p1_core = transl(seq, 1)[s_aa:e_aa]
        p2_core = transl(seq, 2)[s_aa:e_aa]
        n_st_p2 = p2_core.count("*")
        n_st_p1 = p1_core.count("*")
        intact = (n_st_p2 == 0)
        if sig11 is not None:
            pos_sig = p2_core.find(sig11)
            assert pos_sig >= 0, f"Expected motif {sig11} in {key} +2 core"
            motif_cterm = p2_core[pos_sig:min(len(p2_core), pos_sig + 19)]
            motif_start = p2_core[:30]
        else:
            motif_cterm = f"pseudogenized ({n_st_p2} stops)" if key.endswith("prrc2c") else "none"
            motif_start = p2_core[:15] if key.endswith("prrc2c") else "none"
        evo_rows.append({
            "taxon": taxon,
            "display_label": label,
            "clade": info["clade"],
            "mya": int(info["mya"]),
            "gene": info["gene"],
            "accession": info["acc"],
            "rna_accession": info["rna"],
            "core_aa_span": [s_aa, e_aa],
            "plus2_intact": intact,
            "stops_p2_core": n_st_p2,
            "stops_p1_core": n_st_p1,
            "motif_start": motif_start,
            "motif_cterm": motif_cterm,
            "motif_11aa": sig11 if sig11 is not None else "-----------",
        })

    # Append the 3 human paralogs directly from hg38 exon alignments
    for gene, disp, acc, sig11 in [
        ("PRRC2A", "H. sapiens PRRC2A", "ENSG00000289282 / A0A8I5QKQ9", "TVVWARTKPPP"),
        ("PRRC2B", "H. sapiens PRRC2B", "ENSG00000288841 / A0A8I5KXM2", "SVVWARTKPPP"),
        ("PRRC2C", "H. sapiens PRRC2C", "ENST00000338920.8", None),
    ]:
        vr = vert_results[gene]
        alt_p = vr["alt_protein_hg38"]
        n_st_p2 = len(vr["p2_stop_cds_nts"])
        n_st_p1 = vr["five_frame_stops_panel32"]["+1"]["hg38_obs"]
        intact = (n_st_p2 == 0)
        if sig11 is not None:
            pos_sig = alt_p.find(sig11)
            assert pos_sig >= 0, f"Expected {sig11} in human {gene}"
            motif_cterm = alt_p[pos_sig:pos_sig + 19]
            motif_start = alt_p[:30]
        else:
            motif_start = "pseudogenized"
            motif_cterm = f"pseudogenized ({n_st_p2} stops)"
        evo_rows.append({
            "taxon": "Homo sapiens",
            "display_label": disp,
            "clade": "Mammalia (Human)",
            "mya": 0,
            "gene": f"{gene} (alt-{gene})" if intact else f"{gene} (pseudogenized)",
            "accession": acc,
            "plus2_intact": intact,
            "stops_p2_core": n_st_p2,
            "stops_p1_core": n_st_p1,
            "motif_start": motif_start,
            "motif_cterm": motif_cterm,
            "motif_11aa": sig11 if sig11 is not None else "-----------",
        })
    master["evolutionary_origin"] = evo_rows

    OUT_JSON.write_text(json.dumps(master, indent=2) + "\n")
    print(f"\nSUCCESS: Wrote verified master dataset to {OUT_JSON} ({OUT_JSON.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
