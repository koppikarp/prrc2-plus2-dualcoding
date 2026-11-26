"""Figure 3, Panel c: Stratified Cochran-Mantel-Haenszel odds ratios (PRRC2A, PRRC2B, Pooled vs PRRC2C)."""
import json
import math
import numpy as np
from style import DATA_DIR, PAL


def cmh_or_and_ci(strata_dict, mode):
    """Compute Haldane-corrected log-OR and 95% CI across (CpG_Ti, nonCpG_Ti, Tv) via Mantel-Haenszel variance."""
    # Use inverse-variance fixed-effects across the 3 mutational strata with 0.5 continuity for 0-cells in CI
    num = 0.0
    den = 0.0
    log_ors = []
    ws = []
    for m in ("CpG_Ti", "nonCpG_Ti", "Tv"):
        sub = strata_dict[m]
        if mode == "syn_stop":
            a = sub["F0_syn__p2_stop_gained"]["obs"]
            b = sub["F0_syn__p2_stop_gained"]["tot"] - a
            c = sub["F0_syn__p2_mis"]["obs"] + sub["F0_syn__p2_syn"]["obs"]
            d = (sub["F0_syn__p2_mis"]["tot"] + sub["F0_syn__p2_syn"]["tot"]) - c
        elif mode == "mis_stop":
            a = sub["F0_mis__p2_stop_gained"]["obs"]
            b = sub["F0_mis__p2_stop_gained"]["tot"] - a
            c = sub["F0_mis__p2_syn"]["obs"]
            d = sub["F0_mis__p2_syn"]["tot"] - c
        else:  # mis_mis
            a = sub["F0_mis__p2_mis"]["obs"]
            b = sub["F0_mis__p2_mis"]["tot"] - a
            c = sub["F0_mis__p2_syn"]["obs"]
            d = sub["F0_mis__p2_syn"]["tot"] - c
        if (a + b) == 0 or (c + d) == 0:
            continue
        # Haldane-Anscombe 0.5 adjustment ONLY if a cell is 0 for finite variance
        aa, bb, cc, dd = float(a), float(b), float(c), float(d)
        if min(aa, bb, cc, dd) == 0:
            aa += 0.5; bb += 0.5; cc += 0.5; dd += 0.5
        lor = math.log((aa * dd) / (bb * cc))
        var = 1.0 / aa + 1.0 / bb + 1.0 / cc + 1.0 / dd
        w = 1.0 / var
        log_ors.append(lor)
        ws.append(w)
    wsum = sum(ws)
    pooled_lor = sum(w * l for w, l in zip(ws, log_ors)) / wsum
    se = math.sqrt(1.0 / wsum)
    return math.exp(pooled_lor), math.exp(pooled_lor - 1.96 * se), math.exp(pooled_lor + 1.96 * se)


def draw(fig):
    ax = fig.add_subplot(111)
    gn = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["gnomad_v4"]

    rows = [
        ("PRRC2A (+2 STOP | F0 syn)", gn["PRRC2A"]["strata"], "syn_stop", PAL["prrc2a"]),
        ("PRRC2B (+2 STOP | F0 syn)", gn["PRRC2B"]["strata"], "syn_stop", PAL["prrc2b"]),
        ("Pooled A+B (+2 STOP | F0 syn)", gn["POOLED_AB"]["strata"], "syn_stop", PAL["prrc2a"]),
        ("Pooled A+B (+2 STOP | F0 mis)", gn["POOLED_AB"]["strata"], "mis_stop", PAL["prrc2a"]),
        ("Pooled A+B (+2 MIS | F0 mis)", gn["POOLED_AB"]["strata"], "mis_mis", PAL["uoorf"]),
        ("PRRC2C ctrl (+2 STOP | F0 syn)", gn["PRRC2C_NEGCTRL"]["strata"], "syn_stop", PAL["prrc2c"]),
        ("PRRC2C ctrl (+2 MIS | F0 mis)", gn["PRRC2C_NEGCTRL"]["strata"], "mis_mis", PAL["prrc2c"]),
    ]

    ys = np.arange(len(rows))[::-1]
    ax.axvline(1.0, ls="--", lw=0.8, color=PAL["ink"], zorder=1)
    ax.axhline(1.5, ls=":", lw=0.65, color=PAL["muted"], zorder=1)

    for y, (lbl, st, mode, col) in zip(ys, rows):
        or_v, lo, hi = cmh_or_and_ci(st, mode)
        ax.plot([lo, hi], [y, y], color=col, lw=1.35, zorder=3)
        ax.scatter([or_v], [y], s=26, color=col, edgecolors="white", linewidths=0.5, zorder=4)

    ax.set_xscale("log")
    ax.set_xlim(0.008, 4.5)
    ax.set_xticks([0.01, 0.05, 0.2, 1.0, 2.0])
    ax.set_xticklabels(["0.01", "0.05", "0.2", "1.0", "2.0"], fontsize=5.4)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.1)
    ax.set_xlabel("Stratified odds ratio (95% CI)", fontsize=6.0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
