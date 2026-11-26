"""Figure 3, Panel c: Exact conditional stratified CMH odds ratios and 95% CIs (PRRC2A, PRRC2B, Pooled vs PRRC2C)."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    gn = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["gnomad_v4"]

    rows = [
        ("PRRC2A (+2 STOP | F0 syn)", gn["PRRC2A"]["cmh_F0syn_p2STOP_vs_nonSTOP"], PAL["prrc2a"]),
        ("PRRC2B (+2 STOP | F0 syn)", gn["PRRC2B"]["cmh_F0syn_p2STOP_vs_nonSTOP"], PAL["prrc2b"]),
        ("Pooled A+B (+2 STOP | F0 syn)", gn["POOLED_AB"]["cmh_F0syn_p2STOP_vs_nonSTOP"], PAL["prrc2a"]),
        ("Pooled A+B (+2 STOP | F0 mis)", gn["POOLED_AB"]["cmh_F0mis_p2STOP_vs_SYN"], PAL["prrc2a"]),
        ("Pooled A+B (+2 MIS | F0 mis)", gn["POOLED_AB"]["cmh_F0mis_p2MIS_vs_SYN"], PAL["uoorf"]),
        ("PRRC2C ctrl (+2 STOP | F0 syn)", gn["PRRC2C_NEGCTRL"]["cmh_F0syn_p2STOP_vs_nonSTOP"], PAL["prrc2c"]),
        ("PRRC2C ctrl (+2 MIS | F0 mis)", gn["PRRC2C_NEGCTRL"]["cmh_F0mis_p2MIS_vs_SYN"], PAL["prrc2c"]),
    ]

    ys = np.arange(len(rows))[::-1]
    ax.axvline(1.0, ls="--", lw=0.8, color=PAL["ink"], zorder=1)
    ax.axhline(1.5, ls=":", lw=0.65, color=PAL["muted"], zorder=1)

    for y, (_lbl, cmh, col) in zip(ys, rows):
        or_v = cmh["exact_mle_or"]
        lo = cmh["exact_ci_lo"]
        hi = cmh["exact_ci_hi"]
        ax.plot([lo, hi], [y, y], color=col, lw=1.35, zorder=3)
        ax.scatter([or_v], [y], s=26, color=col, edgecolors="white", linewidths=0.5, zorder=4)

    ax.set_xscale("log")
    ax.set_xlim(0.00045, 4.5)
    ax.set_xticks([0.001, 0.01, 0.1, 1.0, 3.0])
    ax.set_xticklabels(["0.001", "0.01", "0.1", "1.0", "3.0"], fontsize=5.3)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.1)
    ax.set_xlabel("Exact conditional OR (95% CI)", fontsize=6.0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
