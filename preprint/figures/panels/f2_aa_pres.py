"""Figure 2, Panel d: +2 amino-acid preservation at Frame-0-divergent vertebrate sites."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["vertebrate_selection"]

    genes = ["PRRC2A", "PRRC2B", "PRRC2C"]
    cols = [PAL["prrc2a"], PAL["prrc2b"], PAL["prrc2c"]]
    rng = np.random.default_rng(7)

    for i, (g, col) in enumerate(zip(genes, cols)):
        vals = 100.0 * np.array(d[g]["p2_aa_pres_at_f0_div_vals"])
        jitter = rng.uniform(-0.14, 0.14, size=len(vals))
        ax.scatter(i + jitter, vals, s=14, color=col, alpha=0.72, edgecolors="none", zorder=3)
        med = float(np.median(vals))
        q25, q75 = np.percentile(vals, [25, 75])
        ax.plot([i, i], [q25, q75], color=PAL["ink"], lw=1.2, zorder=4)
        ax.plot([i - 0.22, i + 0.22], [med, med], color=PAL["ink"], lw=1.6, zorder=5)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(genes, fontsize=5.8)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(14, 97)
    ax.set_yticks([20, 40, 60, 80])
    ax.set_ylabel("+2 AA preserved (%)", fontsize=6.2)
    ax.tick_params(axis="y", labelsize=5.6)
    ax.grid(axis="y", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
