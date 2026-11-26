"""Figure 2, Panel b: Per-species paired inside vs. outside non-stop-risk 4D identity."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["vertebrate_selection"]

    groups = [
        ("PRRC2A", 0.0, PAL["prrc2a"]),
        ("PRRC2B", 1.6, PAL["prrc2b"]),
        ("PRRC2C", 3.2, PAL["prrc2c"]),
    ]
    xticks = []
    xlabels = []
    for gene, x0, col in groups:
        ins = np.array(d[gene]["nsr4d_per_sp_in"])
        outs = np.array(d[gene]["nsr4d_per_sp_out"])
        x_in, x_out = x0, x0 + 0.65
        for yi, yo in zip(ins, outs):
            ax.plot([x_in, x_out], [yi, yo], color=col, alpha=0.28, lw=0.6, zorder=1)
        ax.scatter(np.full_like(ins, x_in), ins, s=11, color=col, alpha=0.75, edgecolors="none", zorder=2)
        ax.scatter(np.full_like(outs, x_out), outs, s=11, color=col, alpha=0.45, edgecolors="none", zorder=2)
        # Bold mean markers + connecting line
        mi, mo = float(np.mean(ins)), float(np.mean(outs))
        ax.plot([x_in, x_out], [mi, mo], color=PAL["ink"], lw=1.35, zorder=4)
        ax.scatter([x_in, x_out], [mi, mo], s=22, color=col, edgecolors=PAL["ink"], linewidths=0.7, zorder=5)
        ax.text(x0 + 0.325, 1.055, gene, fontsize=5.8, fontweight="bold", color=col, ha="center", va="top")
        xticks.extend([x_in, x_out])
        xlabels.extend(["In", "Out"])

    ax.set_xlim(-0.4, 4.25)
    ax.set_ylim(0.28, 1.08)
    ax.set_yticks([0.3, 0.5, 0.7, 0.9])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=5.5)
    ax.set_ylabel("Non-stop-risk 4D identity", fontsize=6.2)
    ax.tick_params(axis="y", labelsize=5.6)
