"""Figure 2, Panel a: Sliding-window strict non-stop-risk 4D conservation across Panel-32 vertebrates."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["vertebrate_selection"]

    # Subtle vertical dashed boundaries marking the +2 dual-coding span (codons 85..227)
    ax.axvline(85, ls="--", lw=0.65, color=PAL["muted"], zorder=1)
    ax.axvline(227, ls="--", lw=0.65, color=PAL["muted"], zorder=1)

    for gene, col, mk, lw in [
        ("PRRC2A", PAL["prrc2a"], "o", 1.5),
        ("PRRC2B", PAL["prrc2b"], "s", 1.5),
        ("PRRC2C", PAL["prrc2c"], "^", 1.15),
    ]:
        wins = [w for w in d[gene]["windows"] if w["nsr4d_mean"] is not None]
        xs = np.array([w["cod_mid"] for w in wins], dtype=float)
        ys = np.array([w["nsr4d_mean"] for w in wins], dtype=float)
        se = np.array([w["nsr4d_se"] for w in wins], dtype=float)
        ax.fill_between(xs, ys - se, ys + se, color=col, alpha=0.14, lw=0, zorder=2)
        ax.plot(xs, ys, marker=mk, ms=3.1, lw=lw, color=col, label=gene, zorder=3)

    ax.set_xlim(0, 375)
    ax.set_ylim(0.40, 1.00)
    ax.set_yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("Frame-0 position (codon)", fontsize=6.2)
    ax.set_ylabel("Non-stop-risk 4D identity", fontsize=6.2)
    ax.tick_params(labelsize=5.6)
    ax.legend(frameon=False, fontsize=5.6, loc="upper right", handlelength=1.4, labelspacing=0.25)
