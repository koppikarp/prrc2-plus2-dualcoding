"""Figure 4, Panel c: GTEx v2 Snaptron splice-junction read support across ~19,000 human samples."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    gt = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["translational_architecture"]["gtex_v2_snaptron"]

    # Plot mean canonical intron reads vs top alternative skipping junction reads on log10 scale!
    x = np.arange(2)
    w = 0.34
    canon = [
        float(np.mean(gt["PRRC2A"]["canonical_introns_reads"])),
        float(np.mean(gt["PRRC2B"]["canonical_introns_reads"])),
    ]
    alt = [
        float(gt["PRRC2A"]["top_alt_internal_junction"]["reads"]),
        float(gt["PRRC2B"]["top_alt_internal_junction"]["reads"]),
    ]

    ax.bar(x - w / 2, canon, width=w, color=[PAL["prrc2a"], PAL["prrc2b"]], label="Canonical", zorder=3)
    ax.bar(x + w / 2, alt, width=w, color=PAL["prrc2c"], label="Top alt. skip", zorder=3)

    # Individual canonical intron dots overlaid
    for xi, gene, col in [(0, "PRRC2A", PAL["prrc2a"]), (1, "PRRC2B", PAL["prrc2b"])]:
        vals = gt[gene]["canonical_introns_reads"]
        ax.scatter(np.full(len(vals), xi - w / 2), vals, s=10, color="white",
                   edgecolors=PAL["ink"], linewidths=0.55, zorder=5)

    for xi, (cv, av) in enumerate(zip(canon, alt)):
        pct = 100.0 * av / (cv + av)
        ax.text(xi + w / 2, av * 1.35, f"{pct:.1f}%", ha="center", va="bottom",
                fontsize=5.2, color=PAL["muted"])

    ax.set_yscale("log")
    ax.set_ylim(1e3, 1e7)
    ax.set_yticks([1e3, 1e4, 1e5, 1e6, 1e7])
    ax.set_xticks(x)
    ax.set_xticklabels(["PRRC2A", "PRRC2B"], fontsize=5.8)
    ax.set_ylabel("GTEx v2 junction reads", fontsize=6.1)
    ax.tick_params(axis="y", labelsize=5.4)
    ax.grid(axis="y", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=5.1, loc="upper right", handlelength=1.1, labelspacing=0.2)
