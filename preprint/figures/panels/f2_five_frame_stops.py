"""Figure 2, Panel c: 5-frame stop-codon purging across Panel-32 vertebrates vs. within-gene synonymous null."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["vertebrate_selection"]

    frames = ["+1", "+2", "-0", "-1", "-2"]
    x = np.arange(len(frames))
    w = 0.25

    for idx, (gene, col, off) in enumerate([
        ("PRRC2A", PAL["prrc2a"], -w),
        ("PRRC2B", PAL["prrc2b"], 0.0),
        ("PRRC2C", PAL["prrc2c"], w),
    ]):
        fs = d[gene]["five_frame_stops_panel32"]
        obs = np.array([fs[fr]["obs_mean"] for fr in frames])
        exp = np.array([fs[fr]["exp_mean"] for fr in frames])
        ax.bar(x + off, obs, width=w * 0.88, color=col, edgecolor="none", label=gene, zorder=3)
        # Clean horizontal dash markers for within-gene flanking-codon synonymous null expectation
        for xi, ei in zip(x + off, exp):
            ax.plot([xi - w * 0.42, xi + w * 0.42], [ei, ei], color=PAL["ink"], lw=1.15, zorder=5)

    # Dummy legend item for null marker
    ax.plot([], [], color=PAL["ink"], lw=1.15, label="Syn. null")

    ax.set_xticks(x)
    ax.set_xticklabels(frames, fontsize=6.0)
    ax.set_ylim(0, 12.2)
    ax.set_yticks([0, 3, 6, 9, 12])
    ax.set_xlabel("Alternative reading frame", fontsize=6.2)
    ax.set_ylabel("Mean stops per species", fontsize=6.2)
    ax.tick_params(labelsize=5.6)
    ax.grid(axis="y", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=5.4, loc="upper right", ncol=2, handlelength=1.2, columnspacing=0.8)
