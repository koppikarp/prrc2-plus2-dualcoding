"""Figure 4, Panel b: 100-vertebrate Kozak context lock at Frame-0 (+4 T) vs. +2 ([A/G]NNATGG)."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    vs = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["vertebrate_selection"]

    # Two side-by-side metrics across PRRC2A, PRRC2B, PRRC2C:
    # 1) % vertebrates with suboptimal +4 T at Frame-0 ATG
    # 2) % vertebrates with optimal purine (-3 A/G) + (+4 G) Kozak at +2 ATG
    genes = ["PRRC2A", "PRRC2B", "PRRC2C"]
    f0_t_pct = []
    p2_opt_pct = []
    for g in genes:
        fp4 = vs[g]["kozak_f0_pos4"]
        tot_f0 = sum(fp4.values())
        f0_t_pct.append(100.0 * fp4.get("T", 0) / max(1, tot_f0))
        kz = vs[g]["kozak_p2_contexts"]
        tot_p2 = sum(kz.values())
        opt = kz.get("A/ATG/G", 0) + kz.get("G/ATG/G", 0)
        p2_opt_pct.append(100.0 * opt / max(1, tot_p2))

    x = np.arange(3)
    w = 0.34
    b1 = ax.bar(x - w / 2, f0_t_pct, width=w, color=PAL["frame0"], label="F0 weak (+4T)", zorder=3)
    b2 = ax.bar(x + w / 2, p2_opt_pct, width=w,
                color=[PAL["prrc2a"], PAL["prrc2b"], PAL["prrc2c"]], label="+2 optimal (R..G)", zorder=3)

    for b, v in zip(b2, p2_opt_pct):
        ax.text(b.get_x() + b.get_width() / 2, v + 2.0, f"{v:.0f}%",
                ha="center", va="bottom", fontsize=5.1, color=PAL["ink"])

    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=5.6)
    ax.set_ylim(0, 130)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Vertebrates (%)", fontsize=6.1)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.grid(axis="y", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=5.0, loc="upper center", ncol=1, handlelength=1.1, labelspacing=0.2)
