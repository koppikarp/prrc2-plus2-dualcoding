"""Figure 3, Panel b: Ordered purifying-selection gradient among Frame-0-missense SNVs."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    gn = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["gnomad_v4"]
    ab = gn["POOLED_AB"]["strata"]
    c2 = gn["PRRC2C_NEGCTRL"]["strata"]

    # 3 groups: A+B non-CpG Ti, A+B Tv, and PRRC2C non-CpG Ti (control)
    groups = [
        ("Non-CpG Ti\n(A+B)", ab["nonCpG_Ti"]),
        ("Tv\n(A+B)", ab["Tv"]),
        ("Non-CpG Ti\n(2C ctrl)", c2["nonCpG_Ti"]),
    ]

    x = np.arange(len(groups))
    w = 0.24
    ax.axvline(1.5, ls=":", lw=0.75, color=PAL["muted"], zorder=1)

    for j, (key, lbl, col, off) in enumerate([
        ("F0_mis__p2_syn", "+2 SYN", PAL["frame0"], -w),
        ("F0_mis__p2_mis", "+2 MIS", PAL["uoorf"], 0.0),
        ("F0_mis__p2_stop_gained", "+2 STOP", PAL["prrc2a"], w),
    ]):
        vals = [g[1][key]["pct"] for g in groups]
        bars = ax.bar(x + off, vals, width=w * 0.88, color=col, label=lbl, zorder=3)
        # Gray out the 3rd group slightly via alpha? No, keep semantic color consistent so reader compares directly!
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 1.0, f"{v:.0f}%",
                    ha="center", va="bottom", fontsize=5.1, color=PAL["ink"])

    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups], fontsize=5.3)
    ax.set_ylim(0, 56)
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.set_ylabel("F0-mis. sites observed (%)", fontsize=6.1)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.grid(axis="y", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=5.2, loc="upper center", ncol=3, handlelength=1.0, columnspacing=0.7)
