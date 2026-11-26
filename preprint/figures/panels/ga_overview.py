"""Graphical Abstract: Minimal-text three-stage visual overview (zero gray boxes)."""
import json
import matplotlib.patches as mpatches
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 181)
    ax.set_ylim(0, 50)
    ax.axis("off")

    # Subtle vertical hairlines separating 3 visual columns (pure white background, zero gray boxes!)
    ax.plot([69, 69], [4, 46], color=PAL["grid"], lw=0.7)
    ax.plot([130, 130], [4, 46], color=PAL["grid"], lw=0.7)

    # ------------------------------------------------------------------
    # STAGE 1 (Left: x=2..66): Translational reinitiation & dual-coding architecture
    # ------------------------------------------------------------------
    ax.text(34, 46.5, "Dual-coding mRNA", fontsize=7.0, fontweight="bold",
            color=PAL["ink"], ha="center", va="top")

    # mRNA backbone line
    y_mrna = 25.0
    ax.plot([4, 65], [y_mrna, y_mrna], color=PAL["ink"], lw=1.0, zorder=1)

    # Frame-0 PRRC2A/B BAT2 coding bar (middle track)
    f0 = mpatches.FancyBboxPatch((15, y_mrna - 2.6), 48, 5.2,
                                 boxstyle="round,pad=0.2,rounding_size=1.0",
                                 fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.9, zorder=2)
    ax.add_patch(f0)
    ax.text(39, y_mrna, "PRRC2A / PRRC2B (Frame 0)", fontsize=5.6, color=PAL["frame0"],
            ha="center", va="center", fontweight="bold", zorder=3)

    # +2 uoORF (upper track, crossing weak F0 ATG)
    uo = mpatches.FancyBboxPatch((9, 32.5), 14, 4.4,
                                 boxstyle="round,pad=0.2,rounding_size=0.9",
                                 fc=PAL["uoorf"], ec="none", zorder=3)
    ax.add_patch(uo)
    ax.text(16, 34.7, "uoORF", fontsize=5.6, color="white", fontweight="bold", ha="center", va="center", zorder=4)

    # +2 alt-PRRC2A/B (lower track, nested inside BAT2)
    alt = mpatches.FancyBboxPatch((33, 13.0), 27, 5.0,
                                  boxstyle="round,pad=0.2,rounding_size=1.0",
                                  fc=PAL["prrc2a"], ec="none", zorder=3)
    ax.add_patch(alt)
    ax.text(46.5, 15.5, "alt-PRRC2A / B (+2)", fontsize=5.7, color="white", fontweight="bold",
            ha="center", va="center", zorder=4)

    # Weak +4 T tick at F0 start & optimal A/G..G tick at +2 start
    ax.plot([15, 15], [y_mrna + 2.6, 31.2], color=PAL["muted"], lw=0.7, ls=":")
    ax.text(15, 19.8, "+4T", fontsize=5.2, color=PAL["muted"], ha="center", va="top")

    # Curved reinitiation arrow from uoORF stop (x=23) to +2 start (x=33)
    ax.annotate("", xy=(33.0, 18.4), xytext=(23.2, 34.5),
                arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.25",
                                color=PAL["uoorf"], lw=1.1,
                                connectionstyle="arc3,rad=-0.28"), zorder=5)
    ax.text(30.2, 31.8, "183 nt", fontsize=5.5, color=PAL["uoorf"], fontweight="bold", ha="center")

    # Bottom scale label
    ax.text(34, 6.2, "130 / 142 aa", fontsize=5.5, color=PAL["prrc2a"], ha="center", fontweight="bold")

    # ------------------------------------------------------------------
    # STAGE 2 (Center: x=72..127): 600-My evolutionary timeline (cladogram)
    # ------------------------------------------------------------------
    ax.text(99.5, 46.5, "600-My evolution", fontsize=7.0, fontweight="bold",
            color=PAL["ink"], ha="center", va="top")

    # Clean horizontal cladogram branches
    # Trunk -> Protostome (top, gray) vs Deuterostome birth (terracotta) -> 2R WGD -> A, B (retained), C (lost)
    y_fly = 38.0
    y_deut = 24.0
    y_urch = 32.0
    y_lamp = 26.0
    y_a = 20.0
    y_b = 14.0
    y_c = 8.0

    # Root split
    ax.plot([73, 77], [30, 30], color=PAL["prrc2c"], lw=1.1)
    ax.plot([77, 77], [y_deut, y_fly], color=PAL["prrc2c"], lw=1.1)
    ax.plot([77, 114], [y_fly, y_fly], color=PAL["prrc2c"], lw=1.1)
    # Deuterostome +2 birth node
    ax.plot([77, 86], [y_deut, y_deut], color=PAL["prrc2a"], lw=1.3)
    ax.scatter([86], [y_deut], s=26, color=PAL["prrc2a"], zorder=5)
    ax.text(85.5, 19.8, "~580 Mya", fontsize=5.2, color=PAL["prrc2a"], ha="center", fontweight="bold")
    # Branches to urchin, lamprey, gnathostomes (2R)
    ax.plot([86, 86], [y_deut, y_urch], color=PAL["prrc2a"], lw=1.2)
    ax.plot([86, 114], [y_urch, y_urch], color=PAL["prrc2a"], lw=1.2)
    ax.plot([86, 94], [y_deut, y_deut], color=PAL["prrc2a"], lw=1.2)
    ax.plot([94, 94], [y_deut, y_lamp], color=PAL["prrc2a"], lw=1.2)
    ax.plot([94, 114], [y_lamp, y_lamp], color=PAL["prrc2a"], lw=1.2)
    # 2R WGD split at x=101 (connected cleanly from y_deut down to y_c)
    ax.plot([94, 101], [y_deut, y_deut], color=PAL["prrc2a"], lw=1.2)
    ax.plot([101, 101], [y_c, y_deut], color=PAL["ink"], lw=1.0)
    ax.plot([101, 114], [y_a, y_a], color=PAL["prrc2a"], lw=1.4)
    ax.plot([101, 114], [y_b, y_b], color=PAL["prrc2b"], lw=1.4)
    ax.plot([101, 114], [y_c, y_c], color=PAL["prrc2c"], lw=1.1, ls="--")
    # Loss cross on PRRC2C
    ax.scatter([106.5], [y_c], marker="x", s=22, color=PAL["prrc2c"], lw=1.3, zorder=5)

    for yv, lbl, col, dot_col in [
        (y_fly, "Fly nocte", PAL["muted"], PAL["prrc2c"]),
        (y_urch, "Urchin", PAL["ink"], PAL["prrc2a"]),
        (y_lamp, "Lamprey", PAL["ink"], PAL["prrc2a"]),
        (y_a, "PRRC2A", PAL["prrc2a"], PAL["prrc2a"]),
        (y_b, "PRRC2B", PAL["prrc2b"], PAL["prrc2b"]),
        (y_c, "PRRC2C", PAL["muted"], PAL["prrc2c"]),
    ]:
        ax.scatter([114], [yv], s=14, color=dot_col, zorder=5)
        ax.text(116.2, yv, lbl, fontsize=5.3, color=col, va="center", ha="left",
                fontweight="bold" if lbl in ("PRRC2A", "PRRC2B") else "normal")

    # ------------------------------------------------------------------
    # STAGE 3 (Right: x=133..179): Human gnomAD v4 purifying selection mini-chart
    # ------------------------------------------------------------------
    ax.text(156, 46.5, "Human selection (gnomAD)", fontsize=7.0, fontweight="bold",
            color=PAL["ink"], ha="center", va="top")

    # Mini inset axes for F0=syn non-CpG SNVs observed (%) — loaded dynamically from master_preprint_data.json
    gn = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["gnomad_v4"]
    t_ab = gn["POOLED_AB"]["tests"]["nonCpG_all__F0syn_p2STOP_vs_nonSTOP"]
    t_2c = gn["PRRC2C_NEGCTRL"]["tests"]["nonCpG_all__F0syn_p2STOP_vs_nonSTOP"]
    v_ns = t_ab["nonstop"]["pct"]
    v_st = t_ab["stop"]["pct"]
    v_2c = t_2c["stop"]["pct"]
    st_lbl = f"{t_ab['stop']['obs']}/{t_ab['stop']['tot']}"

    ax_bar = fig.add_axes([0.755, 0.21, 0.225, 0.56])
    cats = ["A+B\nnon-STOP", "A+B\n+2 STOP", "2C\n+2 STOP"]
    vals = [v_ns, v_st, v_2c]
    cols = [PAL["frame0"], PAL["prrc2a"], PAL["prrc2c"]]
    xs = np.arange(3)
    bars = ax_bar.bar(xs, vals, width=0.58, color=cols, edgecolor="none")
    ax_bar.set_ylim(0, 24)
    ax_bar.set_yticks([0, 10, 20])
    ax_bar.set_xticks(xs)
    ax_bar.set_xticklabels(cats, fontsize=5.1)
    ax_bar.set_ylabel("Observed (%)", fontsize=5.4, labelpad=1.5)
    ax_bar.tick_params(axis="both", labelsize=5.1, length=2.2, pad=1.5)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    for b, v, txt in zip(bars, vals, [f"{v_ns:.1f}%", st_lbl, f"{v_2c:.1f}%"]):
        ax_bar.text(b.get_x() + b.get_width() / 2, v + 0.9, txt,
                    ha="center", va="bottom", fontsize=5.2,
                    fontweight="bold" if txt == st_lbl else "normal",
                    color=PAL["prrc2a"] if txt == st_lbl else PAL["ink"])
