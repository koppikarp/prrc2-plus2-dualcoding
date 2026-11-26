"""Figure 4, Panel a: Clean schematic of the conserved uoORF-mediated reinitiation circuit in PRRC2A/B."""
import matplotlib.patches as mpatches
from style import PAL


def draw(fig):
    ax = fig.add_subplot(111)
    ax.set_xlim(-45, 415)
    ax.set_ylim(-0.35, 3.05)
    ax.axis("off")

    # mRNA line (5' cap -> CDS)
    y_m = 1.35
    ax.plot([-32, 405], [y_m, y_m], color=PAL["ink"], lw=1.0, zorder=1)
    # m7G cap dot
    ax.scatter([-32], [y_m], s=32, color=PAL["ink"], zorder=4)
    ax.text(-32, y_m + 0.25, "5'", fontsize=5.8, fontweight="bold", color=PAL["ink"], ha="center")

    # Frame-0 PRRC2A/B CDS starting at nt 0 (weak +4 T Kozak)
    f0 = mpatches.FancyBboxPatch((0, y_m - 0.22), 400, 0.44,
                                 boxstyle="round,pad=0.0,rounding_size=0.06",
                                 fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.85, zorder=2)
    ax.add_patch(f0)
    ax.text(105, y_m, "PRRC2A / PRRC2B (Frame 0)", fontsize=6.0, fontweight="bold",
            color=PAL["frame0"], ha="center", va="center", zorder=3)

    # +2 uoORF: nt -16 to +71/+74 (upper track)
    uo = mpatches.FancyBboxPatch((-16, 2.02), 89, 0.38,
                                 boxstyle="round,pad=0.0,rounding_size=0.06",
                                 fc=PAL["uoorf"], ec="none", zorder=3)
    ax.add_patch(uo)
    ax.text(28.5, 2.21, "uoORF (+2, 28/29 aa)", fontsize=5.7, fontweight="bold",
            color="white", ha="center", va="center", zorder=4)

    # +2 alt-PRRC2A/B: starting at nt 254/257 (lower track)
    alt = mpatches.FancyBboxPatch((255, 0.36), 145, 0.40,
                                  boxstyle="round,pad=0.0,rounding_size=0.06",
                                  fc=PAL["prrc2a"], ec="none", zorder=3)
    ax.add_patch(alt)
    ax.text(327.5, 0.56, "alt-PRRC2A / B (+2)", fontsize=5.9, fontweight="bold",
            color="white", ha="center", va="center", zorder=4)

    # Start codon vertical ticks & minimal labels
    ax.plot([-16, -16], [y_m, 2.02], color=PAL["uoorf"], lw=0.8, ls=":")
    ax.text(-16, 2.52, "-16 nt", fontsize=5.4, color=PAL["uoorf"], fontweight="bold", ha="center")

    ax.plot([0, 0], [0.92, y_m - 0.22], color=PAL["frame0"], lw=0.8, ls=":")
    ax.text(0, 0.85, "0 nt (+4T)", fontsize=5.4, color=PAL["frame0"], ha="center", va="top")

    ax.plot([255, 255], [0.76, y_m - 0.22], color=PAL["prrc2a"], lw=0.8, ls=":")
    ax.text(255, 0.24, "+254/257 nt (A/G..G)", fontsize=5.4, color=PAL["prrc2a"], fontweight="bold", ha="center", va="top")

    # 183 nt reinitiation spacer bracket & scanning arrow (from uoORF stop at 73 to +2 start at 255)
    ax.annotate("", xy=(255, 2.21), xytext=(75, 2.21),
                arrowprops=dict(arrowstyle="<->", color=PAL["ink"], lw=0.85))
    ax.text(165, 2.34, "183 nt spacer", fontsize=5.8, fontweight="bold", color=PAL["ink"], ha="center", va="bottom")

    # Curved reinitiation arc dropping into +2 start
    ax.annotate("", xy=(254, 0.80), xytext=(74, 2.02),
                arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.28",
                                color=PAL["uoorf"], lw=1.2, ls="--",
                                connectionstyle="arc3,rad=-0.18"), zorder=5)

    # Bottom scale ruler
    ax.plot([0, 400], [-0.18, -0.18], color=PAL["muted"], lw=0.6)
    for t in (0, 100, 200, 300, 400):
        ax.plot([t, t], [-0.24, -0.12], color=PAL["muted"], lw=0.6)
        ax.text(t, -0.29, str(t), fontsize=5.2, color=PAL["muted"], ha="center", va="top")
