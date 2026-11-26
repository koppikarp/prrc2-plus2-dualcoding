"""Figure 1, Panel a: To-scale reading-frame & exon map of human PRRC2A, PRRC2B, and PRRC2C."""
import json
import matplotlib.patches as mpatches
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())
    ax.set_xlim(-60, 760)
    ax.set_ylim(0.2, 3.6)
    ax.axis("off")

    rows = [
        ("PRRC2A\nchr6", 3.0, PAL["prrc2a"], PAL["prrc2a_lt"], -16, 71, 254, 644, 0, "130 aa"),
        ("PRRC2B\nchr9", 1.9, PAL["prrc2b"], PAL["prrc2b_lt"], -16, 74, 257, 683, 0, "142 aa"),
        ("PRRC2C\nchr1", 0.8, PAL["prrc2c"], PAL["prrc2c_lt"], -16, 71, 248, 650, 10, "10 stops"),
    ]

    # Nucleotide coordinate ruler along bottom
    ax.plot([0, 720], [0.28, 0.28], color=PAL["ink"], lw=0.65)
    for tick in (0, 100, 200, 300, 400, 500, 600, 700):
        ax.plot([tick, tick], [0.22, 0.34], color=PAL["ink"], lw=0.6)
        ax.text(tick, 0.16, str(tick), fontsize=5.5, color=PAL["ink"], ha="center", va="top")
    ax.text(360, -0.02, "CDS position (nt from Frame-0 ATG)", fontsize=6.2, color=PAL["ink"], ha="center", va="top")

    # Exact hg38 coding exon boundaries (CDS nt within 0..720) & PRRC2C stop positions from master_preprint_data.json
    vs = d["vertebrate_selection"]
    exon_cuts = {
        "PRRC2A\nchr6": vs["PRRC2A"]["exon_cuts_cds_nt"],
        "PRRC2B\nchr9": vs["PRRC2B"]["exon_cuts_cds_nt"],
        "PRRC2C\nchr1": vs["PRRC2C"]["exon_cuts_cds_nt"],
    }
    p2c_stops = vs["PRRC2C"]["p2_stop_cds_nts"]

    for lbl, y, col, col_lt, uo_s, uo_e, p2_s, p2_e, n_st, tag in rows:
        # Gene label on left
        ax.text(-58, y + 0.05, lbl, fontsize=6.2, fontweight="bold",
                color=col if n_st == 0 else PAL["muted"], ha="left", va="center")

        # 5' UTR leader line + mRNA backbone
        ax.plot([-25, 725], [y, y], color=PAL["muted"], lw=0.7, zorder=1)

        # Frame-0 BAT2 coding bar (nt 0..720)
        f0_bar = mpatches.FancyBboxPatch((0, y - 0.14), 720, 0.28,
                                         boxstyle="round,pad=0.0,rounding_size=0.04",
                                         fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.75, zorder=2)
        ax.add_patch(f0_bar)
        # Subtle exon junction ticks on Frame-0 bar
        for xc in exon_cuts[lbl]:
            ax.plot([xc, xc], [y - 0.14, y + 0.14], color="white", lw=0.9, zorder=3)

        # +2 uoORF bar (above Frame-0)
        uo_bar = mpatches.FancyBboxPatch((uo_s, y + 0.21), uo_e - uo_s, 0.19,
                                         boxstyle="round,pad=0.0,rounding_size=0.03",
                                         fc=PAL["uoorf"], ec="none", zorder=3)
        ax.add_patch(uo_bar)

        # +2 alt-PRRC2 bar (below Frame-0, or dashed gray for PRRC2C)
        if n_st == 0:
            p2_bar = mpatches.FancyBboxPatch((p2_s, y - 0.42), p2_e - p2_s, 0.22,
                                             boxstyle="round,pad=0.0,rounding_size=0.04",
                                             fc=col, ec="none", zorder=3)
            ax.add_patch(p2_bar)
            ax.text(p2_e + 8, y - 0.31, tag, fontsize=5.5, color=col, fontweight="bold", va="center")
        else:
            p2_bar = mpatches.FancyBboxPatch((p2_s, y - 0.42), p2_e - p2_s, 0.22,
                                             boxstyle="round,pad=0.0,rounding_size=0.04",
                                             fc="white", ec=PAL["prrc2c"], lw=0.8, ls="--", zorder=3)
            ax.add_patch(p2_bar)
            # Exact hg38 nucleotide positions of all 10 internal +2 stop codons inside PRRC2C
            for sx in p2c_stops:
                ax.plot([sx, sx], [y - 0.42, y - 0.20], color=PAL["prrc2c"], lw=0.9, zorder=4)
            ax.text(p2_e + 8, y - 0.31, tag, fontsize=5.5, color=PAL["muted"], va="center")

    # Minimal legend strip at top right (no box)
    ax.add_patch(mpatches.Rectangle((460, 3.40), 22, 0.14, fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.7))
    ax.text(487, 3.47, "Frame 0 (BAT2)", fontsize=5.6, color=PAL["ink"], va="center")
    ax.add_patch(mpatches.Rectangle((575, 3.40), 22, 0.14, fc=PAL["uoorf"], ec="none"))
    ax.text(602, 3.47, "+2 uoORF", fontsize=5.6, color=PAL["ink"], va="center")
    ax.add_patch(mpatches.Rectangle((660, 3.40), 22, 0.14, fc=PAL["prrc2a"], ec="none"))
    ax.text(687, 3.47, "+2 alt-PRRC2", fontsize=5.6, color=PAL["ink"], va="center")
