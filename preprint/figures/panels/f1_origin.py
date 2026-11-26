"""Figure 1, Panel b: 600-My evolutionary origin, +2 stop counts, and C-terminal motif conservation."""
import json
import matplotlib.patches as mpatches
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())
    evo = d["evolutionary_origin"]

    # Dynamically populated from master_preprint_data.json (computed from data/metazoan_prrc2_cds.fa + hg38)
    rows = [
        (r["display_label"], r["mya"], r["plus2_intact"], r["stops_p2_core"], r["motif_11aa"])
        for r in evo
    ]
    n = len(rows)
    ys = np.arange(n)[::-1]

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.7, n + 0.25)
    ax.axis("off")

    # Column headers (minimal)
    ax.text(2.0, n - 0.15, "Taxon / paralog", fontsize=5.8, fontweight="bold", color=PAL["ink"], ha="left", va="bottom")
    ax.text(46.0, n - 0.15, "Mya", fontsize=5.8, fontweight="bold", color=PAL["ink"], ha="center", va="bottom")
    ax.text(56.5, n - 0.15, "+2 stops", fontsize=5.8, fontweight="bold", color=PAL["ink"], ha="center", va="bottom")
    ax.text(80.5, n - 0.15, "Core +2 motif", fontsize=5.8, fontweight="bold", color=PAL["ink"], ha="center", va="bottom")
    ax.plot([1, 99], [n - 0.28, n - 0.28], color=PAL["ink"], lw=0.6)

    # Consensus reference for coloring conserved positions
    cons = "SVVWARTKPPP"
    # Group similar biochemical residues for soft highlight
    sim_groups = [{"V", "I", "L"}, {"S", "T"}, {"K", "R"}, {"A", "G"}, {"W"}, {"P"}]

    def is_cons(ch, ref):
        if ch == ref:
            return "exact"
        for g in sim_groups:
            if ch in g and ref in g:
                return "sim"
        return "diff"

    for idx, (lbl, mya, intact, nst, motif) in enumerate(rows):
        y = ys[idx]
        # Subtle row separator between clades
        if idx in (1, 5, 7):
            ax.plot([1, 99], [y + 0.5, y + 0.5], color=PAL["grid"], lw=0.55)

        txt_col = PAL["ink"] if intact else PAL["muted"]
        if "PRRC2A" in lbl:
            txt_col = PAL["prrc2a"]
        elif "PRRC2B" in lbl:
            txt_col = PAL["prrc2b"]

        ax.text(2.0, y, lbl, fontsize=5.4, color=txt_col,
                fontweight="bold" if "H. sapiens PRRC2A" in lbl or "H. sapiens PRRC2B" in lbl else "normal",
                va="center", ha="left")
        ax.text(46.0, y, str(mya) if mya > 0 else "0", fontsize=5.3, color=PAL["muted"], va="center", ha="center")

        # +2 stop badge (clean pill: solid terracotta/teal for 0, gray outline for >0)
        bcol = PAL["prrc2a"] if intact else PAL["prrc2c_lt"]
        if "PRRC2B" in lbl or "prrc2b" in lbl:
            bcol = PAL["prrc2b"]
        tcol = "white" if intact else PAL["muted"]
        pill = mpatches.FancyBboxPatch((52.5, y - 0.30), 8.0, 0.60,
                                       boxstyle="round,pad=0.0,rounding_size=0.18",
                                       fc=bcol, ec="none")
        ax.add_patch(pill)
        ax.text(56.5, y, str(nst), fontsize=5.3, color=tcol, fontweight="bold", va="center", ha="center")

        # Aligned 11-aa motif tiles (x = 64..98)
        x0 = 64.0
        cw = 3.05
        for j, ch in enumerate(motif):
            cx = x0 + j * cw
            if ch == "-":
                ax.text(cx + cw / 2, y, "·", fontsize=6.0, color=PAL["muted"], va="center", ha="center")
            else:
                st = is_cons(ch, cons[j])
                if st == "exact":
                    fc = PAL["prrc2b"] if ("PRRC2B" in lbl or "prrc2b" in lbl) else PAL["prrc2a"]
                    tc = "white"
                elif st == "sim":
                    fc = PAL["prrc2b_lt"] if ("PRRC2B" in lbl or "prrc2b" in lbl) else PAL["prrc2a_lt"]
                    tc = PAL["ink"]
                else:
                    fc = "white"
                    tc = PAL["muted"]
                ax.add_patch(mpatches.Rectangle((cx + 0.12, y - 0.34), cw - 0.24, 0.68, fc=fc, ec="none"))
                ax.text(cx + cw / 2, y, ch, fontsize=5.1, color=tc,
                        fontweight="bold" if st == "exact" else "normal",
                        va="center", ha="center")
