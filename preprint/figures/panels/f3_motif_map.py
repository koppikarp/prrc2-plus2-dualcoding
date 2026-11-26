"""Figure 3, Panel d: Domain/motif map of alt-PRRC2A and alt-PRRC2B with ClinVar variants & silent-F0/+2* switches."""
import json
import matplotlib.patches as mpatches
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "master_preprint_data.json").read_text())
    cv = d["clinvar"]
    sw = d["experimental_switches"]

    ax.set_xlim(-6, 150)
    ax.set_ylim(-0.15, 3.05)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])

    tracks = [
        ("PRRC2A", "alt-PRRC2A", 130, 1.85, PAL["prrc2a"], PAL["prrc2a_lt"],
         [(19, 30), (72, 86), (92, 116)]),
        ("PRRC2B", "alt-PRRC2B", 142, 0.58, PAL["prrc2b"], PAL["prrc2b_lt"],
         [(19, 30), (73, 87), (93, 117)]),
    ]

    for gene, lbl, plen, y, col, col_lt, motifs in tracks:
        # Exact 1-based inclusive residue coordinates: residue r centered at x=r, spanning [0.5 .. plen+0.5]
        ax.add_patch(mpatches.FancyBboxPatch((0.5, y - 0.14), plen, 0.28,
                                             boxstyle="round,pad=0.0,rounding_size=0.05",
                                             fc=col_lt, ec=col, lw=0.8, zorder=2))
        for ms, me in motifs:
            ax.add_patch(mpatches.Rectangle((ms - 0.5, y - 0.14), me - ms + 1, 0.28,
                                            fc=col, ec="none", alpha=0.82, zorder=3))
        ax.text(1, y + 0.45, f"{lbl} ({plen} aa)", fontsize=5.8, fontweight="bold", color=col, ha="left")

        cv_pos = [c["p2_aa_pos"] for c in cv if c["gene"] == gene and c["p2_eff"] in ("MIS", "STOP_GAINED")]
        for p in cv_pos:
            ax.plot([p, p], [y + 0.14, y + 0.31], color=PAL["clinvar"], lw=0.7, zorder=4)
            ax.scatter([p], [y + 0.33], s=11, color=PAL["clinvar"], edgecolors="none", zorder=5)

        clean_cbe_pos = [x["p2_pos"] for x in sw[gene] if x.get("bystander_free_cbe")]
        other_ct_pos = [x["p2_pos"] for x in sw[gene] if x.get("is_ct_transition") and not x.get("bystander_free_cbe")]
        for p in other_ct_pos:
            ax.scatter([p], [y - 0.28], marker="^", s=16, facecolors="white",
                       edgecolors=PAL["motif"], linewidths=0.75, zorder=5)
        for p in clean_cbe_pos:
            ax.scatter([p], [y - 0.28], marker="^", s=18, color=PAL["motif"],
                       edgecolors="none", zorder=6)

    ax.scatter([74], [2.88], s=14, color=PAL["clinvar"])
    ax.text(77.5, 2.88, "ClinVar (+2 mis)", fontsize=5.0, color=PAL["ink"], va="center")
    ax.scatter([74], [2.64], marker="^", s=16, color=PAL["motif"])
    ax.text(77.5, 2.64, "C>T silent-F0/+2* (clean CBE)", fontsize=5.0, color=PAL["ink"], va="center")
    ax.scatter([74], [2.40], marker="^", s=15, facecolors="white", edgecolors=PAL["motif"], linewidths=0.75)
    ax.text(77.5, 2.40, "C>T silent-F0/+2* (PE / bystander)", fontsize=5.0, color=PAL["ink"], va="center")

    ax.set_xticks([1, 25, 50, 75, 100, 125, 142])
    ax.set_xlabel("alt-PRRC2 amino-acid position", fontsize=6.1)
    ax.tick_params(axis="x", labelsize=5.5)
