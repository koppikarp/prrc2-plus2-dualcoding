"""Figure 4, Panel d: Conserved GENCODE Ribo-seq +2 uoORF peptides across PRRC2A, PRRC2B, and PRRC2C."""
import json
import matplotlib.patches as mpatches
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    ta = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["translational_architecture"]

    # Align the 3 replicated uoORF peptides (28, 29, 28 aa -> gap 1 aa in A/C at pos 24 to align C-term QPV)
    seqs = [
        ("PRRC2A", ta["PRRC2A"]["uoorf"]["peptide"][:24] + "-" + ta["PRRC2A"]["uoorf"]["peptide"][24:], PAL["prrc2a"], PAL["prrc2a_lt"]),
        ("PRRC2B", ta["PRRC2B"]["uoorf"]["peptide"], PAL["prrc2b"], PAL["prrc2b_lt"]),
        ("PRRC2C", ta["PRRC2C"]["uoorf"]["peptide"][:24] + "-" + ta["PRRC2C"]["uoorf"]["peptide"][24:], PAL["prrc2c"], PAL["prrc2c_lt"]),
    ]

    # Plot in 2 compact rows (residues 1-15 on top block, 16-29 on bottom block) so every letter is legible!
    ax.set_xlim(-3.8, 15.6)
    ax.set_ylim(-0.6, 7.3)
    ax.axis("off")

    def col_cons(pos):
        chars = [s[1][pos] for s in seqs if s[1][pos] != "-"]
        if len(set(chars)) == 1 and len(chars) == 3:
            return "all3"
        # Check if PRRC2A & PRRC2B identical
        if seqs[0][1][pos] == seqs[1][1][pos] and seqs[0][1][pos] != "-":
            return "ab"
        if len(set(chars)) == 2:
            return "two"
        return "none"

    blocks = [(0, 15, 6.2), (15, 29, 2.2)]
    for b_start, b_end, y_top in blocks:
        # Residue number ticks at start and end of block
        ax.text(0.45, y_top + 0.58, str(b_start + 1), fontsize=5.0, color=PAL["muted"], ha="center")
        ax.text((b_end - b_start) - 0.55, y_top + 0.58, str(b_end), fontsize=5.0, color=PAL["muted"], ha="center")
        for r_i, (gene, seq, col, col_lt) in enumerate(seqs):
            y = y_top - r_i * 0.92
            ax.text(-3.6, y, gene, fontsize=5.3, fontweight="bold", color=col, ha="left", va="center")
            for k, pos in enumerate(range(b_start, b_end)):
                ch = seq[pos]
                st = col_cons(pos)
                if ch == "-":
                    ax.text(k + 0.45, y, "-", fontsize=5.2, color=PAL["muted"], ha="center", va="center")
                    continue
                if st == "all3":
                    fc, tc = PAL["uoorf"], "white"
                elif st == "ab" and gene in ("PRRC2A", "PRRC2B"):
                    fc, tc = PAL["uoorf_lt"], PAL["ink"]
                else:
                    fc, tc = "white", PAL["muted"]
                ax.add_patch(mpatches.Rectangle((k + 0.06, y - 0.36), 0.82, 0.72, fc=fc, ec="none"))
                ax.text(k + 0.47, y, ch, fontsize=5.0, color=tc,
                        fontweight="bold" if st == "all3" else "normal",
                        ha="center", va="center")
