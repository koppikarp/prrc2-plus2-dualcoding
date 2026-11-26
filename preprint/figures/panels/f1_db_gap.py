"""Figure 1, Panel c: Database annotation & homology blind spot for alt-PRRC2A/B."""
import numpy as np
from style import PAL


def draw(fig):
    ax = fig.add_subplot(111)

    labels = [
        "100-way vertebrates\n(this study)",
        "UniRef100 hits\n(MMseqs2)",
        "Ensembl Compara\n(species)",
        "HGNC / RefSeq /\nNCBI / GTEx",
        "Pfam / InterPro\ndomains",
        "Published papers\n/ preprints",
    ]
    vals = [94, 7, 2, 0, 0, 0]
    cols = [PAL["prrc2b"], PAL["frame0"], PAL["uoorf"], PAL["prrc2c"], PAL["prrc2c"], PAL["prrc2c"]]

    ys = np.arange(len(labels))[::-1]
    bars = ax.barh(ys, vals, height=0.56, color=cols, edgecolor="none")

    for y, v in zip(ys, vals):
        ax.text(v + 2.2, y, str(v), va="center", ha="left", fontsize=5.8,
                fontweight="bold" if v in (94, 0) else "normal",
                color=PAL["prrc2b"] if v == 94 else PAL["ink"])

    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=5.4)
    ax.set_xlim(0, 114)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Annotated entries / species", fontsize=6.0)
    ax.tick_params(axis="x", labelsize=5.5)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
