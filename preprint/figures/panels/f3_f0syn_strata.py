"""Figure 3, Panel a: gnomAD v4 Frame-0-synonymous SNVs observed (%) stratified by mutational class."""
import json
import matplotlib.patches as mpatches
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    gn = json.loads((DATA_DIR / "master_preprint_data.json").read_text())["gnomad_v4"]
    ab = gn["POOLED_AB"]["tests"]
    c2 = gn["PRRC2C_NEGCTRL"]["tests"]

    items = [
        ("CpG Ti\n(A+B)", ab["CpG_Ti__F0syn_p2STOP_vs_nonSTOP"]),
        ("Non-CpG Ti\n(A+B)", ab["nonCpG_Ti__F0syn_p2STOP_vs_nonSTOP"]),
        ("Tv\n(A+B)", ab["Tv__F0syn_p2STOP_vs_nonSTOP"]),
        ("Non-CpG\n(A+B)", ab["nonCpG_all__F0syn_p2STOP_vs_nonSTOP"]),
        ("Non-CpG\n(2C ctrl)", c2["nonCpG_all__F0syn_p2STOP_vs_nonSTOP"]),
    ]

    x = np.arange(len(items))
    w = 0.34
    ns_vals = [it[1]["nonstop"]["pct"] for it in items]
    st_vals = [it[1]["stop"]["pct"] for it in items]
    st_lbls = [f"{it[1]['stop']['obs']}/{it[1]['stop']['tot']}" for it in items]

    ax.axvline(3.5, ls=":", lw=0.75, color=PAL["muted"], zorder=1)

    ax.bar(x - w / 2, ns_vals, width=w, color=PAL["frame0"], zorder=3)
    st_cols = [PAL["prrc2a"]] * 4 + [PAL["prrc2c"]]
    b2 = ax.bar(x + w / 2, st_vals, width=w, color=st_cols, zorder=3)

    for rect, v, lbl in zip(b2, st_vals, st_lbls):
        ax.text(rect.get_x() + rect.get_width() / 2, v + 2.0, lbl,
                ha="center", va="bottom", fontsize=5.0,
                fontweight="bold" if lbl.startswith("0/") else "normal",
                color=PAL["prrc2a"] if lbl.startswith("0/") else PAL["ink"])

    handles = [
        mpatches.Patch(fc=PAL["frame0"], label="+2 non-STOP"),
        mpatches.Patch(fc=PAL["prrc2a"], label="+2 STOP (A+B)"),
        mpatches.Patch(fc=PAL["prrc2c"], label="+2 STOP (2C)"),
    ]
    ax.set_xticks(x)
    ax.set_xticklabels([it[0] for it in items], fontsize=5.2)
    ax.set_ylim(0, 108)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("F0-syn. sites observed (%)", fontsize=6.1)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.grid(axis="y", color=PAL["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=handles, frameon=False, fontsize=5.1, loc="upper center",
              bbox_to_anchor=(0.48, 0.99), handlelength=1.1, labelspacing=0.22)
