"""Shared Nature journal styling and matte semantic palette for alt-PRRC2A/B preprint.

Pure white backgrounds, crisp black typography (Arial), zero gray callout boxes,
and restrained color encoding used exclusively for data semantics.
"""
from pathlib import Path
from fignimble import Labels, Style

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR.parent.parent / "data"

PAL = {
    "prrc2a": "#C85A3F",      # Matte terracotta: PRRC2A / alt-PRRC2A
    "prrc2a_lt": "#F6DFD9",
    "prrc2b": "#22859A",      # Matte petrol teal: PRRC2B / alt-PRRC2B
    "prrc2b_lt": "#D7EFF4",
    "frame0": "#4A6FA5",      # Matte slate blue: Frame-0 PRRC2 / BAT2 / non-STOP baseline
    "frame0_lt": "#DCE6F2",
    "uoorf": "#C9933B",       # Matte ochre: Conserved +2 uoORF / reinitiation circuit
    "uoorf_lt": "#F7EBD5",
    "motif": "#4E8A6E",       # Matte sage green: Conserved motif / CBE switches / canonical splicing
    "motif_lt": "#D9ECE3",
    "prrc2c": "#7E848C",      # Matte warm gray: Pseudogenized PRRC2C negative control
    "prrc2c_lt": "#EAECEF",
    "clinvar": "#8E4B78",     # Matte dusty plum: ClinVar variants
    "ink": "#111111",         # Pure black text & axes
    "muted": "#5C6168",       # Subtle secondary text
    "grid": "#E6E8EC",        # Hairline grid
}

STYLE = Style.journal(
    "nature",
    font="Arial",
    size=7.0,
    linewidth=0.65,
    labels=Labels.journal("nature"),
    palette=[
        PAL["prrc2a"],
        PAL["prrc2b"],
        PAL["frame0"],
        PAL["uoorf"],
        PAL["motif"],
        PAL["prrc2c"],
    ],
    despine=True,
    rc={
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    },
)
