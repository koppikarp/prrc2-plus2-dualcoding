"""Figure 4: Conserved uoORF reinitiation circuit, Kozak lock, and constitutive splicing."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=112, style=STYLE)
fig.panel("panels/f4_circuit.py", x=4, y=2, w=175, h=46)
fig.panel("panels/f4_kozak.py", x=4, y=53, w=58, h=56)
fig.panel("panels/f4_splicing.py", x=66, y=53, w=56, h=56)
fig.panel("panels/f4_uoorf_align.py", x=126, y=53, w=53, h=56)
