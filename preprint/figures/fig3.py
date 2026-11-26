"""Figure 3: Active purifying selection across 807,162 human genomes (gnomAD v4) and ClinVar mapping."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=120, style=STYLE)
fig.panel("panels/f3_f0syn_strata.py", x=4, y=2, w=88, h=55)
fig.panel("panels/f3_f0mis_gradient.py", x=96, y=2, w=83, h=55)
fig.panel("panels/f3_forest_cmh.py", x=4, y=62, w=82, h=55)
fig.panel("panels/f3_motif_map.py", x=90, y=62, w=89, h=55)
