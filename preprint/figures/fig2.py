"""Figure 2: Orthogonal vertebrate signatures of purifying selection in Frame +2."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=118, style=STYLE)
fig.panel("panels/f2_sliding_4d.py", x=4, y=2, w=95, h=54)
fig.panel("panels/f2_paired_4d.py", x=104, y=2, w=75, h=54)
fig.panel("panels/f2_five_frame_stops.py", x=4, y=61, w=110, h=54)
fig.panel("panels/f2_aa_pres.py", x=119, y=61, w=60, h=54)
