"""Figure 1: Genomic architecture, 600-My deuterostome origin, and database blind spot."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=118, style=STYLE)
fig.panel("panels/f1_locus.py", x=4, y=2, w=175, h=42)
fig.panel("panels/f1_origin.py", x=4, y=49, w=114, h=66)
fig.panel("panels/f1_db_gap.py", x=122, y=49, w=57, h=66)
