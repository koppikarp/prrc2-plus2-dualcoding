"""Graphical Abstract: Clean three-stage visual overview (minimal text, zero gray boxes)."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=52, style=STYLE)
fig.panel("panels/ga_overview.py", x=1, y=1, w=181, h=50, label="")
