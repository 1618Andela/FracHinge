"""Graphical abstract for the Engineering Structures submission (min. 1328 x 531 px).
Reads docs/tables/numbers.json (run generate_all.py first)."""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import style
style.apply()
P = style.PALETTE
N = json.load(open(os.path.join(ROOT, "docs", "tables", "numbers.json")))
f = lambda k: float(N[k])
rows = [  # label, mean, CoV %
    ("EN 14651, 14 mixes: blind $f_{R2},f_{R4}$", f("valCblindmean"), f("valCblindcov")),
    ("EN 14651, 77 series: blind $f_{R2},f_{R4}$", f("tiCblindmean"), f("tiCblindcov")),
    ("EN 14651, 48 series: blind $f_L$", f("tifLmean"), f("tifLcov")),
    ("29 RC beams: cracking load", f("boPcmean"), f("boPccov")),
    ("29 RC beams: yielding load", f("boPymean"), f("boPycov")),
    ("12 R/SFRC beams: peak load", f("caBmeanMid"), f("caBcovMid")),
]
fig = plt.figure(figsize=(13.28 / 2.54 * 1.6, 5.31 / 2.54 * 1.6), dpi=300)
# left: flow
axl = fig.add_axes([0.01, 0.05, 0.36, 0.90]); axl.axis("off"); axl.set_xlim(0, 10); axl.set_ylim(0, 10)
boxes = [(0.3, 6.6, "EN 14651\n$f_L,\\ f_{R,1..4}$", P[0]), (3.6, 6.6, "inverse\n$\\sigma(w)$", P[1]),
         (6.9, 6.6, "cohesive\ncracked hinge", P[2]), (3.6, 1.6, "ULS / SLS / $\\rho_{min}$\nECOV, Annex L", P[3])]
for x, y, t, c in boxes:
    axl.add_patch(FancyBboxPatch((x, y), 2.8, 2.4, boxstyle="round,pad=0.08,rounding_size=0.25",
                                 fc="white", ec=c, lw=1.6))
    axl.text(x + 1.4, y + 1.2, t, ha="center", va="center", fontsize=7.5, color="#1f1e1c")
for (x0, y0, x1, y1) in ((3.15, 7.8, 3.55, 7.8), (6.45, 7.8, 6.85, 7.8), (8.3, 6.5, 6.5, 4.1)):
    axl.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color="#52514e", lw=1.0))
axl.text(5.0, 0.4, "fibres + bars + size effect in one model", ha="center", fontsize=7, color="#52514e")
# right: validation summary
axr = fig.add_axes([0.66, 0.17, 0.32, 0.72])
for i, (lab, m, cv) in enumerate(rows[::-1]):
    axr.errorbar(m, i, xerr=m * cv / 100.0, fmt="o", color=P[0], ms=4, capsize=2, lw=1.0)
    axr.text(-0.03, i, lab, ha="right", va="center", fontsize=6.5, transform=axr.get_yaxis_transform())
axr.axvline(1.0, color="#52514e", lw=0.8)
axr.set_yticks([])
axr.set_xlim(0.6, 1.4)
axr.set_xlabel("predicted / experimental (mean $\\pm$ 1 CoV)", fontsize=7)
axr.tick_params(labelsize=7)
axr.set_title("validation, no fitting", fontsize=8)
for spine in ("left", "right", "top"):
    axr.spines[spine].set_visible(False)
out = os.path.join(ROOT, "docs", "submission", "graphical_abstract")
fig.savefig(out + ".png", dpi=300)
fig.savefig(out + ".pdf")
print("written", out)
