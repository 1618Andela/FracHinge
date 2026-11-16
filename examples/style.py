"""Common plotting style (validated categorical palette, secondary encoding by line style)."""
import matplotlib as mpl

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (1, 1)), (0, (8, 2))]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def apply():
    mpl.rcParams.update({
        "figure.figsize": (6.0, 3.8),
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "font.family": "serif",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#52514e",
        "axes.labelcolor": "#0b0b0b",
        "xtick.color": "#52514e",
        "ytick.color": "#52514e",
        "axes.grid": True,
        "grid.color": "#e4e3df",
        "grid.linewidth": 0.6,
        "lines.linewidth": 1.6,
        "lines.markersize": 4.5,
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
    })
