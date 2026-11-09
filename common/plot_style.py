"""One colour per substrate for every figure of every use case (import: sys.path.insert(0, <repo>/common); from plot_style import SC, substrate_cmap).
Colours are chosen to be told apart in the same panel: glucose gold, oxygen blue, acetate red, CO2 purple, H2 magenta-pink, methane green, biomass black.
substrate_cmap(name) returns the matching sequential colormap (white -> substrate colour -> dark) for density maps."""
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
SC = dict(glucose="#f1c232", oxygen="#3b82f6", acetate="#c00000", co2="#8e44ad", h2="#6cc3f0", methane="#18a081", biomass="#000000")
_CMAP = {"oxygen": "Blues", "acetate": "Reds", "co2": "Purples"}
_DARK = {"glucose": "#f1c232", "h2": "#12527f", "methane": "#0b4f3d"}   # end of the scale
def substrate_cmap(name):
    if name in _CMAP: return plt.get_cmap(_CMAP[name])
    return LinearSegmentedColormap.from_list(name, ["#ffffff", SC[name], _DARK[name]] if name != "glucose" else ["#ffffff", SC[name]])

# ---- type sizes shared by all the figures -------------------------------------------------------------------------
# Sizes are given in points AS PRINTED. A figure is drawn at any size: apply_style(fig_width_in, print_width_in) scales the sizes by
# fig_width / print_width, so that after the figure is scaled to its width in the manuscript (text width 6.3 in: A4, 2.5 cm margins;
# 0.5 linewidth = 3.15 in) every label, tick label, legend and panel letter has the same size on the page.
import matplotlib as mpl
PRINT_W = 6.3
FS = dict(label=8.0, tick=7.0, legend=7.0, title=8.0, annot=6.5, letter=12.0, lw=1.2, axw=0.6)
def apply_style(fig_w_in, print_w_in=PRINT_W):
    k = fig_w_in / print_w_in
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": FS["label"] * k,
        "axes.labelsize": FS["label"] * k, "axes.titlesize": FS["title"] * k, "xtick.labelsize": FS["tick"] * k, "ytick.labelsize": FS["tick"] * k,
        "legend.fontsize": FS["legend"] * k, "axes.linewidth": FS["axw"] * k, "xtick.major.width": FS["axw"] * k, "ytick.major.width": FS["axw"] * k,
        "xtick.major.size": 3 * k, "ytick.major.size": 3 * k, "lines.linewidth": FS["lw"] * k, "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none"})
    return k
def fs(name, fig_w_in, print_w_in=PRINT_W): return FS[name] * fig_w_in / print_w_in
