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
