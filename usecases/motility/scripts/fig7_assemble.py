"""Assemble Fig. 7: panels A-F (results/fig7/fig7_AF.png, fig7_composite.py) above panel G (results/fig7/panel_g.png, fig7_panel_g.py), with the panel letters. usage: fig7_assemble.py [OUT_PREFIX=results/fig7/fig7]  (run from the use-case directory)"""
import sys, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt, matplotlib.image as mpimg
out = sys.argv[1] if len(sys.argv) > 1 else "results/fig7/fig7"
a, g = mpimg.imread("results/fig7/fig7_AF.png"), mpimg.imread("results/fig7/panel_g.png"); Ha, Wa = a.shape[:2]; Hg, Wg = g.shape[:2]
gw = 0.36; gh = gw * Wa * Hg / Wg / Ha            # panel G: 36% of the width, in proportion
fig = plt.figure(figsize=(9, 9 * (Ha / Wa) * (1 + gh + 0.04)), dpi=300); T = 1 + gh + 0.04
ax = fig.add_axes([0, 1 - 1 / T, 1, 1 / T]); ax.imshow(a); ax.axis("off")
ag = fig.add_axes([0.03, 0.0, gw, gh / T]); ag.imshow(g); ag.axis("off")
for t, x, ytop in (("A", 0.005, 0.005), ("B", 0.385, 0.005), ("C", 0.005, 0.232), ("D", 0.385, 0.232), ("E", 0.005, 0.447), ("F", 0.385, 0.447)):   # ytop = distance of the letter from the top of the figure (fraction)
    fig.text(x, 1 - ytop, t, fontsize=16, fontweight="bold", va="top")
fig.text(0.005, gh / T, "G", fontsize=16, fontweight="bold", va="top")
fig.savefig(out + ".png", facecolor="white"); fig.savefig(out + ".pdf", facecolor="white"); print("wrote", out)
