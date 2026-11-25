"""Supplementary Video S1: acetate exchange of every cell of the quarter-colony run COL-REF_s50 (manuscript Fig. 4), one frame per full-data snapshot (20 min).
usage: make_movie.py RUNDIR OUT.mp4 [FPS=4]      (needs ffmpeg)
Conventions of make_figure.py (Fig. 4): the quadrant x >= 0, y >= 0 is simulated; radius from the volume as a sphere; acetate exchange in amol/s = flux (mmol/gDW/h) x volume x 0.75 x 1.04 / 3.6,
coolwarm_r centred on 0 (negative = uptake, red; positive = secretion, blue); brown = dead, grey = living cell with a necrosis rate > 0.
The colour range is fixed over the whole video: 1st and 99th percentiles of the living cells of all frames (cells beyond it take the end colour)."""
import sys, os, glob, subprocess, tempfile, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Patch
from matplotlib.colors import TwoSlopeNorm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels; from plot_style import apply_style, FS
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; FPS = int(sys.argv[3]) if len(sys.argv) > 3 else 4
DRY = 0.75 * 1.04; W = 6.0; DPI = 200
def load(s):
    root = ET.parse(s).getroot(); lab = labels(root); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: m[lab.index(k)]
    vol = g("total_volume")
    return dict(t=float(root.find(".//current_time").text), x=g("position[0]"), y=g("position[1]"), vol=vol, dead=g("dead") > 0, nec=g("death_rates[1]") > 0, ac=g("R_EX_ac_e") * vol * DRY / 3.6)
snaps = [load(s) for s in sorted(glob.glob(f"{run}/output*[0-9].xml"))]
live = np.concatenate([d["ac"][~d["dead"] & ~d["nec"]] for d in snaps])
lo, hi = float(np.percentile(live, 1)), float(np.percentile(live, 99))
lo, hi = min(lo, -1e-9), max(hi, 1e-9); norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi); cm = plt.cm.coolwarm_r
last = snaps[-1]; VIEW = float(np.ceil((np.hypot(last["x"], last["y"]).max() + 10) / 50) * 50)
print(f"{len(snaps)} snapshots, t = {snaps[0]['t']:.0f}-{snaps[-1]['t']:.0f} min; colour range {lo:.3g} to {hi:.3g} amol/s; view 0-{VIEW:g} um")
apply_style(W, 4.0)   # labels sized as on a 4 in wide figure
tmp = tempfile.mkdtemp(prefix="movie_S1_")
for i, d in enumerate(snaps):
    fig = plt.figure(figsize=(W, W * 1.12), dpi=DPI); ax = fig.add_axes([0.13, 0.245, 0.80, 0.80 / 1.12])
    rad = np.cbrt(3 * d["vol"] / (4 * np.pi)); keep = (d["x"] <= VIEW) & (d["y"] <= VIEW)
    col = [("saddlebrown" if d["dead"][k] else ("grey" if d["nec"][k] else cm(norm(d["ac"][k])))) for k in np.where(keep)[0]]
    ax.add_collection(PatchCollection([Circle((d["x"][k], d["y"][k]), rad[k]) for k in np.where(keep)[0]], facecolors=col, edgecolors="k", linewidths=0.05))
    ax.set_xlim(0, VIEW); ax.set_ylim(0, VIEW); ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.text(0.97, 0.97, f"t = {d['t'] / 60:.1f} h", transform=ax.transAxes, ha="right", va="top", fontsize=FS["title"] * W / 4.0 * 1.3, fontweight="bold")
    ax.text(0.97, 0.90, f"living cells: {int((~d['dead']).sum()):,} (one quadrant)", transform=ax.transAxes, ha="right", va="top", fontsize=FS["label"] * W / 4.0)
    cax = fig.add_axes([0.20, 0.105, 0.45, 0.022]); cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal")
    cb.set_ticks([v for v in (-0.6, -0.3, 0.0, 0.3, 0.6) if lo <= v <= hi]); cb.set_label("Acetate exchange (amol/s)\n< 0 uptake, > 0 secretion")
    fig.legend([Patch(fc="saddlebrown", ec="k", lw=0.3), Patch(fc="grey", ec="k", lw=0.3)], ["Dead cell", "Starving cell"], loc="lower left", bbox_to_anchor=(0.70, 0.07), frameon=False)
    fig.savefig(f"{tmp}/f{i:04d}.png", facecolor="white"); plt.close(fig)
subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", f"{tmp}/f%04d.png", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out], check=True)
print(len(snaps), "frames ->", out)
