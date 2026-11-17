"""Supplementary Video S4: time-lapse of the consortium run CF3-NC50_s50 (manuscript Fig. 6) up to 45 h.
usage: make_movie.py RUNDIR OUT.mp4 [T_END_MIN=2700] [FPS=8]      (needs ffmpeg)
One frame per snapshot (30 min). Cells are drawn at their position with 1.3 times the radius of a sphere of their volume (so that neighbours touch, as in the previous version); colour = growth rate (1/h),
Blues (C. beijerinckii) and Oranges (M. barkeri), colour range 0-0.15 as in the previous version of the video; counters of living cells and of
apoptotic and necrotic cells. The run is cut at 45 h, the last time with no cell beyond the edge of the mesh (Supplementary Methods S7.5)."""
import sys, os, glob, subprocess, tempfile, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; TEND = float(sys.argv[3]) if len(sys.argv) > 3 else 2700.0; FPS = int(sys.argv[4]) if len(sys.argv) > 4 else 8
files = sorted(glob.glob(f"{run}/output*[0-9].xml"))
def tof(f): return float(ET.parse(f).getroot().find(".//current_time").text)
files = [f for f in files if tof(f) <= TEND + 1e-6]
VMAX, L = 0.15, 100.0; cmaps = {0: plt.get_cmap("Blues_r"), 1: plt.get_cmap("Oranges_r")}; names = {0: "C. beijerinckii", 1: "M. barkeri"}
tmp = tempfile.mkdtemp(prefix="movie_")
for i, f in enumerate(files):
    root = ET.parse(f).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
    c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: c[lab.index(k)]
    x, y, v, ty, gr, dead, dm = g("position[0]"), g("position[1]"), g("total_volume"), g("cell_type"), g("growth_rate"), g("dead") > 0, g("current_death_model")
    fig = plt.figure(figsize=(8, 8), dpi=100, facecolor="black"); ax = fig.add_axes([0.16, 0.20, 0.70, 0.70], facecolor="black")
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    pt_per_um = 0.70 * 8 * 72 / (2 * L); rad = 1.3 * (3 * v / (4 * np.pi)) ** (1 / 3)
    for k in (0, 1):
        m = (ty == k) & ~dead
        ax.scatter(x[m], y[m], s=(2 * rad[m] * pt_per_um) ** 2, c=np.clip(gr[m], 0, VMAX), cmap=cmaps[k], vmin=0, vmax=VMAX, linewidths=0, marker="o")
    d = dead; ax.scatter(x[d & (dm == 0)], y[d & (dm == 0)], s=(2 * rad[d & (dm == 0)] * pt_per_um) ** 2, c="#ff3b30", linewidths=0); ax.scatter(x[d & (dm > 0)], y[d & (dm > 0)], s=(2 * rad[d & (dm > 0)] * pt_per_um) ** 2, c="#a855f7", linewidths=0)
    for sp in ax.spines.values(): sp.set_color("white")
    ax.tick_params(colors="white"); ax.set_xlabel("x (µm)", color="white", fontsize=13, fontweight="bold"); ax.set_ylabel("y (µm)", color="white", fontsize=13, fontweight="bold")
    ax.text(0.02, 0.97, f"t = {t / 60:.1f} h", transform=ax.transAxes, color="white", fontsize=16, fontweight="bold", va="top")
    n = [int(((ty == k) & ~dead).sum()) for k in (0, 1)]
    ax.text(0.98, 0.97, f"C. beijerinckii: {n[0]}   M. barkeri: {n[1]}\nApoptotic: {int((dead & (dm == 0)).sum())}   Necrotic: {int((dead & (dm > 0)).sum())}", transform=ax.transAxes, color="white", fontsize=11, fontweight="bold", va="top", ha="right")
    for k, x0 in ((0, 0.10), (1, 0.55)):
        cax = fig.add_axes([x0, 0.08, 0.30, 0.025]); sm = plt.cm.ScalarMappable(cmap=cmaps[k], norm=plt.Normalize(0, VMAX)); cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.ax.tick_params(colors="white", labelsize=8); cb.outline.set_edgecolor("white"); cb.set_label(f"{names[k]} — growth rate (1/h)", color="white", fontsize=9, fontweight="bold")
    fig.savefig(f"{tmp}/f{i:04d}.png", facecolor="black"); plt.close(fig)
subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", f"{tmp}/f%04d.png", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out], check=True)
print(len(files), "frames to t =", t / 60, "h ->", out)
