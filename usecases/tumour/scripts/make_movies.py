"""Supplementary Videos S2 (growth rate) and S3 (uptake stack) of the tumour run TUM5-PFBA_s50 (manuscript Fig. 5, Supplementary Figs. S8 and S9), one frame per full-data snapshot (1 h).
usage: make_movies.py RUNDIR OUTDIR [FPS=8]      (needs ffmpeg; writes OUTDIR/Video_S2.mp4 and OUTDIR/Video_S3.mp4)
S2 is panel A of Fig. 5 (make_figure.py: same style, colour scale 0-0.026803 1/h, Greens) at every hour.
S3 is Fig. S9 (make_supp_figs.py: same rows, dry fraction of the snapshot, O2 row multiplied by the cell volume) at every hour, with the default style of plot_cancer_tissue_spatial.py.
Colour scales of S3 are fixed over the video: 0 (or the minimum) to the 98th percentile of the living cells of all frames, as _living_value_extent does for one frame.
In both, the dashed lines are the proliferative edge and the living edge of each snapshot (95th percentiles of x, as in make_figure.py)."""
import sys, os, glob, io, types, subprocess, tempfile, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels; from plot_style import apply_style, FS
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont
sys.modules["pctk"] = types.ModuleType("pctk"); sys.modules["pctk"].multicellds = types.SimpleNamespace()   # not needed: no reader is used
import figure_paths; figure_paths.ensure_paths = lambda *a, **k: None
import plot_cancer_tissue_spatial as P
from plot_cancer_tissue_spatial_gif import _render_uptake_stack, _pad_to_common_size
from plot_cells_spatial import plot_cells_spatial
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; FPS = int(sys.argv[3]) if len(sys.argv) > 3 else 8; os.makedirs(out, exist_ok=True)

def load(xf):
    r = ET.parse(xf).getroot(); lab = labels(r); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; t = float(r.find(".//current_time").text) / 60
    df = pd.DataFrame({k: m[lab.index(k)] for k in lab if k.startswith("R_EX") and not k.endswith("_lb")})
    for k, n in (("position[0]", "x_position"), ("position[1]", "y_position"), ("position[2]", "z_position"), ("total_volume", "total_volume"), ("growth_rate", "growth_rate"), ("dead", "dead"), ("death_rates[1]", "death_rates_1"), ("fluid_fraction", "fluid_fraction")): df[n] = m[lab.index(k)]
    df["time"] = t; alive = df.dead == 0; grow = alive & (df.growth_rate > 0.5 * 0.0137)
    pe = float(np.percentile(df.x_position[grow], 95)) if grow.any() else 0.0; le = float(np.percentile(df.x_position[alive], 95))
    return t, df, (pe, le), float(1 - df.fluid_fraction[alive].mean())

snaps = [load(f) for f in sorted(glob.glob(f"{run}/output*[0-9].xml"))]
print(f"{len(snaps)} snapshots, t = {snaps[0][0]:.0f}-{snaps[-1][0]:.0f} h")
font = ImageFont.truetype(font_manager.findfont("DejaVu Sans:bold"), 40)

def write(frames, name):
    frames = _pad_to_common_size(frames); tmp = tempfile.mkdtemp(prefix=f"movie_{name}_")
    for i, f in enumerate(frames): Image.fromarray(f).convert("RGB").save(f"{tmp}/f{i:04d}.png")
    subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", f"{tmp}/f%04d.png", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", f"{out}/{name}.mp4"], check=True)
    print(len(frames), "frames ->", f"{out}/{name}.mp4")

# ---- S2: Fig. 5A at every hour ----
apply_style(6.3)
s2 = []
for t, df, zx, dry in snaps:
    sty = dict(X_MIN=0.0, X_MAX=320.0, Y_MIN=-55.0, Y_MAX=55.0, skip_dead=False, vmax=0.026803, vmin=0.0, vessel_style="bar", vessel_label=True, show_zones=True, zone_x=zx, zone_labels=True, zone_bands=False,
               show_nuclei=True, nucleus_scale=1.0, cell_shape="blob", roughness=0.05, shape_seed=0, perspective=0.1, depth_fade=0.05, show_x_ruler=True, show_scale_bar=False, x_major=50.0, x_minor=10.0,
               xlabel="Distance from vessel (µm)", fontsize=FS["tick"], annotation_fontsize=FS["label"], cbar_tick_fontsize=FS["tick"], figsize=(6.3, 2.2))
    fa = plot_cells_spatial(df, time_point=t, color_config={"value_fn": lambda x: x.growth_rate, "cmap": plt.cm.Greens, "label": "Growth rate (1/h)"}, **sty)
    fa.text(0.995, 0.995, f"t = {t:.0f} h", ha="right", va="top", fontsize=FS["title"], fontweight="bold", transform=fa.transFigure)
    b = io.BytesIO(); fa.savefig(b, dpi=250, bbox_inches="tight", pad_inches=0.08, format="png", facecolor="white"); plt.close(fa)
    s2.append(np.array(Image.open(b).convert("RGBA")))
write(s2, "Video_S2")

# ---- S3: Fig. S9 at every hour ----
plt.rcdefaults(); matplotlib.use("Agg")
P._o2_uptake = lambda row: -float(row.R_EX_o2_e) * float(row.total_volume) * P._RHO * P._SOLID / P._MIN_PER_H
specs = {s["stem"]: s for s in P._figure_specs()}; specs["oxygen_uptake"]["required"] = ["R_EX_o2_e", "total_volume"]
specs = [dict(specs[s]) for s in P._DEFAULT_STACK_STEMS]
for sp in specs:   # fixed scale over all frames (living cells without a necrosis rate, as _living_value_extent)
    vals = []
    for t, df, zx, dry in snaps:
        P._SOLID = dry; ok = (df.dead < 1) & ~(df.death_rates_1 > 0)
        vals += [sp["color_config"]["value_fn"](x) for x in df[ok].itertuples()]
    v = np.asarray(vals, float); v = v[np.isfinite(v)]; lo = min(0.0, float(v.min())); hi = float(np.quantile(v, 0.98))
    sp["vmin"], sp["vmax"] = lo, hi; print(f"  {sp['stem']}: {lo:.4g} to {hi:.4g}")
s3 = []
for t, df, zx, dry in snaps:
    P._SOLID = dry; sty = dict(P._DEFAULT_STYLE); sty["zone_x"] = zx
    fr = _render_uptake_stack(df, [t], specs, sty, show_colorbar=True)[0]
    im = Image.fromarray(fr).convert("RGBA"); ImageDraw.Draw(im).text((im.width - 30, 15), f"t = {t:.0f} h", fill="black", font=font, anchor="ra"); s3.append(np.array(im))
    plt.close("all")
write(s3, "Video_S3")
