"""Motility runs: per run the position of the founder cell and of the cells, speed, metabolic mode, first division and cell number. usage: analyse_mot.py [RUNS_DIR=runs] [GLOB=MOT*-*] [OUT_PREFIX=mot]  (e.g. analyse_mot.py runs "MOT2-*" mot2)
Writes results/mot_series.csv (per run and snapshot: time, cells, x of the founder cell (smallest ID), mean x, share of cells in ATP mode (metabolic_mode = 1), founder speed over the interval) and results/mot_summary.csv (per run: x at 10 h, mean speed, share of time in ATP mode, time of the first division, cells at 10 h)."""
import sys, os, glob, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); rd = sys.argv[1] if len(sys.argv) > 1 else f"{UC}/runs"; pat = sys.argv[2] if len(sys.argv) > 2 else "MOT*-*"; pre = sys.argv[3] if len(sys.argv) > 3 else "mot"; ser, summ = [], []
for d in sorted(glob.glob(f"{rd}/{pat}")):
    run = os.path.basename(d); prev = None; fd = None
    for xf in sorted(glob.glob(f"{d}/output[0-9]*[0-9].xml")):
        root = ET.parse(xf).getroot(); lab = labels(root); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; col = lambda k: m[lab.index(k)]; t = float(root.find(".//current_time").text)
        ids = np.asarray(col("ID")); j = int(np.argmin(ids)); x = float(col("position[0]")[j]); n = m.shape[1]; mode = float(np.mean(col("metabolic_mode")))
        sp = abs(x - prev[1]) / (t - prev[0]) if prev and t > prev[0] else np.nan; prev = (t, x)
        if fd is None and n > 1: fd = t
        ser.append(dict(run=run, t_min=t, cells=n, x_founder=x, x_mean=float(np.mean(col("position[0]"))), atp_mode_share=mode, speed_um_min=sp))
    S = pd.DataFrame([r for r in ser if r["run"] == run]); last = S[S.t_min == S.t_min.max()].iloc[0]
    summ.append(dict(run=run, t_end_min=last.t_min, x_founder_end=last.x_founder, x_mean_end=last.x_mean, cells_end=int(last.cells), first_division_min=fd, mean_speed_um_min=S.speed_um_min.mean(), max_speed_um_min=S.speed_um_min.max(), atp_mode_share_time=S.atp_mode_share.mean()))
os.makedirs(f"{UC}/results", exist_ok=True); pd.DataFrame(ser).to_csv(f"{UC}/results/{pre}_series.csv", index=False); D = pd.DataFrame(summ); D.to_csv(f"{UC}/results/{pre}_summary.csv", index=False)
pd.set_option("display.width", 220); pd.set_option("display.float_format", lambda x: f"{x:.4g}"); print(D.to_string(index=False))
