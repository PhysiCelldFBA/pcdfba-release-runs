"""Fig. 7 panel G: time at which the founder cell first passes x = 25 um against Km(glucose), for the two supplies (founder series of the Km scans: results/mot4_series.csv = 1.0 mM, results/mot3_series.csv = 0.325 mM; Km 0.02 = the reference runs). Open symbols = the founder did not reach x = 25 um in 10 h. usage: fig7_panel_g.py [OUT_PREFIX=results/fig7_panel_g]"""
import sys, re, pandas as pd, numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
out = sys.argv[1] if len(sys.argv) > 1 else "results/fig7_panel_g"; X = 25.0
def scan(csv, pat, ref):
    d = pd.read_csv(csv); rows = []
    for run, g in d.groupby("run"):
        m = re.match(pat, run); km = float(m.group(1)) if m else (0.02 if run == ref else None)
        s = re.match(r"MOT3-SENS-LOW_Km([\d.]+)_B1_s50", run)   # Km 0.5 / 2 x 0.02 mM at the reference supply
        if s: km = round(0.02 * float(s.group(1)), 6)
        if km is None: continue
        hit = g[g.x_founder >= X]; rows.append((km, hit.t_min.iloc[0] if len(hit) else np.nan))
    return pd.DataFrame(rows, columns=["Km", "t"]).sort_values("Km")
S = {"HIGH": scan("results/mot4_series.csv", r"MOT4-KM-HIGH_Km([\d.]+)_s50", "MOT4-HIGH_s50"), "LOW": scan("results/mot3_series.csv", r"MOT3-KM-LOW_Km([\d.]+)_s50", "MOT3-LOW_s50")}
fig, ax = plt.subplots(figsize=(3.6, 2.8)); top = max(s.t.max() for s in S.values()) * 1.25
for sup, col, lab, dy in (("HIGH", "#1b7837", "1.0 mM", 1.0), ("LOW", "#b8860b", "0.325 mM", 0.93)):
    s = S[sup]; m = s[s.t.notna()]; n = s[s.t.isna()]
    ax.plot(m.Km, m.t, "-o", color=col, label=f"{lab} at the right wall", ms=4, lw=1.3); ax.plot(n.Km, np.full(len(n), top * dy), "o", mfc="none", color=col, ms=5)
    print(sup, s.to_dict("list"))
ax.set_xscale("log"); ax.set_ylim(0, top * 1.08); ax.set_xlabel("$K_m$ of glucose (mM)"); ax.set_ylabel("Time to reach x = 25 µm (min)")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=7, loc="center left"); fig.tight_layout()
fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".svg"); print("wrote", out)
