"""Benchmark of PhysiCelldFBA against COMETS 2.12.5 and BacArena 1.8.2 on ACE-REF and the consortium of Fig. 6.

BacArena (2026-11-19): ACE v6 (medium and Km scaled to the PhysiCelldFBA box, ATPM, biomass and Vmax converted to the step,
one cell per site on a 17 x 17 lattice, death rate 0, step 0.01 h; data/aceref_bacarena); the Fig. 6 BacArena run
restores the SBML exchange bounds after Bac() and is the closed arena without initial H2 / acetate pool. COMETS: the
Fig. 6 run with Convection 2D. Only the runs shown in the two figures are read.

Writes results/benchmark_summary.csv (every number used in the text), results/benchmark_fig6_criteria.csv (criteria in
results/cf_fig6_qualitative_criteria.md) and the figures results/supp_fig_benchmark_{aceref,crossfeeding}.{png,pdf}.

Run from benchmark_frameworks/:  python scripts/make_benchmark.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
UC = HERE.parent / "usecases"
D = HERE / "data"
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

TOOLS = {"PhysiCelldFBA": dict(color="#1b6ca8", ls="-"),
         "COMETS": dict(color="#d1495b", ls="--"),
         "BacArena": dict(color="#3c8d2f", ls=":")}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "lines.linewidth": 1.4})
rows = []


def rec(case, tool, quantity, value, unit, note=""):
    rows.append(dict(case=case, tool=tool, quantity=quantity, value=value, unit=unit, note=note))


def secs(a, b):
    return (pd.Timestamp(b) - pd.Timestamp(a)).total_seconds()


# ----------------------------------------------------------------------------- ACE-REF
# BacArena ACE (2026-11-19, user decision): v6 with one cell per site (17 x 17 lattice of 2.1 um, bacnum 1), death rate 0
# (PhysiCelldFBA cells do not lose mass when they cannot grow) and a step of 0.01 h. Only the runs of the figures are used.
V6 = "aceref_bacarena/bacarena_aceref_v6_"
ACE_BA = {"BacArena": (V6 + "dt0.01_n17_d0_timeseries.csv", V6 + "dt0.01_n17_d0_timing.json", "mpws")}
ACE_VOL_L = 36e-6 ** 3 * 1e3   # PhysiCelldFBA ACE box, 36 um cube, in L (biomass_gL x volume = dry biomass)


def load_ace():
    p = pd.read_csv(UC / "acetate/results/ace_series.csv")
    p = p[p.run == "ACE-REF_s50"].sort_values("t")
    p = pd.DataFrame(dict(t_h=p.t / 60, fold=p.biomass_gL / p.biomass_gL.iloc[0], pg=p.biomass_gL * ACE_VOL_L * 1e12, n=p.n_alive,
                          glucose=p.glucose_mean, oxygen=p.oxygen_mean, acetate=p.acetate_mean))
    c = pd.read_csv(D / "aceref_comets/comets_aceref_v2_timeseries.csv").sort_values("time_min")
    c = pd.DataFrame(dict(t_h=c.time_min / 60, fold=c.biomass_g / c.biomass_g.iloc[0], pg=c.biomass_g * 1e12,
                          glucose=c.glucose_mM.fillna(0), oxygen=c.oxygen_mM.fillna(0), acetate=c.acetate_mM.fillna(0)))
    def ba(f):
        b = pd.read_csv(D / f).sort_values("time_h")
        b0 = float(b.biomass_pg.iloc[0])
        b = pd.DataFrame(dict(t_h=b.time_h, fold=b.biomass_pg / b0, pg=b.biomass_pg, n=b.n_cells,
                              glucose=b.glucose_mM, oxygen=b.oxygen_mM, acetate=b.acetate_mM))
        b.attrs.update(b0=b0, b0_note="BacArena default cell weight x 20 founders (PhysiCelldFBA / COMETS 6.76 pg)")
        return b
    return {"PhysiCelldFBA": p, "COMETS": c, **{k: ba(v[0]) for k, v in ACE_BA.items()}}


def summarise_ace(S):
    prov = pd.read_csv(UC / "acetate/provenance.csv").set_index("run_id").loc["ACE-REF_s50"]
    timing = {"PhysiCelldFBA": (secs(prov.start, prov.end), int(prov.threads), prov.machine)}
    j = json.load(open(D / "aceref_comets/comets_aceref_v2_timing.json")); timing["COMETS"] = (j["elapsed_s"], 1, "bsclife006")
    for tool, (_, f, host) in ACE_BA.items():
        j = json.load(open(D / f)); timing[tool] = (j["elapsed_s"], 1, host)
    for tool, s in S.items():
        rec("ACE", tool, "biomass fold change 0-10 h", s.fold.iloc[-1], "x")
        rec("ACE", tool, "biomass at 0 / max / 10 h", f"{s.pg.iloc[0]:.3g} / {s.pg.max():.3g} / {s.pg.iloc[-1]:.3g}", "pg")
        if "b0" in s.attrs:
            rec("ACE", tool, "initial biomass", s.attrs["b0"], "pg", s.attrs.get("b0_note", ""))
        i = s.acetate.idxmax()
        rec("ACE", tool, "peak acetate", s.acetate.max(), "mM")
        rec("ACE", tool, "time of peak acetate", s.t_h[i] * 60, "min")
        rec("ACE", tool, "acetate at 10 h", s.acetate.iloc[-1], "mM")
        gl = s[s.glucose < 0.01]
        rec("ACE", tool, "glucose below 0.01 mM from", gl.t_h.iloc[0] * 60 if len(gl) else np.nan, "min",
            "" if len(gl) else "not reached by 10 h")
        rec("ACE", tool, "glucose at 10 h", s.glucose.iloc[-1], "mM")
        if "n" in s:
            rec("ACE", tool, "cells or individuals at 0 / 10 h", f"{int(s.n.iloc[0])} / {int(s.n.iloc[-1])}", "")
        t, th, host = timing[tool]
        rec("ACE", tool, "wall time", t, "s", f"{th} thread(s), {host}")


INK_ = "#222222"
# Figures (2026-11-19, simplified): tools by line style, metabolites / species by colour, one legend per figure
TOOL_LS = {"PhysiCelldFBA": "-", "COMETS": "--", "BacArena": ":"}
# biomass line per tool; PhysiCelldFBA thick and light so that COMETS, which lies on it, stays visible on top
TOOL_BIO = {"PhysiCelldFBA": dict(color="0.72", ls="-", lw=3.2, zorder=1), "COMETS": dict(color=INK_, ls="--", lw=1.3, zorder=3),
            "BacArena": dict(color=INK_, ls=":", lw=1.6, zorder=2)}
MET = {"glucose": ("Glucose", "#2a78d6"), "oxygen": ("O$_2$", "#eb6834"), "acetate": ("Acetate", "#1baf7a")}
INK = INK_


def _save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _ace_axes(a):
    a.set_xlim(0, 10); a.set_xlabel("time (h)")
    a.axvline(3.5, color="0.75", lw=0.6, zorder=0)


def _legend_handles(tools=True, mets=True):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], label=t, **{k: v for k, v in TOOL_BIO[t].items() if k != "zorder"}) for t in TOOL_BIO] if tools else []
    h += [Line2D([], [], color=c, lw=2, label=l) for l, c in MET.values()] if mets else []
    return h


def fig_ace(S):
    """Biomass (pg) of the three tools (left) | one concentration panel per tool, stacked on a shared time axis (right);
    one legend outside the panels."""
    tools = ["PhysiCelldFBA", "COMETS", "BacArena"]
    fig = plt.figure(figsize=(7.0, 3.6))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.15, 1], wspace=0.32, hspace=0.18)
    a0 = fig.add_subplot(gs[:, 0])
    for t in tools:
        a0.plot(S[t].t_h, S[t].pg, **TOOL_BIO[t])
    a0.set_ylabel("dry biomass (pg)"); a0.set_ylim(bottom=0)
    _ace_axes(a0)
    a0.text(-0.2, 1.02, "A", transform=a0.transAxes, fontweight="bold", fontsize=10)
    axs = [fig.add_subplot(gs[k, 1]) for k in range(3)]
    for a in axs[1:]:
        a.sharex(axs[0]); a.sharey(axs[0])
    for a, t in zip(axs, tools):
        for m, (_, c) in MET.items():
            a.plot(S[t].t_h, S[t][m], color=c, ls="-")
        a.text(0.98, 0.9, t, transform=a.transAxes, ha="right", va="top", fontsize=7.5)
        _ace_axes(a); a.set_ylim(-0.3, 10); a.set_yticks([0, 5, 10])
    for a in axs[:2]:
        a.tick_params(labelbottom=False); a.set_xlabel("")
    axs[1].set_ylabel("mean concentration (mM)")
    axs[0].text(-0.2, 1.06, "B", transform=axs[0].transAxes, fontweight="bold", fontsize=10)
    # tool legend under the biomass panel, metabolite legend under the concentration panels
    a0.legend(handles=_legend_handles(mets=False), fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3,
              frameon=False, handlelength=2.4)
    axs[2].legend(handles=_legend_handles(tools=False), fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.55), ncol=3,
                  frameon=False)
    _save(fig, "supp_fig_benchmark_aceref")


# ----------------------------------------------------------------------------- consortium
# Qualitative comparison, same biology in each tool's own representation
SPECIES = {"c_beijerinckii": "C. beijerinckii", "m_barkeri": "M. barkeri"}
# COMETS at 10 um boxes (scripts/comets_cfref_10um.py) with Convection 2D (pushing)
COMETS_CF = {"COMETS, convection": "conv"}
TOOLS["COMETS, convection"] = dict(color="#d1495b", ls="--")
SP_PC = {"CB": "c_beijerinckii", "MB": "m_barkeri"}
# The consortium scenario (criteria in results/cf_fig6_qualitative_criteria.md): the configuration of manuscript Fig. 6
# (CF3-NC50_s50, 50 h).
SCEN = {
    "fig6": dict(case="CF-Fig6", run="CF3-NC50_s50", series="cf3nc50_series.csv", t_end=50, comets_dir="fig6_comets",
                 comets_pre="comets_cfref_fig6_", ba_dir="fig6_bacarena", ba_ts="bacarena_fig6_80x80_closed_nopool_timeseries.csv",
                 ba_x="bacarena_fig6_closed_nopool_exchange_signs.csv", ba_l="bacarena_fig6_closed_nopool_mb_local.csv",
                 ba_timing="bacarena_fig6_80x80_closed_nopool_timing.json",
                 ba_alt={}, pc_note="dx 1.6 um, dt 0.001 min", out="fig6"),
}
EXM = ["glucose", "acetate", "H2", "CO2", "methane"]


def load_cf(sc):
    """Per tool: growth series per species (G), exchange fluxes at 10/20/30 h (X) and M. barkeri rate vs local H2 (L)."""
    G, X, L = {}, {}, {}
    T = sc["t_end"]; XT = [600 * k for k in range(1, T // 10 + 1)]          # exchange snapshots every 10 h (min)
    p = pd.read_csv(UC / "crossfeeding/results" / sc["series"])
    p = p[(p.run == sc["run"]) & (p.t <= 60 * T)]
    P = {}
    for sp, g in p.groupby("species"):
        g = g.sort_values("t")
        # fold_bio: biomass fold change estimated as exp(integral of the mean specific growth rate, 1/h), comparable
        # with the biomass-based COMETS / BacArena curves, unlike the cell count, which only changes at division
        gr = g.growth_mean.fillna(0).to_numpy(); th = (g.t / 60).to_numpy()
        fold_bio = np.exp(np.concatenate([[0], np.cumsum(0.5 * (gr[1:] + gr[:-1]) * np.diff(th))]))
        P[SP_PC[sp]] = pd.DataFrame(dict(t_h=th, fold=g.alive / g.alive.iloc[0], fold_bio=fold_bio, n=g.alive))
    G["PhysiCelldFBA"] = P
    X["PhysiCelldFBA"] = pd.DataFrame([dict(time_h=t / 60, species=SP_PC[r.species], metabolite=m, flux=getattr(r, "flux_" + m))
                                       for t in XT for r in p[p.t == t].itertuples() for m in EXM])
    mb = p[(p.species == "MB") & (p.t > 0)].sort_values("t")
    L["PhysiCelldFBA"] = pd.DataFrame(dict(time_h=mb.t / 60, growth_rate=mb.growth_mean, H2=mb.conc_H2))
    for tool, tag in COMETS_CF.items():
        d = D / sc["comets_dir"]; tag = sc["comets_pre"].replace("comets_cfref_", "") + tag
        f = d / f"comets_cfref_{tag}_biomass.csv"
        if not f.exists():
            print(f"{sc['case']}: {tool} skipped ({f.name} missing)"); continue
        c = pd.read_csv(f).sort_values("time_min")
        G[tool] = {sp: pd.DataFrame(dict(t_h=c.time_min / 60, fold=c[sp] / c[sp].iloc[0])) for sp in SPECIES}
        fd = d / f"comets_cfref_{tag}_domain_biomass.csv"
        if fd.exists():
            # biomass inside the domain only (Convection 2D pushes biomass into the fixed-concentration ring once the
            # colony reaches the edge); edge = first log with biomass in the ring
            c2 = pd.read_csv(fd).sort_values("time_min")
            ring = c2.c_beijerinckii_ring + c2.m_barkeri_ring
            edge = c2.time_min[ring > 1e-15].min() / 60
            G[tool] = {sp: pd.DataFrame(dict(t_h=c2.time_min / 60, fold=c2[sp + "_domain"] / c2[sp + "_domain"].iloc[0],
                                             fold_total=c[sp].to_numpy() / c[sp].iloc[0])) for sp in SPECIES}
            for sp in SPECIES:
                G[tool][sp].attrs["edge_h"] = edge if np.isfinite(edge) else None
        X[tool] = pd.read_csv(d / f"comets_cfref_{tag}_exchange_signs.csv")
        X[tool] = X[tool][X[tool].time_h.isin([t / 60 for t in XT])]
        q = pd.read_csv(d / f"comets_cfref_{tag}_mb_local.csv")
        q = q[q.species == "m_barkeri"].sort_values("time_min")
        # population rate of M. barkeri between logs (amended 2026-11-19): under Convection 2D the per-box ln-ratio
        # also counts biomass pushed between boxes
        rate = np.log(q.biomass_g.shift(-1) / q.biomass_g) / (q.time_min.diff().shift(-1) / 60)
        L[tool] = pd.DataFrame(dict(time_h=q.time_min / 60, growth_rate=rate, H2=q.H2_mM)).dropna()
        L[tool] = L[tool][L[tool].time_h > 0]
    bd = D / sc["ba_dir"]
    if not (bd / sc["ba_ts"]).exists():
        print(f"{sc['case']}: BacArena skipped ({sc['ba_ts']} missing)"); return G, X, L
    b = pd.read_csv(bd / sc["ba_ts"])
    G["BacArena"] = {sp: pd.DataFrame(dict(t_h=g.sort_values("step_h").step_h, fold=g.sort_values("step_h").biomass_pg / g.sort_values("step_h").biomass_pg.iloc[0],
                                           n=g.sort_values("step_h").n_cells)) for sp, g in b.groupby("species")}
    X["BacArena"] = pd.read_csv(bd / sc["ba_x"]).rename(columns={"step_h": "time_h"})
    q = pd.read_csv(bd / sc["ba_l"])
    q = q[q.species == "m_barkeri"].sort_values("step_h")
    # population rate between hourly steps (individuals cannot be followed across steps)
    rate = np.log(q.biomass_pg.shift(-1) / q.biomass_pg)
    L["BacArena"] = pd.DataFrame(dict(time_h=q.step_h, growth_rate=rate, H2=q.H2_mM)).dropna()
    L["BacArena"] = L["BacArena"][L["BacArena"].time_h > 0]
    for name, v in sc.get("ba_alt", {}).items():
        b = pd.read_csv(bd / sc["ba_ts"].replace("closed_nopool", v))
        G[name] = {sp: pd.DataFrame(dict(t_h=g.sort_values("step_h").step_h, fold=g.sort_values("step_h").biomass_pg / g.sort_values("step_h").biomass_pg.iloc[0],
                                         n=g.sort_values("step_h").n_cells)) for sp, g in b.groupby("species")}
        X[name] = pd.read_csv(bd / sc["ba_x"].replace("closed_nopool", v)).rename(columns={"step_h": "time_h"})
        q = pd.read_csv(bd / sc["ba_l"].replace("closed_nopool", v)); q = q[q.species == "m_barkeri"].sort_values("step_h")
        L[name] = pd.DataFrame(dict(time_h=q.step_h, growth_rate=np.log(q.biomass_pg.shift(-1) / q.biomass_pg), H2=q.H2_mM)).dropna()
        L[name] = L[name][L[name].time_h > 0]
    return G, X, L


def _spearman(a, b):
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def criteria_cf(sc, G, X, L):
    T = sc["t_end"]; nx = T // 10; tl = "/".join(str(10 * k) for k in range(1, nx + 1))
    out = []
    def add(tool, crit, value, ok):
        out.append(dict(tool=tool, criterion=crit, value=value, passed=bool(ok)))
    for tool in G:
        g = G[tool]
        fb = {sp: (g[sp].fold_bio if "fold_bio" in g[sp] else g[sp].fold).iloc[-1] for sp in SPECIES}
        add(tool, f"1 both species grow (fold 0-{T} h, C.b. / M.b.)", f"{fb['c_beijerinckii']:.3g} / {fb['m_barkeri']:.3g}",
            all(v > 1 for v in fb.values()))
        x = X[tool]
        sign = lambda sp, m: x[(x.species == sp) & (x.metabolite == m)].sort_values("time_h").flux.to_numpy()
        want = {"c_beijerinckii": {"glucose": -1, "acetate": 1, "H2": 1, "CO2": 1}, "m_barkeri": {"H2": -1, "acetate": -1, "methane": 1}}
        for k, sp in (("2", "c_beijerinckii"), ("3", "m_barkeri")):
            vals = {m: sign(sp, m) for m in want[sp]}
            ok = all(len(v) == nx and np.all(np.sign(v) == want[sp][m]) for m, v in vals.items())
            add(tool, f"{k} {SPECIES[sp]} exchange signs at {tl} h",
                "; ".join(f"{m} " + "/".join(f"{u:+.2g}" for u in v) for m, v in vals.items()), ok)
        l = L[tool].sort_values("time_h")
        r5 = l[l.time_h <= 5].growth_rate.iloc[-1]; rT = l.growth_rate.iloc[-1]
        add(tool, f"4a M. barkeri rate {T} h > 5 h (1/h)", f"{r5:.3g} -> {rT:.3g} (at {l.time_h.iloc[-1]:g} h)", rT > r5)
        rho = _spearman(l.growth_rate, l.H2)
        add(tool, "4b Spearman rho(M. barkeri rate, H2 at its locations)", f"{rho:.2f} (n={len(l)})", rho > 0)
        if tool.startswith("BacArena") and "n" in g["c_beijerinckii"]:
            # descriptive, chosen after the runs (not a criterion): lattice occupancy and rho while < 50% occupied
            occ = (g["c_beijerinckii"].set_index("t_h").n + g["m_barkeri"].set_index("t_h").n)
            nsites = 6400
            h50 = occ[occ >= nsites / 2].index.min()
            pre = l[l.time_h < h50] if np.isfinite(h50) else l
            out.append(dict(tool=tool, criterion="(descriptive, post hoc) lattice 50% full at (h); rho while < 50% full",
                            value=f"{h50:g} h; {_spearman(pre.growth_rate, pre.H2):.2f} (n={len(pre)})", passed=None))
    C = pd.DataFrame(out); C.to_csv(OUT / f"benchmark_{sc['out']}_criteria.csv", index=False)
    print(C.to_string(index=False))
    return C


def summarise_cf(sc, G):
    T, case = sc["t_end"], sc["case"]
    prov = pd.read_csv(UC / "crossfeeding/provenance.csv").set_index("run_id").loc[sc["run"]]
    timing = {"PhysiCelldFBA": (secs(prov.start, prov.end), f"{int(prov.threads)} thread(s), {prov.machine}, {sc['pc_note']}")}
    for tool, tag in COMETS_CF.items():
        f = D / sc["comets_dir"] / f"{sc['comets_pre']}{tag}_timing.json"
        if tool in G:
            j = json.load(open(f))
            timing[tool] = (j["elapsed_s"], f"{j['threads']} thread(s), bsclife006, {j['grid_um']:g} um boxes, dt {j['timestep_min']:g} min")
    if "BacArena" in G:
        j = json.load(open(D / sc["ba_dir"] / sc["ba_timing"]))
        timing["BacArena"] = (j["elapsed_s"], "1 thread(s), bsclife006, 2.5 um sites, dt 1 h")
        for name, v in sc.get("ba_alt", {}).items():
            j = json.load(open(D / sc["ba_dir"] / sc["ba_timing"].replace("closed_nopool", v)))
            timing[name] = (j["elapsed_s"], "1 thread(s), bsclife006, 2.5 um sites, dt 1 h")
    for tool, S2 in G.items():
        for sp, s in S2.items():
            basis = "biomass" if tool != "PhysiCelldFBA" else "cells"
            rec(case, tool, f"{SPECIES[sp]} fold change 0-{T} h ({basis})", s.fold.iloc[-1], "x")
            if "fold_total" in s:
                rec(case, tool, f"{SPECIES[sp]} fold change 0-{T} h (biomass incl. ring outside the domain)", s.fold_total.iloc[-1], "x")
                if s.attrs.get("edge_h") and sp == "c_beijerinckii":
                    rec(case, tool, "colony reaches the domain edge", s.attrs["edge_h"], "h", "first biomass in the ring")
                if s.attrs.get("edge_h"):
                    e = s[s.t_h <= s.attrs["edge_h"]].iloc[-1]  # noqa
                    rec(case, tool, f"{SPECIES[sp]} fold change at the edge time", e.fold, "x", f"{e.t_h:g} h")
            if "fold_bio" in s:
                rec(case, tool, f"{SPECIES[sp]} fold change 0-{T} h (biomass, from integrated mean growth rate)", s.fold_bio.iloc[-1], "x")
            if "n" in s:
                rec(case, tool, f"{SPECIES[sp]} cells or individuals at 0 / {T} h", f"{int(s.n.iloc[0])} / {int(s.n.iloc[-1])}", "")
        t, note = timing[tool]
        rec(case, tool, "wall time", t, "s", note)


SP_COL = {"c_beijerinckii": "#2a78d6", "m_barkeri": "#eb6834"}
CF_PANELS = {"PhysiCelldFBA": "PhysiCelldFBA", "COMETS, convection": "COMETS (Convection 2D)", "BacArena": "BacArena (closed arena)"}


def cf_biomass_pg(sc, G):
    """Dry biomass (pg) per tool and species. PhysiCelldFBA: founder dry mass (that given to COMETS) x fold change from the
    integrated mean growth rate; COMETS: biomass in the domain; BacArena: summed individual weights."""
    d = D / sc["comets_dir"]
    c0 = pd.read_csv(d / f"{sc['comets_pre']}conv_domain_biomass.csv").iloc[0]
    out = {}
    for tool in CF_PANELS:
        if tool not in G:
            continue
        out[tool] = {}
        for sp in SPECIES:
            s = G[tool][sp]
            if tool == "PhysiCelldFBA":
                y = s.fold_bio * c0[sp + "_domain"] * 1e12
            elif tool.startswith("COMETS"):
                y = s.fold * c0[sp + "_domain"] * 1e12
            else:
                b = pd.read_csv(D / sc["ba_dir"] / sc["ba_ts"]); b = b[b.species == sp].sort_values("step_h")
                y = b.biomass_pg.to_numpy()
            out[tool][sp] = (s.t_h.to_numpy(), np.asarray(y))
    return out


def fig_cf(sc, G, L):
    """One column per tool, both species in each panel; top row linear, bottom row log scale (shows M. barkeri)."""
    B = cf_biomass_pg(sc, G)
    fig, ax = plt.subplots(2, len(B), figsize=(2.3 * len(B), 4.2), sharex=True, sharey="row")
    for k, (tool, S2) in enumerate(B.items()):
        for r in range(2):
            a = ax[r, k]
            for sp, (t, y) in S2.items():
                a.plot(t, y, color=SP_COL[sp], label=f"${SPECIES[sp].replace(' ', '~')}$")
            e = G[tool]["c_beijerinckii"].attrs.get("edge_h")
            if e:
                a.axvline(e, color="0.75", lw=0.6, zorder=0)
            a.set_xlim(0, sc["t_end"])
        ax[0, k].set_ylim(bottom=0); ax[0, k].set_title(CF_PANELS[tool], fontsize=8)
        ax[1, k].set_yscale("log"); ax[1, k].set_xlabel("time (h)")
    ax[1, 0].yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10))
    ax[1, 0].yaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext(base=10))
    ax[0, 0].set_ylabel("dry biomass (pg)"); ax[1, 0].set_ylabel("dry biomass (pg, log scale)")
    ax[0, 0].legend(loc="upper left", fontsize=7)
    fig.tight_layout(); _save(fig, "supp_fig_benchmark_crossfeeding")


if __name__ == "__main__":
    A = load_ace(); summarise_ace(A); fig_ace(A)
    sc = SCEN["fig6"]
    G, X, L = load_cf(sc); summarise_cf(sc, G); criteria_cf(sc, G, X, L); fig_cf(sc, G, L)
    pd.DataFrame(rows).to_csv(OUT / "benchmark_summary.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
