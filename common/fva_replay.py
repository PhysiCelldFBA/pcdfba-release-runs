"""FVA replay of recorded cells. usage: fva_replay.py RUNDIR [--snap N|last] [--n 50 (per cell type)] [--seed 1] [--fractions 0.99999,0.9999] [--tol 1e-2] [--solver glpk] [--out PREFIX]
For a sample of living cells of one snapshot, rebuild the linear problem that the cell solved (the SBML model of its cell type with the lower bound of every recorded exchange set to the value in the *_lb column), solve FBA, and run flux variability analysis on the exchanges at fraction_of_optimum of the objective.
Per cell and fraction: the FBA objective of the replay (compared with the recorded growth_rate column), and for each exchange the recorded flux, the lower bound, the FVA range, whether it is determined (range <= tol mmol/gDW/h; the range at fraction 0.99999 of a unique optimum is itself about 1e-3, and it scales with 1 - fraction, which the summary reports as width_ratio) and whether the recorded flux is inside the range.
Cells with identical recorded bounds of the same type give the same problem and are solved once (column n_same). Writes PREFIX_cells.csv (one row per sampled cell and exchange) and PREFIX_summary.csv; default PREFIX = RUNDIR/fva_replay_<snapshot>.
Needs cobra; reads the model files named in RUNDIR/PhysiCell_settings.xml (relative to the use-case folder, then to RUNDIR). The replay does not use the PhysiCell binary or the Clp solver of the run: it is an independent check of the optimum and of its uniqueness."""
import sys, os, re, glob, argparse, hashlib, warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from pcxml import labels
import cobra


def rid(mod, k):
    """Reaction of the model for a recorded exchange name (cobra strips the 'R_' prefix of SBML ids), or None."""
    for n in (k, k[2:] if k.startswith("R_") else k):
        if n in mod.reactions: return mod.reactions.get_by_id(n)
    return None


def snapshots(d):
    return sorted(glob.glob(f"{d}/output[0-9]*[0-9].xml"))


def cell_models(settings, d):
    """{cell type index: model path}: the sbml_filename of each <cell_definition> that has a dfba <intracellular>."""
    root = ET.parse(settings).getroot(); out = {}
    for i, cd in enumerate(root.findall(".//cell_definitions/cell_definition")):
        s = cd.find(".//intracellular/sbml_filename") if cd.find(".//intracellular") is not None else None
        if s is None: s = cd.find(".//sbml_filename")
        if s is None or not s.text: continue
        name = s.text.strip(); uc = os.path.dirname(os.path.dirname(os.path.abspath(d)))
        for base in (uc, d, os.path.join(uc, "config", "shipped")):
            p = os.path.join(base, name)
            if os.path.exists(p): out[int(cd.get("ID", i))] = os.path.normpath(p); break
        else: raise FileNotFoundError(f"{name} (cell type {cd.get('name')})")
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--snap", default="last"); ap.add_argument("--n", type=int, default=50); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fractions", default="0.99999,0.9999"); ap.add_argument("--tol", type=float, default=1e-2); ap.add_argument("--solver", default="glpk"); ap.add_argument("--out", default=None)
    a = ap.parse_args(); d = a.rundir.rstrip("/"); fr = [float(x) for x in a.fractions.split(",")]
    xs = snapshots(d); xf = xs[-1] if a.snap == "last" else f"{d}/output{int(a.snap):08d}.xml"
    root = ET.parse(xf).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
    m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; col = lambda k: m[lab.index(k)]
    ex = [k for k in lab if k.startswith("R_EX") and not k.endswith("_lb")]; assert all(k + "_lb" in lab for k in ex), "no _lb columns: the run is not on the pinned code"
    types = np.asarray(col("cell_type")).astype(int); alive = np.asarray(col("dead")) == 0
    models = cell_models(f"{d}/PhysiCell_settings.xml", d); print(f"{xf}: t = {t:g} min, {m.shape[1]} cells, {int(alive.sum())} living, exchanges {len(ex)}, models {sorted(set(models.values()))}")
    rng = np.random.default_rng(a.seed); idx = np.sort(np.concatenate([(lambda ix: rng.choice(ix, min(a.n, len(ix)), replace=False))(np.flatnonzero(alive & (types == ty))) for ty in sorted(set(models))]))
    cache_model = {}
    def load(p):
        if p not in cache_model:
            mod = cobra.io.read_sbml_model(p); mod.solver = a.solver; cache_model[p] = mod
        return cache_model[p]
    cache = {}; rows = []
    for j in idx:
        ty = int(types[j]); mp = models.get(ty) or next(iter(models.values())); mod = load(mp)
        rec = {k: float(col(k)[j]) for k in ex}; lb = {k: float(col(k + "_lb")[j]) for k in ex}; gr = float(col("growth_rate")[j])
        key = (mp, tuple(round(v, 10) for v in lb.values()))
        if key not in cache:
            res = {}
            with mod as M:
                R = {k: rid(M, k) for k in ex}; R = {k: v for k, v in R.items() if v is not None}; res["keys"] = list(R)
                for k in R: R[k].lower_bound = lb[k]
                sol = M.slim_optimize(error_value=np.nan); res["obj"] = sol
                if not np.isnan(sol):
                    fl = M.optimize().fluxes
                    for f in fr:
                        r = cobra.flux_analysis.flux_variability_analysis(M, reaction_list=list(R.values()), fraction_of_optimum=f, processes=1)
                        r.index = list(R); res[f] = r; res["fba"] = pd.Series({k: fl[R[k].id] for k in R})
            res.setdefault("keys", list(ex)); cache[key] = [res, 0]
        cache[key][1] += 1; res = cache[key][0]
        for f in fr:
            for k in ex:
                if k not in res.get("keys", ex): continue  # exchange of the other species
                r = res.get(f)
                lo, hi = (r.loc[k, "minimum"], r.loc[k, "maximum"]) if r is not None and k in r.index else (np.nan, np.nan)
                rows.append(dict(run=os.path.basename(d), t=t, cell=int(col("ID")[j]), type=ty, x=float(col("position[0]")[j]), y=float(col("position[1]")[j]), z=float(col("position[2]")[j]), growth_rec=gr, obj_fba=res["obj"],
                                 fraction=f, exchange=k, lb=lb[k], rec=rec[k], fva_min=lo, fva_max=hi, width=hi - lo, determined=bool(hi - lo <= a.tol) if r is not None else None,
                                 rec_in_range=bool(lo - a.tol <= rec[k] <= hi + a.tol) if r is not None else None, fba_flux=float(res["fba"][k]) if "fba" in res and k in res["fba"].index else np.nan, key=hashlib.md5(repr(key).encode()).hexdigest()[:8]))
    df = pd.DataFrame(rows); df["n_same"] = df.key.map({hashlib.md5(repr(k).encode()).hexdigest()[:8]: v[1] for k, v in cache.items()})
    out = a.out or f"{d}/fva_replay_{os.path.basename(xf)[6:14]}"; df.to_csv(out + "_cells.csv", index=False)
    ok = df.drop_duplicates("cell"); feas = ok.obj_fba.notna()
    print(f"cells {len(ok)}, distinct problems {len(cache)}, feasible {int(feas.sum())}; replay objective (1/h): median {np.nanmedian(ok.obj_fba):.4g}; |objective - recorded growth_rate| / recorded: median {np.nanmedian((ok.obj_fba - ok.growth_rec).abs() / ok.growth_rec.abs().replace(0, np.nan)):.2g}, max {np.nanmax((ok.obj_fba - ok.growth_rec).abs() / ok.growth_rec.abs().replace(0, np.nan)):.2g}")
    sm = []; W = df.pivot_table(index=["cell", "exchange"], columns="fraction", values="width") if len(fr) > 1 else None
    for (f, k), g in df.groupby(["fraction", "exchange"]):
        g = g[g.obj_fba.notna()]; used = g[(g.rec.abs() > a.tol) | (g.width > a.tol)]
        sm.append(dict(run=os.path.basename(d), t=t, fraction=f, exchange=k, cells=len(g), active=len(used), determined_pct=100 * g.determined.mean() if len(g) else np.nan, median_width=g.width.median(), max_width=g.width.max(),
                       median_width_rel_pct=100 * (g.width / g.rec.abs().replace(0, np.nan)).median(), rec_in_range_pct=100 * g.rec_in_range.mean() if len(g) else np.nan, median_rec=g.rec.median(),
                       width_ratio=float(np.nanmedian((W[fr[1]] / W[fr[0]].replace(0, np.nan)).xs(k, level="exchange"))) if W is not None and f == fr[0] and len(g) else np.nan,
                       infeasible=int((df[df.fraction == f].drop_duplicates("cell").obj_fba.isna()).sum()), infeasible_and_zero_growth=int(((df[df.fraction == f].drop_duplicates("cell").obj_fba.isna()) & (df[df.fraction == f].drop_duplicates("cell").growth_rec == 0)).sum())))
    S = pd.DataFrame(sm); S.to_csv(out + "_summary.csv", index=False); pd.set_option("display.width", 200); print(S[S.fraction == fr[0]].drop(columns=["run", "t", "fraction"]).round(5).to_string(index=False))
    print("wrote", out + "_cells.csv", out + "_summary.csv")


if __name__ == "__main__": main()
