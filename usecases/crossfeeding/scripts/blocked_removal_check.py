"""Is the removal of blocked reactions and gap metabolites safe in the cross-feeding medium? Full model (exchanges of the trimmed file, other exchanges closed, recorded lower bounds, ATPM lower bound 0) against the trimmed model with the capped internal bounds relaxed, same recorded bounds, first 40 living cells of each species at 300, 900, 1800 min of CF-REF_s50. Run from usecases/crossfeeding."""
import sys,os,numpy as np,scipy.io as sio,xml.etree.ElementTree as ET,warnings,cobra; warnings.filterwarnings("ignore")
sys.path.insert(0,'../../common'); from pcxml import labels; from fva_replay import cell_models, rid
d='runs/CF-REF_s50'; models=cell_models(f'{d}/PhysiCell_settings.xml',d)
full={'barkeri':'../unit_test/models_full/M_barkeri_model_iMG746.xml','beijerinckii':'../unit_test/models_full/iCB925.xml'}
for ty,p in models.items():
    key=[k for k in full if k in p][0]; Mt=cobra.io.read_sbml_model(p); Mf=cobra.io.read_sbml_model(full[key]); Mt.solver=Mf.solver='gurobi'
    # relaxed trimmed (caps removed: changed internal bounds widened to the wider of the two); ATPM lower bound 0 in all
    low={r.id.upper():r for r in Mf.reactions}
    for r in Mf.reactions:
        if r.boundary: r.bounds=(0,0)
    for r in Mt.reactions:
        f=low.get(r.id.upper())
        if f is None: continue
        if r.boundary: f.bounds=r.bounds
        else: r.bounds=(min(r.lower_bound,f.lower_bound),max(r.upper_bound,f.upper_bound))
    for M in (Mt,Mf):
        if 'ATPM' in M.reactions: M.reactions.ATPM.lower_bound=0
    for s in (10,30,60):
        xf=f'{d}/output{s:08d}.xml'; root=ET.parse(xf).getroot(); lab=labels(root); m=sio.loadmat(xf.replace('.xml','_cells.mat'))['cells']; col=lambda k:m[lab.index(k)]
        ex=[k for k in lab if k.startswith('R_EX') and not k.endswith('_lb')]; ok=np.flatnonzero((np.asarray(col('cell_type')).astype(int)==ty)&(np.asarray(col('dead'))==0))[:40]
        g=[]
        for j in ok:
            v=[]
            for M in (Mt,Mf):
                with M as MM:
                    for k in ex:
                        r=rid(MM,k) or rid(MM,k.replace('R_EX','R_ex'))
                        if r is not None: r.lower_bound=float(col(k+'_lb')[j])
                    v.append(MM.slim_optimize(error_value=np.nan))
            g.append(v)
        g=np.array(g); r=g[:,1]/g[:,0]
        print(f'RES {key} t={s*30} min cells {len(ok)}: full model / trimmed-uncapped growth median {np.nanmedian(r):.6f} (min {np.nanmin(r):.6f}, max {np.nanmax(r):.6f}); trimmed growth median {np.nanmedian(g[:,0]):.4f}')
