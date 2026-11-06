"""Which of the capped internal bounds of the shipped M. barkeri limit growth in the cross-feeding medium (companion of mb_caps.py): relax one reaction or a group of the changed bounds to the wider of the shipped and the full iMG746 value and report the growth ratio for the first 20 M. barkeri cells of the 30 and 60 snapshots of CF-REF_s50. Run from usecases/crossfeeding."""
import sys,os,numpy as np,scipy.io as sio,xml.etree.ElementTree as ET,warnings,cobra; warnings.filterwarnings("ignore")
sys.path.insert(0,'../../common'); from pcxml import labels; from fva_replay import cell_models, rid
d='runs/CF-REF_s50'; models=cell_models(f'{d}/PhysiCell_settings.xml',d); ty=[t for t,p in models.items() if 'barkeri' in p][0]
Mt=cobra.io.read_sbml_model(models[ty]); Mf=cobra.io.read_sbml_model('../unit_test/models_full/M_barkeri_model_iMG746.xml'); Mt.solver='gurobi'
chg=[r.id for r in Mt.reactions if not r.boundary and r.bounds!=Mf.reactions.get_by_id(r.id).bounds]
for s in (30,60):
    xf=f'{d}/output{s:08d}.xml'; root=ET.parse(xf).getroot(); lab=labels(root); m=sio.loadmat(xf.replace('.xml','_cells.mat'))['cells']; col=lambda k:m[lab.index(k)]
    ex=[k for k in lab if k.startswith('R_EX') and not k.endswith('_lb')]; ok=np.flatnonzero((np.asarray(col('cell_type')).astype(int)==ty)&(np.asarray(col('dead'))==0))[:20]
    def solve(relax):
        out=[]
        for j in ok:
            with Mt as M:
                for i in relax:
                    r=M.reactions.get_by_id(i); f=Mf.reactions.get_by_id(i); r.bounds=(min(r.lower_bound,f.lower_bound),max(r.upper_bound,f.upper_bound))
                for k in ex:
                    r=rid(M,k)
                    if r is not None: r.lower_bound=float(col(k+'_lb')[j])
                out.append(M.slim_optimize())
        return np.array(out)
    base=solve([]); print('RES t',s,'min; base growth median',round(float(np.median(base)),5))
    for name,rel in [('all caps',chg),('SO3t',['SO3t']),('GF4GL_6',['GF4GL_6']),('GF4GL_5',['GF4GL_5']),('OORr',['OORr']),('CF3Sa+CF3Sg',['CF3Sa','CF3Sg']),('all six active',['GF4GL_6','SO3t','OORr','CF3Sa','CF3Sg','GF4GL_5']),('all caps except the six',[i for i in chg if i not in ['GF4GL_6','SO3t','OORr','CF3Sa','CF3Sg','GF4GL_5']])]:
        r=solve(rel)/base; print(f'RES   relax {name}: growth x{np.median(r):.3f} (min {r.min():.3f}, max {r.max():.3f})')
