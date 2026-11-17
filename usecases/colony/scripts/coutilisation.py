"""Cells that take up glucose and acetate at the same time, colony reference run, from one snapshot.
usage: coutilisation.py RUNDIR T_H   (e.g. runs/COL-REF_s50 9)
A cell takes up a substrate when its exchange flux (mmol/gDW/h) is below -THR (THR = 0.01); living cells only.
Used for the sentence 'all of them also took up glucose' (main text, colony) and the Discussion."""
import sys, glob, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET
run, T, THR = sys.argv[1].rstrip('/'), float(sys.argv[2]), 0.01
files = sorted(glob.glob(f'{run}/output*[0-9].xml'))
tof = lambda f: float(ET.parse(f).getroot().find('.//current_time').text) / 60
f = min(files, key=lambda f: abs(tof(f) - T)); root = ET.parse(f).getroot()
idx = {x.text: int(x.get('index')) for x in root.iter('label')}
m = sio.loadmat(f.replace('.xml', '_cells.mat'))['cells']
al = m[idx['dead']] == 0; glc = m[idx['R_EX_glc__D_e']] < -THR; ac = m[idx['R_EX_ac_e']] < -THR
print(f'{f.split("/")[-1]}, t = {tof(f):.2f} h, living {al.sum()}, glucose uptake {(glc&al).sum()}, acetate uptake {(ac&al).sum()} '
      f'({(ac&al).sum()/al.sum():.1%}), both {(glc&ac&al).sum()}, acetate without glucose {(ac&~glc&al).sum()}')
