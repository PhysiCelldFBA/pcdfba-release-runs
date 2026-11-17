"""Mapping of the supplementary tables of the revised manuscript.

OLD_NUMBER: number printed before the Supplementary Data workbook was introduced
(tracked PDF of 2026-11-17).  LABEL: LaTeX label.  SHEET: sheet name in the workbook
(None = stays only in the PDF).  WHERE: 'pdf', 'workbook' or 'both'.
The new PDF numbers are computed from the order, skipping tables that left the PDF.
"""
# (old number, label, where, sheet name)
TABLES = [
 (1,'unit_test','pdf',None),
 (2,'ut_design','workbook','D1'),
 (3,'ut_dxdt','workbook','D2'),
 (4,'ut_voxel_cell','workbook','D3'),
 (5,'ledger','pdf',None),
 (6,'dfba_unit_test_simulation','pdf',None),
 (7,'ecoli_acetic_switch_simulation','pdf',None),
 (8,'bacterial_colony_simulation','pdf',None),
 (9,'cancer_metabolism_simulation','pdf',None),
 (10,'crossfeeding_simulation','pdf',None),
 (11,'metabolic_driven_motility_simulation','pdf',None),
 (12,'ecoli_transport_summary','pdf',None),
 (13,'mcf7_cancer_summary','pdf',None),
 (14,'crossfeeding_transport_summary','pdf',None),
 (15,'substrate_diffusion_summary','pdf',None),
 (16,'param_provenance','pdf',None),
 (17,'free_uptake','pdf',None),
 (18,'tum_dt_dx','workbook','D4'),
 (19,'cf_dt_dx','workbook','D5'),
 (20,'cf_nc_dt','workbook','D6'),
 (21,'mot_conv','workbook','D7'),
 (22,'cap_summary','both','D8'),
 (23,'tum_grid','workbook','D9'),
 (24,'cf_sens','workbook','D10'),
 (25,'cost_usecases','pdf',None),
 (26,'cost_scaling','workbook','D11'),
 (27,'cost_feasible','workbook','D12'),
 (28,'zs','pdf',None),
 (29,'use_case_map','pdf',None),
]
def mapping():
    out={};n=0
    for old,label,where,sheet in TABLES:
        if where in('pdf','both'): n+=1
        out['tab:'+label]=dict(old_number='S%d'%old,where=where,sheet=sheet,
                               new_pdf_number=('S%d'%n if where in('pdf','both') else None))
    return out
if __name__=='__main__':
    import json,sys
    m=mapping(); json.dump(m,open('table_map.json','w'),indent=1)
    for k,v in m.items(): print(k,v)
