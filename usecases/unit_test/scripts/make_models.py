"""Create the model files of UT-ATPM and UT-IML1515 in models/ and record their checksums in models/MODELS.tsv.
usage: make_models.py IML1515_SOURCE.xml
  Ecoli_core_ATPM8.39.xml  = config/shipped/Ecoli_core.xml (the pinned code's model) with the lower bound of ATPM changed from 0 to 8.39 mmol/gDW/h (the value of the original BiGG core model); nothing else is changed.
  iML1515.xml              = a copy of the curated model used for the iML1515 repetition (source: dfba_sbml_models, see PhysiCelldFBA/config/iML1515.xml); ATPM lower bound 6.86 as supplied."""
import sys, re, hashlib, os, shutil
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
core = open(f"{UC}/config/shipped/Ecoli_core.xml").read()
pat = r'(<parameter id="R_ATPM_lower_bound" value=")0\.0(")'; assert len(re.findall(pat, core)) == 1
open(f"{UC}/models/Ecoli_core_ATPM8.39.xml", "w").write(re.sub(pat, r"\g<1>8.39\2", core))
shutil.copy(sys.argv[1], f"{UC}/models/iML1515.xml")
with open(f"{UC}/models/MODELS.tsv", "w") as f:
    f.write("file\tsha256\torigin\n")
    f.write(f"config/shipped/Ecoli_core.xml\t{sha(UC + '/config/shipped/Ecoli_core.xml')}\tpinned code (code/sample_projects_intracellular/fba/dfba_unit_test/config/Ecoli_core.xml), ATPM lower bound 0\n")
    f.write(f"models/Ecoli_core_ATPM8.39.xml\t{sha(UC + '/models/Ecoli_core_ATPM8.39.xml')}\tscripts/make_models.py: as above with the ATPM lower bound 8.39\n")
    f.write(f"models/iML1515.xml\t{sha(UC + '/models/iML1515.xml')}\tcurated iML1515 (dfba_sbml_models) [[Zenodo record]]; ATPM lower bound 6.86 as supplied\n")
print(open(f"{UC}/models/MODELS.tsv").read())
