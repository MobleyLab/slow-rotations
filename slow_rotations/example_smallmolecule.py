# an example for comparing simulations repeats of the same system

import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path

from scipy.special import kl_div


tf_list = []
for rpt in range(2):
	lmda=0
	topf_bnd = '/Users/megosato/Desktop/traj.gro'
	trajf_bnd = f'/Users/megosato/Desktop/edge_ZINC922_ZINC337835/step1/{rpt+1}/lambda_{lmda}/Production_MD/traj.xtc'

	smiles = "[H]c1c(c(c2c(c1[H])C(C(C2([H])[H])([H])C(=O)O[H])([H])[H])[H])[H]" #ZINC337835

	ligcode1 = "LIG"
	ligcode2 = "UNL"
	ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode2,smiles)

	tf_list.append(ligtor_bnd)



ligcomp = compare.LigandTorsionComparator(tf_list)

torsions = ligcomp.get_torsions()

results = {}
for idx,t in enumerate(torsions):
	t_result = ligcomp.plot_all_distributions(t,save_path=f"/Users/megosato/Desktop/torsions/{t}.png")
	results[f't{idx}'] = t_result

import json
with open("/Users/megosato/Desktop/torsions/torsiondata.json", "w") as f:
	json.dump(results, f)

