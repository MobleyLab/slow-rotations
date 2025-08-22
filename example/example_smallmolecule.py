# an example for comparing simulations repeats of the same system

from slow_rotations import torsions as tor
from slow_rotations import rdkit_wrapper as rdw
from slow_rotations import molconverter as mc
from slow_rotations import mappings
from slow_rotations import compare

import warnings
import json

tf_list = []
for rpt in range(3):
	print("Loading repeat {rpt}")
	lmda=0
	topf_bnd = 'traj.gro'
	trajf_bnd = f'traj_{rpt+1}.xtc'

	smiles = "[H]c1c(c(c(c(c1C(=O)O[H])O[H])[H])N([H])[H])[H]" #ZINC922

	ligcode1 = "LIG"
	ligcode2 = "UNL"
	ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode1,smiles)

	tf_list.append(ligtor_bnd)



ligcomp = compare.LigandTorsionComparator(tf_list)

torsions = ligcomp.get_torsions()

results = {}
for idx,t in enumerate(torsions):
	imgname = f'{"_".join(map(str, t))}.png'
	t_result = ligcomp.plot_all_distributions(t,save_path=f"../example/{imgname}")
	results[f't{idx}'] = t_result


with open("../example/torsiondata.json", "w") as f:
	json.dump(results, f)

