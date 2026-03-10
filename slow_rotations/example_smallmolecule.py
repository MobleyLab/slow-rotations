# an example for comparing simulations repeats of the same system

from slow_rotations import torsions as tor
from slow_rotations import rdkit_wrapper as rdw
from slow_rotations import molconverter as mc
from slow_rotations import mappings
from slow_rotations import compare

import warnings
import json

tf_list = []
for rpt in range(1,3):
	print("Loading repeat {rpt}")
	topf_bnd = '/Users/megosato/Desktop/lig3_flex_move/complex_wat.prmtop'
	trajf_bnd = f'/Users/megosato/Desktop/lig3_flex_move/lig3_{rpt+1}.nc'

	smiles = "C1(=NC(=CC(=C1C#N)N([H])[H])N(C(CC2=CC(=CC=C2C)C)=O)[H])OCC"

	ligcode1 = "LIG"
	ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode1,smiles)

	tf_list.append(ligtor_bnd)



ligcomp = compare.LigandTorsionComparator(tf_list)

torsions = ligcomp.get_torsions()

results = {}
for idx,t in enumerate(torsions):
	if t[0] != 1:
		continue
	imgname = f'{"_".join(map(str, t))}.png'
	t_result = ligcomp.plot_all_distributions(t,save_path=f"../example/{imgname}")
	results[f't{idx}'] = t_result


with open("../example/torsiondata.json", "w") as f:
	json.dump(results, f)

