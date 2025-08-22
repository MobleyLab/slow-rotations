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
	topf_bnd = 'traj.gro'
	trajf_bnd = f'traj_{rpt+1}.xtc'

	smiles = "[H]c1c(c(c(c(c1C(=O)O[H])O[H])[H])N([H])[H])[H]" #ZINC922

	ligcode1 = "LIG"
	ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode1,smiles)


	tf_list.append(pro_tf)



comp = compare.ProteinTorsionComparator(tf_list, 5)

torsions = comp.get_torsions()

results = {}
for idx,t in enumerate(torsions):
	imgname = f'{"_".join(map(str, t))}.png'
	t_result = ligcomp.plot_all_distributions(t,save_path=f"{imgname}")
	results[f't{idx}'] = t_result

import json
with open("side_chain_torsiondata.json", "w") as f:
	json.dump(results, f)