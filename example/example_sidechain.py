# an example for comparing simulations repeats of the same system

from slow_rotations import torsions as tor
from slow_rotations import rdkit_wrapper as rdw
from slow_rotations import molconverter as mc
from slow_rotations import mappings
from slow_rotations import compare

import warnings
import json


tf_list = []
for rpt in range(2):
	topf_bnd = 'traj.gro'
	trajf_bnd = f'traj_{rpt+1}.xtc'

	ligcode1 = "LIG"
	ligtor_bnd = tor.ProteinTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode1)

	tf_list.append(pro_tf)



comp = compare.ProteinTorsionComparator(tf_list, 5)

torsions = comp.get_torsions()

results = {}
for idx,t in enumerate(torsions):
	imgname = f'{"_".join(map(str, t))}.png'
	t_result = ligcomp.plot_all_distributions(t,save_path=f"{imgname}")
	results[f't{idx}'] = t_result

import json
with open("sidechain_torsiondata.json", "w") as f:
	json.dump(results, f)