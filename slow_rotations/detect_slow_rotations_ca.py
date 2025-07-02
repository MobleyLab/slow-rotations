import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path
import json

from scipy.special import kl_div

# bound is inside the protein or in this case, restrained
topf_bnd = 'mobley_2913224/traj.gro'
trajf_bnd = 'mobley_2913224/traj.xtc'

smiles = "[H]c1c(c(c(c(c1[H])C(=O)O[H])OC(=O)C([H])([H])[H])[H])[H]" 

lig='mobley_2913224'

ligcode1 = "UNK"
ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode1,smiles)

imgs_save_path = Path('mobley_2913224/')
imgs_save_path.mkdir(exist_ok=True)

aa_img_save_path = imgs_save_path / lig
aa_img_save_path.mkdir(exist_ok=True)

results_dict = {}

for t in ligtor_bnd.get_torsions():
	try:
		t_tup = str(tuple(t))
		a1,a2,a3,a4 = tuple(t)

		tor_aa_img_save_path = aa_img_save_path / f'{lig}-{a1}_{a2}_{a3}_{a4}'

		ligtor_bnd.plot_kde(t, save_path=str(tor_aa_img_save_path)+"kde.png")

		ligtor_bnd.plot_dihedral_scatter(t, save_path=str(tor_aa_img_save_path)+"scatter.png", angle_min=-180)

		min_max, transition_ctr, states_list, symmetry, populations = ligtor_bnd.make_torsion_img(t, save_path=str(tor_aa_img_save_path))
		
		results_dict[t_tup] = dict()
		results_dict[t_tup]['nstates'] = transition_ctr.num_states
		results_dict[t_tup]['min_transitions'] = transition_ctr.min_transitions()
		results_dict[t_tup]['max_transitions'] = transition_ctr.max_transitions()
		results_dict[t_tup]['has_symmetry'] = populations
		results_dict[t_tup]['populations'] = symmetry

		for i in range(len(states_list)):
			results_dict[t_tup][f'{i}_in'] = transition_ctr.count_transitions_into_state(i)
			results_dict[t_tup][f'{i}_out'] = transition_ctr.count_transitions_out_of_state(i)

	

	except Exception as e:
		print(f"ERROR: {t}\n{e}")

imgs_save_path = Path(f'mobley_2913224/results.json')
with open(str(imgs_save_path), "w") as json_file:
	json.dump(results_dict, json_file, indent=4)