import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
from simulation_db import *
import mappings
import warnings
import compare
from pathlib import Path
import numpy as np

from pymbar import timeseries
from scipy.special import kl_div


imgs_save_path = Path('/Users/megosato/Desktop/images/')
imgs_save_path.mkdir(exist_ok=True)


results_list = list()

for ligname in AllSims_DB:

	print(ligname)

	smiles = AllSims_DB[ligname]['smiles']
	topf = AllSims_DB[ligname]['topf']
	trajf = AllSims_DB[ligname]['trajf']

	ligcode = AllSims_DB[ligname]['ligcode']
	ligtor = tor.LigandTorsionFinder(str(trajf),str(topf),ligcode,smiles)

	if ligname.count('_') == 4:

		system = ligname.split('_')[0]

		lig_img_save_path = imgs_save_path / system
		lig_img_save_path.mkdir(exist_ok=True)

		t = [int(i) - 1 for i in ligname.split('_')[1:]]
		a1,a2,a3,a4 = tuple(t)

		tor_lig_img_save_path = lig_img_save_path / f'{system}-{a1}_{a2}_{a3}_{a4}_restr'

		am, sa = ligtor.shift_torsion_angles(t)
		#np.savetxt(f'/Users/megosato/Desktop/{ligname}-{a1}_{a2}_{a3}_{a4}_restr.dat', sa)

		g = timeseries.statisticalInefficiency(sa.flatten())

		min_max, transition_ctr, states_list = ligtor.make_torsion_img(t, show=False, save_path = str(tor_lig_img_save_path))
		ligtor.plot_dihedral_scatter(t, save_path = str(tor_lig_img_save_path)+"scatter", show=False)
		ligtor.plot_dihedral_histogram(t, angle_min=-180, save_path=str(tor_lig_img_save_path)+"hist", show=False)

		results_list.append(
			{	
				'name': ligname,
				'num_states': len(states_list),
				'min_transitions': transition_ctr.min_transitions(),
				'max_transitions': transition_ctr.max_transitions(),
				'tot_transitions': transition_ctr.total_transitions(), 
				'g': g,
				'check_symmetry': int(mappings.check_symmetry(ligtor.oemol, t)),
			}
		)

	else:
		for t in ligtor.get_torsions():

			system = ligname
			lig_img_save_path = imgs_save_path / system
			lig_img_save_path.mkdir(exist_ok=True)

			a1,a2,a3,a4 = tuple(t)

			am, sa = ligtor.shift_torsion_angles(t)

			tor_lig_img_save_path = lig_img_save_path / f'{system}-{a1}_{a2}_{a3}_{a4}'

			g = timeseries.statisticalInefficiency(sa.flatten())

			print("G:", g)

			min_max, transition_ctr, states_list = ligtor.make_torsion_img(t, show=False, save_path = str(tor_lig_img_save_path))
			ligtor.plot_dihedral_scatter(t, save_path = str(tor_lig_img_save_path)+"scatter", show=False)
			ligtor.plot_dihedral_histogram(t, angle_min=-180, save_path=str(tor_lig_img_save_path)+"hist", show=False)

			results_list.append(
				{	
					'name': ligname,
					'num_states': len(states_list),
					'min_transitions': transition_ctr.min_transitions(),
					'max_transitions': transition_ctr.max_transitions(),
					'tot_transitions': transition_ctr.total_transitions(), 
					'g': g,
					'check_symmetry': int(mappings.check_symmetry(ligtor.oemol, t)),
				}
			)



import csv

with open('names.csv', 'w', newline='') as csvfile:
    fieldnames = [ k for k in results_list[0] ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()

    for result in results_list:
    	writer.writerow(result)
