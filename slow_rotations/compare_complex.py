import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import matplotlib.pyplot as plt
from pathlib import Path

# for protein we just care that there are sufficient transitions 

ANGSTROM_CUTOFF = 5.0
base_path = Path('/Users/megosato/Desktop/mpro/md_traj/')

ligname = 'z4ylu'





imgs_save_path = Path('/Users/megosato/Desktop/protein_images_compare/')
imgs_save_path.mkdir(exist_ok=True)

aa_img_save_path = imgs_save_path / ligname
aa_img_save_path.mkdir(exist_ok=True)

sysnames = []
ptfs = []
torsions = [[331,332,333,334], [335,336,337,338], [335,336,337,338], [335,336,337,338]]

for pose in ['0A', '0B', '0C', '0D']:
	system = f'{ligname}_{pose}_HIE167'
	# bound is inside the protein or in this case, restrained
	topf = base_path / f'complex/{system}/npt.gro'
	trajf = base_path / f'complex/{system}/prod.xtc'
	ligcode = "UNK"
	ptf = tor.ProteinTorsionFinder(str(trajf), str(topf), ligcode)
	sysnames.append(system)
	ptfs.append(ptf)



for i,ptf in enumerate(ptfs):
	# X, scores, angle_min = ptf.get_kde(t)
	# num_peaks = ptf.get_kde_num_peaks(scores, smoothing_window=100, peak_prominence=0.008)
	# ptf.plot_dihedral_histogram(t, show=False, title=f"# peaks = {num_peaks}")
	# min_max = ptf.get_bounds_knn(X, num_peaks)
	# angles = ptf.shift_torsion_angles(t, angle_min=angle_min)[1].flatten()

	# transition_matrix = ptf.transition_counter(angles, min_max)

	# print(transition_matrix)

	# ptf.check_transitions(transition_matrix, 10)

	t = torsions[i]

	a1,a2,a3,a4 = tuple(t)
	tor_aa_img_save_path = aa_img_save_path / f'{sysnames[i]}-{a1}_{a2}_{a3}_{a4}'

	if i == 0:
		angle_min = ptf.make_aa_torsions_img(t, save_path=str(tor_aa_img_save_path))
	else:
		print("in else")
		ptf.make_aa_torsions_img(t, save_path=str(tor_aa_img_save_path), angle_min=angle_min)