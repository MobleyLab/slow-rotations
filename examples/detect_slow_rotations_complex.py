import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import matplotlib.pyplot as plt

# for protein we just care that there are sufficient transitions 

ANGSTROM_CUTOFF = 5.0

topf = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Meghan_Traj/pEH_lGVG_2_0_c0-prod00-centered.pdb'
trajf = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Meghan_Traj/pEH_lGVG_2_0_c0-prod00-centered.dcd'
ligcode = "LIG"

ptf = tor.ProteinTorsionFinder(trajf, topf, ligcode)


torsions = ptf.get_chi_x_torsions(1, a_cutoff = 5.0)

for t in torsions:
	X, scores, angle_min = ptf.get_kde(t)
	num_peaks = ptf.get_kde_num_peaks(scores, smoothing_window=100, peak_prominence=0.008)
	ptf.plot_dihedral_histogram(t, show=False, title=f"# peaks = {num_peaks}")
	min_max = ptf.get_bounds_knn(X, num_peaks)
	angles = ptf.shift_torsion_angles(t, angle_min=angle_min)[1].flatten()

	transition_matrix = ptf.transition_counter(angles, min_max)

	print(transition_matrix)

	ptf.check_transitions(transition_matrix, 10)


plt.show()