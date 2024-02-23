import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path

from scipy.special import kl_div

base_path = Path('/Users/megosato/Desktop/simulations/bace_a/')



# gas is considered truth
topf_gas = base_path / f'solvent/confout.gro'
trajf_gas = base_path / f'solvent/traj_comp.xtc'

# bound is inside the protein or in this case, restrained
topf_bnd = base_path / f'confout.gro'
trajf_bnd = base_path / f'traj_comp.xtc'

smiles = "[N+](=C1N[C@](C(=O)N1C)(C2=CC=CC=C2)C3=CC(=CC=C3)C4=CN=CC=C4)([H])[H]"

# notes CL should be Cl NOT CL
# remove x from atom names of gmx topology
# atom numbering cannot have no space next to atom type

ligcode1 = "LIG"
ligcode2 = "LIG"
ligtor_gas = tor.LigandTorsionFinder(str(trajf_gas),str(topf_gas),ligcode1,smiles)
ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode2,smiles)

lc = compare.LigandComparator(ligtor_gas, ligtor_bnd)

imgs_save_path = Path('/Users/megosato/Desktop/images_bace_a/compare_solv')
imgs_save_path.mkdir(exist_ok=True)

lig_img_save_path = imgs_save_path
lig_img_save_path.mkdir(exist_ok=True)

with open(str(lig_img_save_path/'comparable.txt'), 'w') as comp_f:

	for g in lc.get_gas_torsions():

		print("TORSION:", g)

		X, score, angle_min = ligtor_gas.get_kde(g, num_bins=60, angle_min=None)
		nump=ligtor_gas.get_kde_num_peaks(score)

		result = lc.compare_torsions(g)
		a1,a2,a3,a4 = tuple(g)

		tor_lig_img_save_path = lig_img_save_path / f'bace_a-{a1}_{a2}_{a3}_{a4}'

		lc.make_compare_torsions_img(result, show=False, save_path = str(tor_lig_img_save_path))
		print(result.to_dict())
