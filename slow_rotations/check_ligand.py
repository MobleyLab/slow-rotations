import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path

from scipy.special import kl_div

base_path = Path('/Users/megosato/Desktop/simulations/mpro/md_traj/')

mpro_ligs = {
    "x11641_0A": "CC(C(=O)NC=1C=NC=CC1C)C=2C=C(Cl)C=C(OC3CC(=O)N3)C2"
}

imgs_save_path = Path('/Users/megosato/Desktop/images/')
imgs_save_path.mkdir(exist_ok=True)

for ligname in mpro_ligs:
	system = f'{ligname}'


	# gas is considered truth
	topf = base_path / f'ligand/{ligname}/nvt_prod.gro'
	trajf = base_path / f'ligand/{ligname}/prod.xtc'

	smiles = mpro_ligs[ligname]

	# notes CL should be Cl NOT CL
	# remove x from atom names of gmx topology
	# atom numbering cannot have no space next to atom type

	ligcode = "UNL"
	ligtor = tor.LigandTorsionFinder(str(trajf),str(topf),ligcode,smiles)

	lig_img_save_path = imgs_save_path / ligname
	lig_img_save_path.mkdir(exist_ok=True)

	for tor in ligtor.get_torsions():

		a1,a2,a3,a4 = tuple(tor)

		#X, score, angle_min = ligtor.get_kde(tor, num_bins=60, angle_min=None)

		tor_lig_img_save_path = lig_img_save_path / f'{system}-{a1}_{a2}_{a3}_{a4}'

		ligtor.make_torsion_img(tor, show=False, save_path = str(tor_lig_img_save_path))
		ligtor.plot_dihedral_scatter(tor, save_path = str(tor_lig_img_save_path)+"scatter", show=False)
