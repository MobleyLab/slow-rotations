import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path
import numpy as np

from scipy.special import kl_div

base_path = Path('/Users/megosato/Desktop/butylbenzene')

# gas is considered truth
topf = base_path / f'trj.pdb' #"/Users/megosato/Desktop/trj.pdb" #
trajf = base_path / f'trj.xtc'

smiles = "CCCCC1=CC=CC=C1"

# notes CL should be Cl NOT CL
# remove x from atom names of gmx topology
# atom numbering cannot have no space next to atom type

ligcode = "MOL0"
ligtor = tor.LigandTorsionFinder(str(trajf),str(topf),ligcode,smiles)

lig_img_save_path = base_path / "images"
lig_img_save_path.mkdir(exist_ok=True)

for t in ligtor.get_torsions() + [[9,8,6,17]]:

	print('TORSION', t)

	a1,a2,a3,a4 = tuple(t)

	am, sa = ligtor.shift_torsion_angles(t)
	#np.savetxt(f'/Users/megosato/Desktop/{ligname}-{a1}_{a2}_{a3}_{a4}.dat', sa)


	#X, score, angle_min = ligtor.get_kde(tor, num_bins=60, angle_min=None)

	tor_lig_img_save_path = lig_img_save_path / f'butylbenzene-{a1}_{a2}_{a3}_{a4}'

	ligtor.make_torsion_img(t, show=False, save_path = str(tor_lig_img_save_path))
	ligtor.plot_dihedral_scatter(t, save_path = str(tor_lig_img_save_path)+"scatter", show=False)
	ligtor.plot_dihedral_histogram(t, angle_min=-180, save_path=str(tor_lig_img_save_path)+"hist", show=False)
