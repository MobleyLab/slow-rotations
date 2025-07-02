import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path
import numpy as np

from scipy.special import kl_div
base_path = Path('/Users/megosato/Desktop/slow-rotations/simulations/hipen/unrestrained')

ligname = 'lig_163'

# bound is inside the protein or in this case, restrained
topf = base_path / f'nvt_4344392.gro'
trajf = base_path / f'nvt_4344392.xtc'
ligcode = "UNK"

smiles = "[H]O[C@@]1([H])c2c(c(Oc3c([H])c(F)c([H])c(C#N)c3[H])c([H])c([H])c2S(=O)(=O)C([H])([H])[H])C([H])([H])C1(F)F"

# notes CL should be Cl NOT CL
# remove x from atom names of gmx topology
# atom numbering cannot have no space next to atom type

ligtor = tor.LigandTorsionFinder(str(trajf),str(topf),ligcode,smiles)

lig_img_save_path = base_path / "images_4344392"
lig_img_save_path.mkdir(exist_ok=True)

for t in ligtor.get_torsions():

	print('TORSION', t)

	a1,a2,a3,a4 = tuple(t)

	am, sa = ligtor.shift_torsion_angles(t)
	#np.savetxt(f'/Users/megosato/Desktop/{ligname}-{a1}_{a2}_{a3}_{a4}.dat', sa)


	#X, score, angle_min = ligtor.get_kde(tor, num_bins=60, angle_min=None)

	tor_lig_img_save_path = lig_img_save_path / f'butylbenzene-{a1}_{a2}_{a3}_{a4}'

	ligtor.make_torsion_img(t, show=False, save_path = str(tor_lig_img_save_path))
	ligtor.plot_dihedral_scatter(t, angle_min=0, save_path = str(tor_lig_img_save_path)+"scatter", show=False)
	ligtor.plot_dihedral_histogram(t, angle_min=0, save_path=str(tor_lig_img_save_path)+"hist", show=False)
