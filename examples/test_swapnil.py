import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path

from scipy.special import kl_div

base_path = Path('/Users/megosato/Desktop/simulations/bace_a')


# bound is inside the protein or in this case, restrained
topf_bnd = base_path / f'confout.gro'
trajf_bnd = base_path / f'traj_comp.xtc'


# notes CL should be Cl NOT CL
# remove x from atom names of gmx topology
# atom numbering cannot have no space next to atom type

smiles = "[N+](=C1N[C@](C(=O)N1C)(C2=CC=CC=C2)C3=CC(=CC=C3)C4=CN=CC=C4)([H])[H]"
system="bace_a"

ligcode2 = "LIG"
ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd), ligcode2, smiles)


imgs_save_path = Path('/Users/megosato/Desktop/images_bace_a/')
imgs_save_path.mkdir(exist_ok=True)

lig_img_save_path = imgs_save_path 
lig_img_save_path.mkdir(exist_ok=True)


for t in ligtor_bnd.get_torsions():

	a1,a2,a3,a4 = tuple(t)

	tor_lig_img_save_path = lig_img_save_path / f'{system}-{a1}_{a2}_{a3}_{a4}'

	ligtor_bnd.make_torsion_img(t, save_path = str(tor_lig_img_save_path))

	sys_t = [ligtor_bnd.convert_ligidx_to_sysidx(i) for i in t]


	ligtor_bnd.plot_dihedral_scatter(sys_t, save_path = str(tor_lig_img_save_path)+"scatter", show=False)

