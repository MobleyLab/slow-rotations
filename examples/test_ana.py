import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare
from pathlib import Path

from scipy.special import kl_div


base_path = Path('/Users/megosato/Desktop/simulations/hif2a/complex/2_bad')


# bound is inside the protein or in this case, restrained
topf_bnd = base_path / f'../2/em.gro'
trajf_bnd = base_path / f'traj.xtc'


# notes CL should be Cl NOT CL
# remove x from atom names of gmx topology
# atom numbering cannot have no space next to atom type

smiles = "CS(=O)(=O)C1=CC=C(C2=C1C(C(C2)(F)F)N)OC3=CC(=CC(=C3)C#N)F"
#smiles = "C1=C(C(=C(C(=C1)S(=O)(=O)C(F)F)CO)Cl)OC2=CC(=CC(=C2)C#N)F"
system= "hif2a"

ligcode2 = "LGA"
ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd), ligcode2, smiles)


imgs_save_path = Path('/Users/megosato/Desktop/images_hif2a/2_bad/')
imgs_save_path.mkdir(exist_ok=True)

lig_img_save_path = imgs_save_path 
lig_img_save_path.mkdir(exist_ok=True)


for t in ligtor_bnd.get_torsions():

	a1,a2,a3,a4 = tuple(t)

	tor_lig_img_save_path = lig_img_save_path / f'{system}-{a1}_{a2}_{a3}_{a4}'

	ligtor_bnd.make_torsion_img(t, save_path = str(tor_lig_img_save_path))

	sys_t = [ligtor_bnd.convert_ligidx_to_sysidx(i) for i in t]

	ligtor_bnd.plot_dihedral_scatter(sys_t, save_path = str(tor_lig_img_save_path)+"scatter", show=False)

