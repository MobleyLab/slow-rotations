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
    'z4ylu_0A': 'CCC(=O)Nc1ccc(cc1)N(Cc2ccsc2)C(=O)Cn3c4ccccc4nn3 z4ylu',
    'P2203_0A': 'CN(C=1C=CC=2C=NC=C(NC(=O)CC=3C=CC=C(Cl)C3)C2C1)S(=O)(=O)C P2203',
    'x10876_0A': 'CN(C)C1=CC=C(C=C1)N(CC1=CSC=C1)C(=O)CN1N=NC2=C1C=CC=C2 x10876',
    'x10871_0A': 'CCC(=O)NC=1C=CC(=CC1)N(CC=2C=CSC2)C(=O)CN3N=NC=4C=CC=CC34 x10871',
    'x10870_0A': 'O=C(CN1N=NC=2C=CC=CC12)N(CC=3C=CSC3)C=4C=CC(NC(=O)C5CC5)=CC4 x10870',
    "x11798_0A": "CNC(=O)NC=1C=CC(=CC1)N(CC=2C=CSC2)C(=O)CN3N=NC=4C=CC=CC34",
    "x11790_0A": "CN(C)C=1C=CC(=CC1)N(CC=2C=CSC2)C(=O)CC=3C=NC=C4C=CC=CC34",
    "x10889_0A": "CS(=O)(=O)NCCC(C(=O)NC=1C=CC=NC1)C=2C=CC=CC2",
    "x3305_0A": "CN1C=C(CNC(=O)N(CCC=2C=CC=CC2)CC=3C=CC=NC3)N=N1",
    "x3305_1A": "CN1C=C(CNC(=O)N(CCC=2C=CC=CC2)CC=3C=CC=NC3)N=N1",
    "P0045_0A": "CC=1C=CN=CC1NC(=O)CC=2C=C(Cl)C=C(OC3CC(=O)N3)C2",
    "x10322_0A": "CC=1C=C(NC(=O)CC=2C=NC=NC2)C=C(OC3CC(=O)N3)C1",
    "x10371_0A": "O=C(CC=1C=CC=NC1)NC=2C=CC=C(OC3CC(=O)N3)C2",
    "x10387_0A": "FC=1C=C(NC(=O)CC=2C=NC=NC2)C=C(OC3CC(=O)N3)C1",
    "x10756_0A": "CC=1C=CN=CC1NC(=O)CC=2C=CC=C(OC3CC(=O)N3)C2",
    "x10789_0A": "CC=1C=CN=CC1NC(=O)CC=2C=C(Cl)C=C(OC3CC(=O)N3)C2",
    "x11641_0A": "CC(C(=O)NC=1C=NC=CC1C)C=2C=C(Cl)C=C(OC3CC(=O)N3)C2"
}

for ligname in mpro_ligs:

	system = f'{ligname}_0A'
	if ligname == 'z4ylu':
		system += "_HIE167"


	# bound is inside the protein or in this case, restrained
	topf_bnd = base_path / f'complex/{system}/npt.gro'
	trajf_bnd = base_path / f'complex/{system}/prod.xtc'


	smiles = mpro_ligs[ligname]

	# notes CL should be Cl NOT CL
	# remove x from atom names of gmx topology
	# atom numbering cannot have no space next to atom type

	ligcode2 = "UNK"
	ligtor_bnd = tor.LigandTorsionFinder(str(trajf_bnd),str(topf_bnd),ligcode2,smiles)


	imgs_save_path = Path('/Users/megosato/Desktop/images_bnd_only/')
	imgs_save_path.mkdir(exist_ok=True)

	lig_img_save_path = imgs_save_path / ligname
	lig_img_save_path.mkdir(exist_ok=True)


	for t in ligtor_bnd.get_torsions():


		a1,a2,a3,a4 = tuple(t)

		tor_lig_img_save_path = lig_img_save_path / f'{system}-{a1}_{a2}_{a3}_{a4}'

		sys_t = [ligtor_bnd.convert_ligidx_to_sysidx(i) for i in t]

		ligtor_bnd.make_torsion_img(t, save_path = str(tor_lig_img_save_path), show=False)