# rdkit wrapper

from utils import NotImplementedError
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit import Molecule

def load_rdmol_from_file(molfile: str):
	''' molfile: sdf file or mol2 file representing bond order
				 in the trajectory

		Returns: rdkit mol
	'''
	# loads in molecule as openff Molecule object
	# returns molecule as an rdkit mol object
	if molfile.endswith(".sdf") or molfile.endswith(".mol2"):
		offmol = Molecule(molfile)
		rdmol = offmol.to_rdkit()

	elif  molfile.endswith(".pdb"):
		rdmol = Chem.MolFromPDBFile(molfile, removeHs=False)

	else:
		raise NotImplementedError

	return rdmol

def assign_bond_order_from_smiles(smiles: str, molfile: str):
	# the smiles will contain the bond order information
	# the molfile contains the positions and index order
	# the returned mol will contain the bond order from the 
	# smiles but the index order from the molfile

	# need to check to make sure I do actually need this
	lig_mol_wo_bond_orders = load_rdmol_from_file(molfile)
	smi_mol = Chem.MolFromSmiles(smiles)
	smi_mol = Chem.AddHs(smi_mol)

	lig_mol = AllChem.AssignBondOrdersFromTemplate(smi_mol, lig_mol_wo_bond_orders)
	return lig_mol

def sanitize_rdmol(mol):
	# From Jeffrey Wagner, sanitize RDKIT mol
	Chem.SanitizeMol(mol, Chem.SANITIZE_ALL)
	Chem.AssignStereochemistryFrom3D(mol)
	Chem.Kekulize(mol, clearAromaticFlags=True)
	Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
	return mol

	



	
