from rdkit import Chem
from rdkit.Chem import rdFMCS
from openeye import oechem

def get_atom_by_index(mol, idx):
	''' returns atom object with index for the given mol
		works for both rdkit and openeye (ducktyped)
	'''
	for a in mol.GetAtoms():
		if a.GetIdx() == idx:
			return a

def convert_dihedral(mapping, mol1_dih):
	return [ mapping[i] for i in mol1_dih ]

def rd_map_mols(rdmol1, rdmol2):
	mcs = rdFMCS.FindMCS([rdmol1, rdmol2])
	patt = Chem.MolFromSmarts(mcs.smartsString)
	rdmol1_match = rdmol1.GetSubstructMatch(patt)
	rdmol2_match = rdmol2.GetSubstructMatch(patt)

	mapping = dict()
	for m1,m2 in zip(rdmol1_match, rdmol2_match):
		mapping[m1] = m2
	return mapping

def oe_map_mols(oemol1, oemol2):
	# MO TODO: create way to map molecules using openeye MCS
	# see rmsd_calculator.py from sampl containers
	raise utils.NotImplementedError

def map_mols(mol1, mol2):
	''' returns a dictionary of mapped atoms for the SAME ligand
		requires rdkit mols
	'''
	# MO TODO: check number of atoms
	# MO TODO: check number of each atom type
	# MO TODO: ensure they are the SAME ligand

	rd_type = Chem.rdchem.Mol
	oe_type = oechem.OEMol

	if type(mol1) == rd_type and type(mol2) == rd_type:
		mapping = rd_map_mols(mol1, mol2)

	elif type(mol1) == oe_type and type(mol2) == oe_type:
		mapping = oe_map_mols(mol1, mol2)

	else: 
		raise utils.NotImplementedError

