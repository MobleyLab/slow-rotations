from rdkit import Chem
from rdkit.Chem import rdFMCS
#from openeye import oechem, oeomega
from openff.toolkit import Molecule
import numpy as np


# do we need to account for symmetry

def get_atom_by_index(mol, idx):
	''' 
  	finds RDKit Atom objects of interest based on provided atom

	Args:
		mol: (Chem.Mol) RDKit Mol
		idx: int; index of atom of interest

	Returns:
		Chem.Atom: RDKit Atom
	'''
	for a in mol.GetAtoms():
		if a.GetIdx() == idx:
			return a

def torsions_equivalent(t1, t2):
	'''
  	checks if 2 torsions are equivalent

	Args:
		t1: (int, int, int, int); indices of atoms in torsion of interest
		t2: (int, int, int, int); indices of atoms in torsion of interest

	Returns:
		bool
	'''
	t1_a1, t1_a2, t1_a3, t1_a4 = t1
	t2_a1, t2_a2, t2_a3, t2_a4 = t2

	if t1_a1 == t2_a1 and t1_a2 == t2_a2 and t1_a3 == t2_a3 and t1_a4 == t2_a4:
		return True

	if t1_a1 == t2_a4 and t1_a2 == t2_a3 and t1_a3 == t2_a2 and t1_a4 == t2_a1:
		return True
	return False

def convert_dihedral(mapping, dihedral):
	'''
  	converts indices in dihedral based on the provided mapping

	Args:
		mapping: dict; mapping dihedral atoms
		dihedral: 

	Returns:
		bool
	'''
	return [ mapping[i] for i in dihedral ]

def rd_map_mols(rdmol1, rdmol2):
	'''
	Given 2 RDKit Mols, determines the mapping based on the MCS

	Args:
		rdmol1: Chem.Mol
		rdmol2: Chem.Mol

	Returns:
		dict: mapping of atoms
	'''

	mcs = rdFMCS.FindMCS([rdmol1, rdmol2])

	patt = Chem.MolFromSmarts(mcs.smartsString)

	rdmol1_match = rdmol1.GetSubstructMatch(patt)
	rdmol2_match = rdmol2.GetSubstructMatch(patt)

	mapping = dict()
	for m1,m2 in zip(rdmol1_match, rdmol2_match):
		mapping[m1] = m2
	return mapping

# def oe_map_mols(oemol1, oemol2):
# 	# oemol1 = pattern
# 	# opemol2 = target

# 	# has trouble with charged ligand

# 	atomexpr = oechem.OEExprOpts_AtomicNumber
# 	bondexpr = 0

# 	if oemol1.NumAtoms() != oemol2.NumAtoms():
# 		raise ValueError(f"Molecules not same, {oemol1.NumAtoms()} v. {oemol2.NumAtoms()}")

# 	ss = oechem.OESubSearch(oemol1, atomexpr, bondexpr)
# 	oechem.OEPrepareSearch(oemol2, ss)
	   
# 	for count, match in enumerate(ss.Match(oemol2)):
# 		oemol1_match = np.array([ma.pattern.GetIdx() for ma in match.GetAtoms()], dtype=int)
# 		oemol2_match = np.array([ma.target.GetIdx() for ma in match.GetAtoms()], dtype=int)

# 		mapping = dict()
		
# 		for m1,m2 in zip(oemol1_match, oemol2_match):
# 			mapping[m1] = m2
# 		return mapping


# def reindex_smiles_from_pdb(smiles, pdb):
# 	''' 
# 	Assigns the indices to the smiles based on the atom indices of the small molecule in the pdb. 
# 	(Does the opposite of rdkit_wrapper.assign_bond_order_from_smiles)

# 	Args:
# 		smiles: str; smiles string representing small molecule in pdb
# 		pdb: str; pdbfile representing small molecule in smiles

# 	Returns:
# 		Chem.Mol
# 	'''

# 	smimol = oechem.OEMol()
# 	oechem.OESmilesToMol(smimol, smiles)
# 	omega = oeomega.OEOmega()
# 	omega.SetMaxConfs(100)
# 	omega.SetStrictStereo(False)
# 	omega(smimol)

# 	pdbmol = Molecule.from_file(pdb, allow_undefined_stereo=True).to_openeye()

# 	mapping = oe_map_mols(smimol, pdbmol)

# 	rd_smimol = Molecule.from_openeye(smimol, allow_undefined_stereo=True).to_rdkit()
# 	smiles_match = np.zeros(len(mapping))
# 	pdb_match = np.zeros(len(mapping))

# 	for i,kv in enumerate(mapping.items()):
# 		smiles_match[i] = kv[0]
# 		pdb_match[i] = kv[1]

# 	pdb_argsort = np.argsort(pdb_match)
# 	smiles_reorder = [ int(i) for i in smiles_match[pdb_argsort] ]

# 	rd_smimol_reorder = Chem.RenumberAtoms(rd_smimol, list(smiles_reorder))

# 	return rd_smimol_reorder


def map_mols(mol1, mol2):
	''' 
	Creates dictionary mapping of indices between atoms of mol1 and mol2 

	Args:
		mol1: RDKit Mol
		mol2: RDKit Mol

	Returns:
		dict: {int: int}; dictionary mapping of indices {mol1: mol2}
	'''

	rd_type = Chem.rdchem.Mol

	if type(mol1) == rd_type and type(mol2) == rd_type:
		mapping = rd_map_mols(mol1, mol2)

	# elif type(mol1) == oe_type and type(mol2) == oe_type:
	# 	mapping = oe_map_mols(mol1, mol2)

	else: 
		raise utils.NotImplementedError

	return mapping


def check_symmetry(mol, torsion):
	''' 
	Checks symmetry of central bond atoms in a torsion using RDKit.

	Args:
		mol: RDKit Mol
		torsion: (int, int, int, int) - atom indices defining the torsion

	Returns:
		bool: True if torsion has symmetry, False otherwise
	'''
	
	# Compute symmetry classes
	symmetry_classes = Chem.CanonicalRankAtoms(mol, breakTies=False)

	# The central bond atoms (b and c in a-b-c-d torsion)
	cb_indices = [1, 2]

	for i, cb_aidx in enumerate(cb_indices):
		cb_atm = mol.GetAtomWithIdx(torsion[cb_aidx])

		nbr_ct = 0
		symmetry_numbers = []

		for nbr_atm in cb_atm.GetNeighbors():
			nbr_idx = nbr_atm.GetIdx()
			nbr_ct += 1

			if nbr_idx == torsion[cb_indices[(i+1) % 2]]:
				# skip the other atom in the central bond
				continue

			symmetry_numbers.append(symmetry_classes[nbr_idx])

		# Check for symmetry condition
		if nbr_ct == 4 and len(set(symmetry_numbers)) == 1:
			return True
		elif nbr_ct == 3 and len(set(symmetry_numbers)) == 1:
			return True

	return False

# def oe_check_symmetry(mol, torsion):
#	 ''' 
# 	Creates dictionary mapping of indices between atoms of mol1 and mol2 

#	 Args:
#		 mol1: RDKit Mol
#		 torsion: (int, int, int, int)

#	 Returns:
#		 bool: True if torsion has symmetry, False otherwise
#	 '''
#	 oechem.OEPerceiveSymmetry(mol)
#	 cb_indices = [1,2]
#	 for i,cb_aidx in enumerate(cb_indices):
#		 cb_atm = get_atom_by_index(mol, torsion[cb_aidx])

#		 nbr_ct = 0
#		 symmetry_numbers = []

#		 for nbr_atm in cb_atm.GetAtoms():
#			 nbr_ct += 1
			
#			 if nbr_atm.GetIdx() == torsion[cb_indices[(i+1)% 2]]:
#				 # it is the other atom in the central bond
#				 # skip this atom
#				 continue
#			 symmetry_numbers.append(nbr_atm.GetSymmetryClass())
			
#		 if nbr_ct == 4 and len(set(symmetry_numbers)) == 1:
#			 return True

#		 elif nbr_ct == 3 and len(set(symmetry_numbers)) == 1:
#			 return True

#	 return False 
