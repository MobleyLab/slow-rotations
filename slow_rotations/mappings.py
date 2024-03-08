from rdkit import Chem
from rdkit.Chem import rdFMCS
from openeye import oechem, oeomega
from openff.toolkit import Molecule
import numpy as np


# do we need to account for symmetry

def get_atom_by_index(mol, idx):
	''' returns atom object with index for the given mol
		works for both rdkit and openeye (ducktyped)
	'''
	for a in mol.GetAtoms():
		if a.GetIdx() == idx:
			return a

def convert_dihedral(mapping, mol1_dih):
	''' takes a mapping dictionary and returns the dihedral mapping
		for the molecule based on the mapping
		dictionary in the form of {mol1_idx: mol2_idx}
	'''
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
	# oemol1 = pattern
	# opemol2 = target

	# has trouble with charged ligand

	atomexpr = oechem.OEExprOpts_AtomicNumber
	bondexpr = 0

	if oemol1.NumAtoms() != oemol2.NumAtoms():
		raise ValueError(f"Molecules not same, {oemol1.NumAtoms()} v. {oemol2.NumAtoms()}")

	ss = oechem.OESubSearch(oemol1, atomexpr, bondexpr)
	oechem.OEPrepareSearch(oemol2, ss)
	   
	for count, match in enumerate(ss.Match(oemol2)):
	    oemol1_match = np.array([ma.pattern.GetIdx() for ma in match.GetAtoms()], dtype=int)
	    oemol2_match = np.array([ma.target.GetIdx() for ma in match.GetAtoms()], dtype=int)

	    mapping = dict()
	    
	    for m1,m2 in zip(oemol1_match, oemol2_match):
	    	mapping[m1] = m2
	    return mapping


def reindex_smiles_from_pdb(smiles, pdb):
	''' does the "opposite" of assign bond order from smiles
	'''

	smimol = oechem.OEMol()
	oechem.OESmilesToMol(smimol, smiles)
	omega = oeomega.OEOmega()
	omega.SetMaxConfs(100)
	omega.SetStrictStereo(False)
	omega(smimol)

	pdbmol = Molecule.from_file(pdb, allow_undefined_stereo=True).to_openeye()

	mapping = oe_map_mols(smimol, pdbmol)

	rd_smimol = Molecule.from_openeye(smimol, allow_undefined_stereo=True).to_rdkit()
	smiles_match = np.zeros(len(mapping))
	pdb_match = np.zeros(len(mapping))

	for i,kv in enumerate(mapping.items()):
		smiles_match[i] = kv[0]
		pdb_match[i] = kv[1]

	pdb_argsort = np.argsort(pdb_match)
	smiles_reorder = [ int(i) for i in smiles_match[pdb_argsort] ]

	rd_smimol_reorder = Chem.RenumberAtoms(rd_smimol, list(smiles_reorder))

	return rd_smimol_reorder


def map_mols(mol1, mol2):
	''' returns a dictionary of mapped atoms for the SAME ligand
		requires rdkit mols (for now)
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

	return mapping


def check_symmetry(mol, torsion):
    oechem.OEPerceiveSymmetry(mol)
    cb_indices = [1,2]
    for i,cb_aidx in enumerate(cb_indices):
        cb_atm = get_atom_by_index(mol, torsion[cb_aidx])

        nbr_ct = 0
        symmetry_numbers = []

        for nbr_atm in cb_atm.GetAtoms():
            nbr_ct += 1
            
            if nbr_atm.GetIdx() == torsion[cb_indices[(i+1)% 2]]:
                # it is the other atom in the central bond
                # skip this atom
                continue
            symmetry_numbers.append(nbr_atm.GetSymmetryClass())
            
        if nbr_ct == 4 and len(set(symmetry_numbers)) == 1:
            return True

        elif nbr_ct == 3 and len(set(symmetry_numbers)) == 1:
            return True

    return False 
