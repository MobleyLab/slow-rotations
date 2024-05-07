# convert rdkit to oemol
# ensure atom numberings are the same

from openff.toolkit import Molecule

def get_oemol_from_rdmol(rdmol):
	''' converts an rdkit Molecule to an OEMol object 
		using openff toolkit
	'''
	off_rdmol = Molecule.from_rdkit(rdmol)
	return Molecule.to_openeye(off_rdmol)

def get_rdmol_from_oemol(oemol):
	off_oemol = Molecule.from_rdkit(oemol)
	return Molecule.to_openeye(off_oemol)
	

