# convert rdkit to oemol
# ensure atom numberings are the same

from openff.toolkit import Molecule

# def get_oemol_from_rdmol(rdmol):
# 	'''
#   	converts RDKit Mol to OpenEye OEMol

#     Args:
#         rdmol: Chem.Mol RDKit Mol

#     Returns:
#         OEMol: oechem.OEMol OpenEye Mol
# 	'''
# 	off_rdmol = Molecule.from_rdkit(rdmol)
# 	return Molecule.to_openeye(off_rdmol)

# def get_rdmol_from_oemol(oemol):
# 	'''
#   	converts RDKit Mol to OpenEye OEMol

#     Args:
#         rdmol: oechem.OEMol OpenEye Mol

#     Returns:
#         RDMol: Chem.Mol RDKit Mol
# 	'''
# 	off_oemol = Molecule.from_rdkit(oemol)
# 	return Molecule.to_openeye(off_oemol)
	

