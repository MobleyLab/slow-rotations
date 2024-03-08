import parmed as pmd
import tempfile
SYMBOL_TO_ATMNUM = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
                        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'K': 19, 'Ar': 18,
                        'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Ni': 28, 'Co': 27,
                        'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
                        'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
                        'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
                        'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
                        'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
                        'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
                        'Bi': 83, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
                        'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
                        'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113,
                        'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}
ELEMENT_LIST = sorted([k for k in SYMBOL_TO_ATMNUM.keys()], key=len, reverse=True)


def get_atom_type(atype):
	for atm in ELEMENT_LIST:
		if atm.lower() in atype.lower():
			return atm

def make_standard_atmname(pmd_struct):
	''' Converts third column pdb atom name from any atom name
		to atom type naming convention
		In gmx to prevent forcefield issues, we rename ligand
		atom types either as lower case, or with extra characters
		to prevent clashes with protein and other forcefield 
		naming conventions
		That atom naming conventions are the same as the atomic
		symbol on the periodic table
	'''
	new_structure = pmd_struct

	first_residue = None
	for r in new_structure.residues:
		# Assumes only 1 residue intended in PDB file
		# Get the residue name of the first atom
		# Renames all atom residues to the same residue name
		# Fixes issue with UNK1 and UNK2 where the H's are a 
		# different residue name from the 
		if not first_residue:
			first_residue = r
		for a in new_structure.atoms:
			a.name = get_atom_type(a.name)
			a.atomic_number = SYMBOL_TO_ATMNUM[a.name]
			a.residue = first_residue

	return new_structure

def rename_lig_pdb_atoms(ilig_pdb: str, olig_pdb: str):
	lig_struct = pmd.load_file(ilig_pdb)
	nox_lig_struct = make_standard_atmname(lig_struct)

	pdb_correct_atype = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)

	nox_lig_struct.save(pdb_correct_atype.name, overwrite=True)

	with open(ilig_pdb, 'r') as i:
		with open(pdb_correct_atype.name, 'r') as intermed:
			with open(olig_pdb, 'w') as o:
				for line in intermed:
					if not line.startswith("END"):
						o.write(line)
				for line in i:
					if line.startswith('CONECT'):
						o.write(line)

