# test script for slow_rotations
import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import rdkit.Chem
from rdkit.Chem import rdFMCS
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis.dihedrals import Dihedral

#trajf = "/Users/megosato/Desktop/pseudo_binding/pseudo_trajectories/61095_6_10_20_16.xtc"
#topf = "/Users/megosato/Desktop/pseudo_binding/pseudo_trajectories/61095_6_10_20_16.gro"
topf = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_61095.gro'
trajf = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_61095.xtc'

smiles = "CCOc1ccc2nc(/N=C\c3ccccc3O)sc2c1"
molfile = "../TD_Project/Gas_Phase_Sim/Ligand_mol2_files/ligand_61095.mol2"


lig_tor_find = tor.LigandTorsionFinder(trajf,topf,"UNK", smiles)

for a in lig_tor_find.rdmol.GetAtoms():
	print(a.GetAtomicNum())

torsions = lig_tor_find.get_torsions()
print(torsions)
print(torsions[0])

angles = lig_tor_find.get_torsion_angles([6,9,19,15])
print(angles)



u = lig_tor_find.mda_universe
ags = u.select_atoms("resname UNK")
for a in ags:
	print(lig_tor_find.convert_ligidx_to_sysidx(a.index))



# smiles_mol = rdkit.Chem.MolFromSmiles(smiles)
# smiles_mol = rdkit.Chem.AddHs(smiles_mol)


# mol2_mol = rdkit.Chem.MolFromMol2File(molfile)

# assigned_mol = rdw.assign_bond_order_from_smiles(smiles, molfile)
# assigned_mol = rdw.sanitize_rdmol(mol2_mol)

# for a, b, c in zip(smiles_mol.GetAtoms(), mol2_mol.GetAtoms(), assigned_mol.GetAtoms()):
# 	print(f"{a.GetAtomicNum()}\t{b.GetAtomicNum()}\t{c.GetAtomicNum()}")


# u = mda.Universe(topf, trajf)
# print(u)
# atms = u.select_atoms("all")
# for r,j in zip(lig_tor_find.rdmol.GetAtoms(),atms):
# 	print(f"{r.GetAtomicNum()}\t{j.name}")

# ags = [mda.AtomGroup(u.atoms[([2,7,13,17])])]
# Run = Dihedral(ags).run()
# shape = (Run.results.angles.shape)
# angles = Run.results.angles

# print(angles)

# atms = u.select_atoms("all")
# atms.write("/Users/megosato/Desktop/test.pdb", frames=u.trajectory[[0,]])


# rdmol_smiles = rdkit.Chem.MolFromSmiles(smiles)
# rdmol_smiles = rdkit.Chem.AddHs(rdmol_smiles)
# rdmol_mol2 = rdw.load_rdmol_from_file(molfile)

# i = 0
# for a, b in zip(rdmol_smiles.GetAtoms(), rdmol_mol2.GetAtoms()):
# 	print(a.GetAtomicNum(), b.GetAtomicNum())
# 	i += 1
# print(i)


# mcs = rdFMCS.FindMCS([rdmol_smiles, rdmol_mol2])

# print(mcs.numAtoms)

# patt = rdkit.Chem.MolFromSmarts(mcs.smartsString)
# query_match = rdmol_smiles.GetSubstructMatch(patt)
# template_match = rdmol_mol2.GetSubstructMatch(patt)

# print("query:", query_match)
# print("templ:", template_match)

# def get_atom_by_index(mol, idx):
# 	for a in mol.GetAtoms():
# 		if a.GetIdx() == idx:
# 			return a

# print(type(rdmol_smiles))

# for q,t in zip(query_match,template_match):
# 	print(get_atom_by_index(rdmol_smiles,q).GetAtomicNum(), get_atom_by_index(rdmol_mol2,t).GetAtomicNum())

# rdmol_smiles = rdw.sanitize_rdmol(rdmol_smiles)
# oe_smiles = mc.get_oemol_from_rdmol(rdmol_smiles)
# rdmol_mol2 = rdw.sanitize_rdmol(rdmol_mol2)
# oe_mol2 = mc.get_oemol_from_rdmol(rdmol_mol2)
# for q,t in zip(query_match,template_match):
# 	print(get_atom_by_index(oe_smiles, q).GetAtomicNum(), get_atom_by_index(oe_mol2, t).GetAtomicNum())

# print(type(oe_smiles))
