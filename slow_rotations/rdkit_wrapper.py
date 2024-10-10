# rdkit wrapper

from utils import NotImplementedError
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS
from openff.toolkit import Molecule
import numpy as np

def load_rdmol_from_file(molfile: str, removeHs=False):
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
        rdmol = Chem.RWMol(Chem.MolFromPDBFile(molfile, removeHs=removeHs, sanitize=False, proximityBonding=True))
        conf = rdmol.GetConformer()

        # Define typical bond distances (in Angstroms)
        bond_cutoff = {
            'C': 1.2,  # Typical C-H bond length
            'N': 1.2,  # Typical N-H bond length
            'O': 1.1,  # Typical O-H bond length
        }

        # Loop over atoms and infer missing bonds based on distances
        for atom in rdmol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Only deal with hydrogens
                pos_h = np.array(conf.GetAtomPosition(atom.GetIdx()))

                for neighbor in rdmol.GetAtoms():
                    if neighbor.GetIdx() == atom.GetIdx():  # Skip the same atom
                        continue

                    # Get the element symbol of the neighbor atom
                    neighbor_symbol = neighbor.GetSymbol()
                    if neighbor_symbol in bond_cutoff:
                        pos_neighbor = np.array(conf.GetAtomPosition(neighbor.GetIdx()))
                        distance = np.linalg.norm(pos_h - pos_neighbor)

                        # If the distance matches a typical bond length, create a bond
                        if distance < bond_cutoff[neighbor_symbol]:
                            try:    
                                rdmol.AddBond(atom.GetIdx(), neighbor.GetIdx(), Chem.BondType.SINGLE)
                            except RuntimeError:
                                pass

        rdmol = Chem.Mol(rdmol)

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

    if smi_mol.GetNumAtoms() < lig_mol_wo_bond_orders.GetNumAtoms():
        smi_mol = Chem.AddHs(smi_mol)
        Chem.MolToMolFile(smi_mol, "smi_mol.mol")
        sanitize_rdmol(smi_mol)

        Chem.MolToMolFile(lig_mol_wo_bond_orders, "pdb_mol.mol")
        sanitize_rdmol(lig_mol_wo_bond_orders)

        for i,a in enumerate(smi_mol.GetAtoms()):
            print(i, a.GetAtomicNum())

        print()

        for i,a in enumerate(lig_mol_wo_bond_orders.GetAtoms()):
            print(i, a.GetAtomicNum())

    lig_mol = AllChem.AssignBondOrdersFromTemplate(smi_mol, lig_mol_wo_bond_orders)
    return lig_mol

def sanitize_rdmol(mol):
    # From Jeffrey Wagner, sanitize RDKIT mol
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL)
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    return mol

def get_mapped_heavy_atom_indices(mol, mapped_atoms):
    mapped_heavy_atoms = []
    for ai in mapped_atoms:
        if mol.GetAtomWithIdx(ai).GetAtomicNum() > 1:
            mapped_heavy_atoms.append(ai)
    return mapped_heavy_atoms
    
    
def get_mapped_bonds(mol, mapped_atoms):
    hit_bonds = []
    for bond in mol.GetBonds():
        aid1 = bond.GetBeginAtomIdx()
        aid2 = bond.GetEndAtomIdx()
        if aid1 in mapped_atoms and aid2 in mapped_atoms:
            hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
    return hit_bonds


def highlight_dihedral(mol, dihedral, save_path=None):  
    mol = Chem.Mol(mol)
    mol_wo_H = Chem.RemoveHs(mol)

    mcs = rdFMCS.FindMCS([mol, mol_wo_H])
    patt = Chem.MolFromSmarts(mcs.smartsString)

    query_match = mol.GetSubstructMatch(patt)
    template_match = mol_wo_H.GetSubstructMatch(patt)

    index_convert = {query_match[i]: template_match[i] for i in range(len(query_match))}

    new_dihedral = list()

    for aid in dihedral:
        try:
            new_dihedral.append(index_convert[aid])

        except KeyError:
            pass

    AllChem.Compute2DCoords(mol_wo_H)
        
    highlightAtoms = get_mapped_heavy_atom_indices(mol_wo_H, new_dihedral)
    highlightBonds = get_mapped_bonds(mol_wo_H, new_dihedral)

    Draw.MolToImageFile(
        mol_wo_H,
        save_path,
        highlightBonds=highlightBonds,
        highlightAtoms=highlightAtoms,
        size=(1000, 500),
        wedgeBonds=True, kekulize=True, wedgeFontSize=0, wedgeLineWidth=0
    )



















