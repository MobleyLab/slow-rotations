# rdkit wrapper

from slow_rotations.utils import NotImplementedError
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS
from openff.toolkit import Molecule
import numpy as np

def load_rdmol_from_file(molfile: str, removeHs=False):
    """
    load molecule in molfile as an RDKit Mol
    
    Args:
        molfile: str; sdf or mol2 filename of small molecule
        removeHs: bool; strip Hs from the rdkit molecule

    Returns:
        Chem.Mol: rdkit molecule 
    """

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
    """
    Assigns the bond order of a rdkit molecule lacking proper bond information based on its smiles
    
    Args:
        smiles: str; smiles str representing the correct bond orders
        molfile: str; pdb, sdf, or mol2 filename of molecule that lacks correct bond order information

    Returns:
        Chem.Mol: RDKit Mol with 
    """

    # the smiles will contain the bond order information
    # the molfile contains the positions and index order
    # the returned mol will contain the bond order from the 
    # smiles but the index order from the molfile
    lig_mol_wo_bond_orders = load_rdmol_from_file(molfile)
    Chem.SanitizeMol(lig_mol_wo_bond_orders, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)


    smi_mol = Chem.MolFromSmiles(smiles)

    if smi_mol.GetNumAtoms() < lig_mol_wo_bond_orders.GetNumAtoms():
        smi_mol = Chem.AddHs(smi_mol)
        AllChem.EmbedMolecule(smi_mol)
        AllChem.UFFOptimizeMolecule(smi_mol)
        Chem.MolToMolFile(smi_mol, "smi_mol.mol")
        sanitize_rdmol(smi_mol)

        Chem.MolToMolFile(lig_mol_wo_bond_orders, "pdb_mol.mol")
        sanitize_rdmol(lig_mol_wo_bond_orders)


    lig_mol = AllChem.AssignBondOrdersFromTemplate(smi_mol, lig_mol_wo_bond_orders)
    Chem.MolToMolFile(lig_mol_wo_bond_orders, "bondorder_mol.mol")
    return lig_mol

def sanitize_rdmol(mol) :
    """
    Checks and corrects the chemical structure of molecule. (Checks valence, aromaticity, kekulization, and stereochemistry)
    
    Args:
        mol: Chem.Mol; rdkit molecule

    Returns:
        Chem.Mol: santized rdkit molecule
    """
    # From Jeffrey Wagner, sanitize RDKIT mol
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL)
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    return mol

def get_mapped_heavy_atom_indices(mol, mapped_atoms):
    """
    Finds RDKit Atom objects of interest based on mapped_atoms indices list

    Args:
        mol: (Chem.Mol) RDKit Mol
        mapped_atoms: [int]; indices of atoms of interest

    Returns:
        list: list of RDKit Atoms 
    """
    mapped_heavy_atoms = []
    for ai in mapped_atoms:
        if mol.GetAtomWithIdx(ai).GetAtomicNum() > 1:
            mapped_heavy_atoms.append(ai)
    return mapped_heavy_atoms
    
    
def get_mapped_bonds(mol, mapped_atoms):
    """
    Finds RDKit bonds objects of interest based on mapped_atoms indices list

    Args:
        mol: (Chem.Mol) RDKit Mol
        mapped_atoms: [int]; indices of atoms of interest

    Returns:
        list: list of RDKit bonds
    """
    hit_bonds = []
    for bond in mol.GetBonds():
        aid1 = bond.GetBeginAtomIdx()
        aid2 = bond.GetEndAtomIdx()
        if aid1 in mapped_atoms and aid2 in mapped_atoms:
            hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
    return hit_bonds

def get_rotatable_bonds(mol):
    """
    Finds the rotatable bonds in molecule. A rotatable bond is a single, non-ring bond between two nonterminal, non–triple-bonded atoms. 

    Args:
        mol: (Chem.Mol) RDKit Mol

    Returns:
        list: list of RDKit Mol bonds that are rotatable
    """
    rotbond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    rotbonds = mol.GetSubstructMatches(rotbond)
    bonds = []
    for i1,i2 in rotbonds:
        bonds.append(mol.GetBondBetweenAtoms(i1,i2))
    return list(bonds)


def highlight_dihedral(mol, dihedral, save_path=None):  
    """
    Creates an image of the molecule with the dihedral of interest highlighted in red
    
    Args:
        mol: Chem.Mol; rdkit molecule
        dihedral: (int, int, int, int); 4 consecutive atom indices that represent the torsion of interest 
        save_path: str; file path to save image to

    Returns:
        None
    """

    mol = Chem.Mol(mol)
    mol_wo_H = Chem.RemoveHs(mol)

    mcs = rdFMCS.FindMCS([mol, mol_wo_H])
    patt = Chem.MolFromSmarts(mcs.smartsString)

    query_match = mol.GetSubstructMatch(patt)
    template_match = mol_wo_H.GetSubstructMatch(patt)

    index_convert = {query_match[i]: template_match[i] for i in range(len(query_match))}
    index_convert_rev = {template_match[i]: query_match[i] for i in range(len(query_match))}
    print(index_convert)
    print(index_convert_rev)

    new_dihedral = list()

    for aid in dihedral:
        try:
            new_dihedral.append(index_convert[aid])

        except KeyError:
            pass

    AllChem.Compute2DCoords(mol_wo_H)

    #for atm in mol_wo_H.GetAtoms():
    #    atm.SetProp("atomNote", str(index_convert_rev[atm.GetIdx()]))
        
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



















