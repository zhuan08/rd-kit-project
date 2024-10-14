from functools import reduce
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdmolops import CombineMols, SanitizeMol, Kekulize, FragmentOnBonds, GetMolFrags, RemoveHs
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolFragmentToSmiles

def adjust_iridium_bonds(mol, target_elem = "Ir", chemdraw = False):
    '''Given an organometallic molecule downloaded from the Cambridge
    Structural Database, modify the bonds with iridium to match conventions and
    return'''
    editable_mol = RWMol(mol)
    
    # If you don't make a list, it loops infinitely over the bonds it's creating
    for bond in list(editable_mol.GetBonds()):
        iridium = None
        nitrogen = None
        carbene = None
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() in ["N", "P"] and \
                bond.GetEndAtom().GetFormalCharge() == 1:
            iridium = bond.GetBeginAtom()
            nitrogen = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() in ["N", "P"] and \
                bond.GetBeginAtom().GetFormalCharge() == 1:
            iridium = bond.GetEndAtom()
            nitrogen = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() == "C" and \
                bond.GetEndAtom().GetTotalValence() == 3:
            iridium = bond.GetBeginAtom()
            carbene = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() == "C" and \
                bond.GetBeginAtom().GetTotalValence() == 3:
            iridium = bond.GetEndAtom()
            carbene = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            
        if not chemdraw:
            if nitrogen is not None:
                # Replace N+ - Ir with N -> Ir
                nitrogen.SetFormalCharge(0)
                
            if iridium is not None and nitrogen is not None:
                editable_mol.RemoveBond(start_idx, end_idx)
                editable_mol.AddBond(start_idx, end_idx, Chem.rdchem.BondType.DATIVE)

        if iridium is not None and carbene is not None:
            editable_mol.RemoveBond(start_idx, end_idx)
            editable_mol.AddBond(start_idx, end_idx, Chem.rdchem.BondType.DOUBLE)
            
    outmol = editable_mol.GetMol()
    Chem.SanitizeMol(outmol)
    
    return outmol

def ligate(ligands, metal_atom_element = "Ir", metal_atom = None):
    '''Given a list of RDKit molecule objects representing organometallic
    molecules, return a molecule in which all are attached to a shared metal
    atom'''
    for ligand in ligands:
        ligand.RemoveAllConformers()
        Kekulize(ligand)
    if metal_atom is None:
        metal_atom = MolFromSmiles(f"[{metal_atom_element}]")
    # Create a molecule that contains the metal atom as well as all the
    # ligands, but without bonds between the metal and the ligands
    mol = reduce(CombineMols, ligands, metal_atom)
    # Record the index of the metal atom, so that bonds can be created to it
    # Based on testing, I assume this is always the first atom
    metal_atom_index = 0
    # Create an editable molecule and begin batch editing, which makes
    # functions available for adding and removing bonds as well as removing the
    # dummy atoms, and should keep the indexes stable as we do so
    editable = Chem.EditableMol(mol)
    editable.BeginBatchEdit()
    # Check each bond to see if it is a coordination site, and if so, bond the
    # atom to the metal. Coordination sites are recognized as bonds to the
    # target element, which should not include the added atom since it does not
    # have any bonds
    for bond in mol.GetBonds():
        # Figure out which, if any, of the atoms in the bond is the coordinating atom
        if bond.GetBeginAtom().GetSymbol() == metal_atom_element:
            dummy_atom_index = bond.GetBeginAtomIdx()
            coordination_atom_index = bond.GetEndAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == metal_atom_element:
            dummy_atom_index = bond.GetEndAtomIdx()
            coordination_atom_index = bond.GetBeginAtomIdx()
        else:
            continue
        # Add a bond to metal of the same type as the bond to the dummy atom
        # Has to have metal second for the dative bonds to work
        # Order matters for dative bonds
        # https://www.rdkit.org/docs/RDKit_Book.html#dative-bonds
        editable.AddBond(coordination_atom_index, metal_atom_index, bond.GetBondType())
        # Remove the bond between the coordinating atom and the dummy atom
        editable.RemoveBond(dummy_atom_index, coordination_atom_index)
    # Remove all the metals from the original molecules
    for i, atom in enumerate(editable.GetMol().GetAtoms()):
        if i > 0 and atom.GetSymbol() == metal_atom_element:
            editable.RemoveAtom(i)
    # Apply the changes and get the molecule
    editable.CommitBatchEdit()
    outmol = editable.GetMol()
    # Shouldn't have any conformers but just in case
    outmol.RemoveAllConformers()
    # Probably already done by GetMol but just in case
    SanitizeMol(outmol)
    return outmol

def fragment(mol, metal_atom_element = 'Ir'):
    '''Given an RDKit moelcule object representing an organometallic molecule, return a list of smiles strings of ligands, each bonded to a metal atom'''
    # Index of the metal atom
    metal_atom_index = [i for i, atom in enumerate(mol.GetAtoms()) \
                        if atom.GetSymbol() == metal_atom_element][0]

    # Bonds with the metal atom
    metal_bond_indices = [i for i, bond in enumerate(mol.GetBonds()) \
            if bond.GetBeginAtom().GetSymbol() == metal_atom_element or \
               bond.GetEndAtom().GetSymbol() == metal_atom_element]

    # Next line can give "ValueError: empty bond indices"
    fragmented_mol = FragmentOnBonds(mol, metal_bond_indices, addDummies = False)

    # The indices of the fragments, which includes the ligands, as well as the
    # metal atom
    fragment_index_tuples = GetMolFrags(fragmented_mol)

    # List of ligands to return
    ligand_smiles_strings = []

    for fragment_indices in fragment_index_tuples:
        # Since the output ligand should include the metal, add its index
        ligand_indices = [metal_atom_index] + list(fragment_indices)

        # If it's just the metal index, skip
        if len(set(ligand_indices)) == 1:
            continue

        # Create the ligand by extracting the targeted atoms
        ligand_smiles_string = MolFragmentToSmiles(mol, atomsToUse =  ligand_indices)

        ligand_smiles_strings.append(ligand_smiles_string)

    return ligand_smiles_strings

