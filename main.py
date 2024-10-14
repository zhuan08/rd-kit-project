from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

ligand_l = Chem.MolFromSmiles('SC1=C(N=C(C2=N1)N3C(N2C)=[Ir]4=C5N(C)C6=NC(S)=C(N=C6N5C7=C(Cl)C(I)=C(Cl)C3=C74)S)S')
# Draw.MolToImage(Chem.AddHs(ligand_l)).show()
# ligand_r = Chem.MolFromSmiles('[Ir]24C1=NC(=NC=C1OC3=CC(=CC(=[N]->23)C5=NC(=N[N]45)C(F)(F)F)C(C)(C)C)C(F)(F)F')

# -------- script for substituting CF3 group with list of substitutions (Zelin's Molecules)--------


def smiles_to_rdkit_mol(sub_list):      # converts list of smiles strings to rdkit molecules
    rdkit_list = []                     # list of rdkit molecules
    for sub in sub_list:                # goes through lists and does conversion
        rdkit_mol = Chem.MolFromSmiles(sub)
        rdkit_list.append(rdkit_mol)
    return rdkit_list


# list of substitutions
substitutions = ['[CH3]', '[CH0](C)(C)C', '[CH0](F)(F)F', '[CH0](C(F)(F)F)(C(F)(F)F)C(F)(F)F', '[NH0]1C=CC=C1',
                 '[SiH0](C)(C)C', '[H]']


def sub_combinations(ligand, sub1, sub2, sub3,              # returns a set of all combinations of substitutions
                     p1=Chem.MolFromSmiles('S'),            # limitation, 3 lists, 3 groups to substitute
                     p2=Chem.MolFromSmiles('Cl'),
                     p3=Chem.MolFromSmiles('I')):
    substituted_ligands = set()
    for x in sub1:
        for y in sub2:
            for z in sub3:
                ligand_w_sub1 = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
                ligand_w_sub2 = AllChem.ReplaceSubstructs(ligand_w_sub1, p2, Chem.MolFromSmiles(y), True)[0]
                ligand_w_sub3 = AllChem.ReplaceSubstructs(ligand_w_sub2, p3, Chem.MolFromSmiles(z), True)[0]
                substituted_ligands.add(Chem.MolToSmiles(ligand_w_sub3, True))
    return substituted_ligands


sub_1 = ["[NH2]", "[OH]"]
sub_2 = ["[CH3]", "[H]"]
sub_3 = ["[NH2]", "[OH]"]

set1 = sub_combinations(ligand_l, sub_1, sub_2, sub_3)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set1], molsPerRow=4, subImgSize=(200, 200)).show()
