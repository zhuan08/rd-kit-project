import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Load substitutions from a file
substitutions = pd.read_csv('substitutions.csv')
# Column 'smiles' as a list
substitution_smiles = substitutions['smiles'].tolist()
sub_1 = substitution_smiles
sub_2 = ["[CH3]", "[H]"]
sub_3 = substitution_smiles

ligand_l = Chem.MolFromSmiles('SC1=C(N=C(C2=N1)N3C(N2C)=[Ir]4=C5N(C)C6=NC(S)=C(N=C6N5C7=C(Cl)C(I)=C(Cl)C3=C74)S)S')
ligand_top = Chem.MolFromSmiles('SC1=CN=C(C2=N1)N3C(N2C)=[Ir]4=C5N(C)C6=NC(S)=CN=C6N5C7=C(Cl)C(I)=C(Cl)C3=C74')
ligand_bot = Chem.MolFromSmiles('CN1C(N2C3=NC(S)=CN=C31)=[Ir]4=C5N(C)C6=NC=C(S)N=C6N5C7=C(Cl)C(I)=C(Cl)C2=C74')


def smiles_to_rdkit_mol(sub_list):  # converts list of smiles strings to rdkit molecules
    rdkit_list = []                 # list of rdkit molecules
    for sub in sub_list:            # goes through lists and does conversion
        rdkit_mol = Chem.MolFromSmiles(sub)
        rdkit_list.append(rdkit_mol)
    return rdkit_list


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
                if ligand_w_sub3 is not None:
                    substituted_ligands.add(Chem.MolToSmiles(ligand_w_sub3, True))
    return substituted_ligands


def sub_one_group(ligand, sub, p1=Chem.MolFromSmiles('S')):
    substituted_ligands = set()
    for x in sub:
        ligand_w_sub = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
        if ligand_w_sub is not None:
            substituted_ligands.add(Chem.MolToSmiles(ligand_w_sub, True))
    return substituted_ligands


# Show images of sub_combinations
# # all combinations for ligand_l
# set_l = sub_combinations(ligand_l, sub_1, sub_2, sub_3)
# set_fl = set(smiles for smiles in set_l if Chem.MolFromSmiles(smiles) is not None)
# Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set_fl], molsPerRow=8, subImgSize=(200, 200)).show()

# all combinations for ligand_top
set_t = sub_combinations(ligand_top, sub_1, sub_2, sub_3)
set_ft = set(smiles for smiles in set_t if Chem.MolFromSmiles(smiles) is not None)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set_ft], molsPerRow=8, subImgSize=(200, 200)).show()

# all combinations for ligand_bot
set_b = sub_combinations(ligand_bot, sub_1, sub_2, sub_3)
set_fb = set(smiles for smiles in set_b if Chem.MolFromSmiles(smiles) is not None)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set_fb], molsPerRow=8, subImgSize=(200, 200)).show()

# Show images of sub_one_group
# sub_one_set = sub_one_group(ligand_bot, sub_1)
# sos_set = set(smiles for smiles in sub_one_set if Chem.MolFromSmiles(smiles) is not None)
# Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in sos_set], molsPerRow=5, subImgSize=(200, 200)).show()

# # Save SMILES strings of substituted ligands
# output_smiles_filename = 'substituted_ligands.txt'
# with open(output_smiles_filename, 'w') as f:
#     for x in set_fl:
#         f.write(x + '\n')
