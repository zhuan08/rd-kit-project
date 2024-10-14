from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd
import re

# data = pd.read_csv("ligands_rdkit.csv", sep=",")    # reads data
#
# mol_id = []                 # adds mol_id from csv file to a new list
# k = 0
# for li in data['mol_id']:
#     if not pd.isna(li):     # filters for empty entries
#         mol_id.append(li)
#         k += 1
#
# mol_lig_a = []              # adds ligand a from csv file to new list as rdkit molecules
# i = 0
# for li in data['ligand_a']:
#     if not pd.isna(li):
#         mol_lig_a.append(Chem.MolFromSmiles(li))
#         i += 1
#     # print(*mol_lig_a, sep='\n')
#     # Draw.MolsToGridImage(mol_lig_a).show()

# smiles = 'C1=CC(=CC=C1C(CC(=O)O)CN)Cl'
# x = re.search("CC", smiles)
# if x:
#     print("Found CC")
# else:
#     print("Not found")


# img = Draw.MolToImage(Chem.MolFromSmiles('C1=CC(=CC=C1C(CC(=O)O)CN)Cl')).show()


# deleting structures (DeleteSubstructs module)
# m = Chem.MolFromSmiles('CCOC')
# patt = Chem.MolFromSmarts('C')
# rm = AllChem.DeleteSubstructs(m, patt)
# Draw.MolsToGridImage((m, patt, rm), subImgSize=(250, 250)).show()
#
# Chem.MolToSmiles(rm)

# # replacing substructures (ReplaceSubstructs module)
# repl = Chem.MolFromSmiles('OC')                         # replaced with
# patt = Chem.MolFromSmarts('[$(NC(=O))]')                # pattern to replace
# m = Chem.MolFromSmiles('CC(=O)N')                       # base structure
# rms = AllChem.ReplaceSubstructs(m, patt, repl)
# Draw.MolToImage(rms[0]).show()
#
# # removing side chains (ReplaceSidechains Module)
# m1 = Chem.MolFromSmiles('BrCCc1cncnc1C(=O)O')           # base structure
# core = Chem.MolFromSmiles('c1cncnc1')                   # keep core/ring
# tmp = Chem.ReplaceSidechains(m1, core)            # removes side chains
# Chem.MolToSmiles(tmp)
# Draw.MolsToGridImage((m1, core, tmp), subImgSize=(250, 250)).show()
#
# # removing cores (ReplaceCore module)
# tmp = Chem.ReplaceCore(m1, core)                  # removes core, keeps side chains
# Chem.MolToSmiles(tmp)
# Draw.MolToImage(tmp).show()
#
# # side chains are labeled based on the order they are found
# m1 = Chem.MolFromSmiles('c1c(CCO)ncnc1C(=O)O')          # CCO will be found first, then C(=O)O
# tmp = Chem.ReplaceCore(m1, core, labelByIndex=True)
# Chem.MolToSmiles(tmp)
# Draw.MolToImage(tmp).show()
#
# # splitting side chains into separate molecules (GetMolFrags module)
# rs = Chem.GetMolFrags(tmp, asMols=True)
# print(len(rs))
# Chem.MolToSmiles(rs[0])
# Draw.MolToImage(rs[0]).show()
# Chem.MolToSmiles(rs[1])
# Draw.MolToImage(rs[1]).show()


def smiles_to_rdkit_mol(sub_list):  # converts list of smiles strings to rdkit molecules
    rdkit_list = []                 # list of rdkit molecules
    for sub in sub_list:            # goes through lists and does conversion
        rdkit_mol = Chem.MolFromSmiles(sub)
        rdkit_list.append(rdkit_mol)
    return rdkit_list


ligand_top = Chem.MolFromSmiles('SC1=CN=C(C2=N1)N3C(N2C)[Ir]4C5N(C)C6=NC(S)=CN=C6N5C7=CC(C(F)(F)F)=CC3=C74')
ligand_bot = Chem.MolFromSmiles('SC1=CN=C2C(N3C(N2C)[Ir]4C5N(C)C6=NC=C(N=C6N5C7=CC(C(F)(F)F)=CC3=C74)S)=N1')
sub1 = ["[H]", "[CH3]", "[CH0](F)(F)F", "[CH0]#N", "[CH2]C(C)(C)C", "[SiH0](-C)(-C)-C", "[NH0]1C=CC=C1",
        "[cH0]1ccccc1", "[cH0]1c(C)cc(C)cc(C)1"]


def sub_one_group(ligand, sub, p1=Chem.MolFromSmiles('S')):
    substituted_ligands = set()
    for x in sub:
        ligand_w_sub = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
        if ligand_w_sub is not None:
            substituted_ligands.add(Chem.MolToSmiles(ligand_w_sub, True))
    return substituted_ligands


set1_unfiltered = sub_one_group(ligand_top, sub1)
set1 = set(smiles for smiles in set1_unfiltered if Chem.MolFromSmiles(smiles) is not None)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set1], molsPerRow=8, subImgSize=(200, 200)).show()
