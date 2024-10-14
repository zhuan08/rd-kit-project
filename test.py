from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# adj = ["red", "big", "tasty"]
# fruits = ["apple", "cherry", "tomato"]
# ripeness = ["ripe", "unripe"]
# for x in adj:
#     for y in fruits:
#         for z in ripeness:
#             print(x, y, z)

#
# lst1 = []
# mol = []
# smiles = []
# for i in lst1:
#     k = Chem.MolFromSmiles(i, True)
#     j = Chem.AddHs(k)
#     mol.append(j)
# for a in mol:
#     b = Chem.MolToSmiles(a)
#     smiles.append(b)
# print(*smiles, sep=', \n')

ligand_top = Chem.MolFromSmiles('SC1=CN=C(C2=N1)N3C(N2C)=[Ir]4=C5N(C)C6=NC(S)=CN=C6N5C7=C(Cl)C(I)=C(Cl)C3=C74')
ligand_bot = Chem.MolFromSmiles('CN1C(N2C3=NC(S)=CN=C31)=[Ir]4=C5N(C)C6=NC=C(S)N=C6N5C7=C(Cl)C(I)=C(Cl)C2=C74')
Draw.MolsToGridImage((Chem.AddHs(ligand_top), Chem.AddHs(ligand_bot)),
                     legends=['(a)', '(b)'], molsPerRow=1, subImgSize=(200,200)).show()
