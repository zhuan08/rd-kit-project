from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Left ligand for Ir complex
l_ligand = Chem.MolFromSmiles('CN(C(N=CC=N1)C1N2C3=CC(C(F)(F)F)=C4)C2=[Ir](C3=C4N5C6C7N=CC=N6)=C5N7C')
Draw.MolToImage(l_ligand).show() #displays left ligand
