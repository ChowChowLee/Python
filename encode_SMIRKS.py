
from __future__ import annotations
from collections.abc import Iterable
from collections import defaultdict
from typing import Optional
from typing import Any
from pathlib import Path

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from rdkit.Chem.Draw import IPythonConsole

def predict_product(molecule, reaction_smirks):
    """
    Predicts the product(s) of a reaction given a molecule and a SMIRKS string.

    Args:
        molecule (rdkit.Chem.Mol): The reactant molecule.
        reaction_smirks (str): The SMIRKS notation representing the reaction.

    Returns:
        list of str: List of SMILES strings representing the predicted product(s).
    """
    rxn = AllChem.ReactionFromSmarts(reaction_smirks)
    product = rxn.RunReactants([Chem.MolFromSmiles(x) for x in ('NCCc1ccccc1','C1CC1C(=O)')])
    final_smiles_list = []
    try:
        for i in range(len(product)):
            final_smiles = Chem.MolToSmiles(product[i][0])
            final_smiles_list.append(final_smiles)
    except:
        print("We got bug!")
    return final_smiles_list, product


def main() -> None:

    new_file = Path.cwd()
    print(f"Current working directory: {new_file}")
    print((new_file.exists()))
    full_path = new_file.resolve() # Switch from the relative path into abs path
    print(f"Full path: {full_path}")
    print(f"Parent: {full_path.parent}")
    print(f"Stem: {full_path.stem}")
    print(f"Is directory: {full_path.is_dir()}")
    print((f"Is file: {full_path.is_file()}"))
    # new_dir = new_file.parent / "new_folder_1"
    # new_dir.mkdir()
    # new_dir.rmdir()



    rxn = AllChem.ReactionFromSmarts(
        '[cH1:1]1:[c:2](-[CH2:7]-[CH2:8]-[NH2:9]):[c:3]:[c:4]:[c:5]:[c:6]:1.[#6:11]-[CH1;R0:10]=[OD1]>>[c:1]12:[c:2](-[CH2:7]-[CH2:8]-[NH1:9]-[C:10]-2(-[#6:11])):[c:3]:[c:4]:[c:5]:[c:6]:1')

    reacts = [Chem.MolFromSmiles(x) for x in ('NCCc1ccccc1', 'C1CC1C(=O)')]
    ps = rxn.RunReactants(reacts)
    ps0 = ps[0]
    for p in ps0:
        Chem.SanitizeMol(p)
    Draw.MolsToGridImage(ps0)

    # ethanol = Chem.MolFromSmiles('CCO')
    # isopropanol = Chem.MolFromSmiles('CC(C)O')
    # phenol = Chem.MolFromSmiles('c1cccc(O)c1')
    # t_butanol = Chem.MolFromSmiles('C(C)(C)(C)O')
    # glycerol = Chem.MolFromSmiles('OCC(O)CO')
    #
    # alcohols = [ethanol, isopropanol, phenol, t_butanol, glycerol]
    #
    # Draw.MolsToGridImage(alcohols, molsPerRow=5)


if __name__ == "__main__":
    main()