import rdkit.Chem as Chem
from rdkit.Chem import AllChem


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
    product = rxn.RunReactants([molecule, ])
    final_smiles_list = []
    try:
        for i in range(len(product)):
            final_smiles = Chem.MolToSmiles(product[i][0])
            final_smiles_list.append(final_smiles)
    except:
        pass
    return final_smiles_list


# Example usage:
reactant_smiles = "CC(=O)C"  # Replace with your reactant SMILES
smirks = "[*:1]-[C:2]=[C:3]-[N;H1:4]>>[*:1]-[C:2]=[C:3]-[C:5]=[N:4]"  # Example SMIRKS
reactant_molecule = Chem.MolFromSmiles(reactant_smiles)
predicted_products = predict_product(reactant_molecule, smirks)

print("Predicted product(s):")
for product_smiles in predicted_products:
    print(product_smiles)
