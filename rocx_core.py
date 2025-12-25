import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
import random


# -------------------------------
# Parse NPZ dataset
# -------------------------------
def parse_npz(npz_data):
    """
    Extract SMILES strings from an NPZ file object
    """
    entries = npz_data["smiles"]
    smiles_list = []

    for row in entries:
        row = list(row)
        if len(row) >= 2:
            smiles_list.append(row[1])

    return smiles_list


# -------------------------------
# Convert SMILES → NetworkX graph
# -------------------------------
def mol_to_nx(smiles):
    """
    Convert a SMILES string to a NetworkX graph
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            label=atom.GetSymbol()
        )

    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType().name
        )

    return G


# -------------------------------
# Main function used by Streamlit
# -------------------------------
def build_causal_graph(uploaded_file):
    """
    Build ONE causal graph from an uploaded NPZ file
    """
    # Load uploaded file object (no local path!)
    npz_data = np.load(uploaded_file, allow_pickle=True)

    smiles_list = parse_npz(npz_data)

    if not smiles_list:
        raise ValueError("No valid SMILES found in NPZ file.")

    # Pick one molecule at random
    smiles = random.choice(smiles_list)

    G = mol_to_nx(smiles)

    if G is None:
        raise ValueError("Failed to generate graph from SMILES.")

    return G


# -------------------------------
# Plot graph
# -------------------------------
def plot_graph(G):
    """
    Plot the causal graph using NetworkX + Matplotlib
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    pos = nx.spring_layout(G, seed=42)

    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=1200,
        node_color="lightblue",
        edge_color="gray",
        font_size=8,
        ax=ax
    )

    return fig
