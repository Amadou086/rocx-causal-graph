# ==============================
# RoCX Streamlit App + Core
# ==============================
import streamlit as st
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
import random
import math

# -------------------------------
# Core functions
# -------------------------------

def parse_npz(npz_data):
    """Extract SMILES strings from an NPZ file object"""
    entries = npz_data["smiles"]
    smiles_list = []

    for row in entries:
        row = list(row)
        if len(row) >= 2:
            smiles_list.append(row[1])
    return smiles_list

def mol_to_nx(smiles):
    """Convert SMILES string to NetworkX graph"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), label=atom.GetSymbol())
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType().name
        )
    return G

def generate_ba_graph(n_nodes=20, m_edges=2):
    """Generate synthetic BA graph"""
    return nx.barabasi_albert_graph(n_nodes, m_edges)

def compute_causal_scores(G, model=None):
    """
    Placeholder: Compute node-level causal contributions
    Replace with actual GNN-NCM inference
    """
    return {node: random.random() for node in G.nodes()}

def apply_counterfactual_robustness(causal_scores, lambda_decay=1.0):
    """Apply robustness weighting to causal scores"""
    return {node: score * math.exp(-lambda_decay * random.random())
            for node, score in causal_scores.items()}

def extract_robust_subgraph(G, robust_scores, threshold=0.5):
    """Extract robust explanatory subgraph (Gamma_r)"""
    selected_nodes = [n for n, s in robust_scores.items() if s >= threshold]
    return G.subgraph(selected_nodes).copy()

def plot_graph(G, node_scores=None):
    """Plot a graph with optional node coloring by scores"""
    fig, ax = plt.subplots(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)

    if node_scores:
        # normalize scores to [0,1] for coloring
        scores = np.array([node_scores.get(n,0) for n in G.nodes()])
        if scores.max() - scores.min() < 1e-6:
            scores_norm = np.ones_like(scores) * 0.5
        else:
            scores_norm = (scores - scores.min()) / (scores.max() - scores.min())
        node_colors = [plt.cm.Reds(s) for s in scores_norm]
    else:
        node_colors = "lightblue"

    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=1200,
        node_color=node_colors,
        edge_color="gray",
        font_size=8,
        ax=ax
    )
    return fig

def build_causal_graph(dataset_type="molecular", uploaded_file=None,
                       n_nodes=20, m_edges=2, lambda_decay=1.0, threshold=0.5):
    """
    Build a causal graph with counterfactual reasoning
    Supports:
    - Molecular dataset (.npz)
    - Synthetic BA graph
    Returns: robust explanatory subgraph Gamma_r and node scores
    """
    if dataset_type == "molecular":
        if uploaded_file is None:
            raise ValueError("No NPZ file uploaded.")
        npz_data = np.load(uploaded_file, allow_pickle=True)
        smiles_list = parse_npz(npz_data)
        smiles = random.choice(smiles_list)
        G = mol_to_nx(smiles)
        if G is None:
            raise ValueError("Failed to generate graph from SMILES.")
    else:
        # Synthetic BA graph
        G = generate_ba_graph(n_nodes=n_nodes, m_edges=m_edges)

    # 1. Apply GNN-NCM (placeholder)
    causal_scores = compute_causal_scores(G)

    # 2. Counterfactual robustness
    robust_scores = apply_counterfactual_robustness(causal_scores, lambda_decay)

    # 3. Extract robust subgraph
    Gamma_r = extract_robust_subgraph(G, robust_scores, threshold)

    return Gamma_r, robust_scores

# ==============================
# Streamlit App
# ==============================
st.set_page_config(page_title="RoCX Causal Graph", layout="wide")
st.title("RoCX Causal Graph Explorer")

# Dataset selection
dataset_type = st.radio("Select dataset type:", ["Synthetic BA Graph", "Molecular Dataset"])

# Initialize graph
G = None
node_scores = None

if dataset_type == "Synthetic BA Graph":
    n_nodes = st.number_input("Number of nodes", min_value=5, max_value=100, value=20)
    m_edges = st.number_input("Edges per new node", min_value=1, max_value=5, value=2)

    # Build graph
    G, node_scores = build_causal_graph(dataset_type="ba", n_nodes=n_nodes, m_edges=m_edges)

elif dataset_type == "Molecular Dataset":
    uploaded_file = st.file_uploader("Upload NPZ molecular dataset", type=["npz"])
    if uploaded_file:
        G, node_scores = build_causal_graph(dataset_type="molecular", uploaded_file=uploaded_file)
    else:
        st.info("Please upload a .npz file to generate a molecular causal graph.")

# Plot graph
if G:
    fig = plot_graph(G, node_scores=node_scores)
    st.pyplot(fig)
