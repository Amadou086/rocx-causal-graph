import streamlit as st
import numpy as np
from rocx_core import (
    build_causal_graph,       # builds graph from SMILES or BA graph
    apply_gnn_ncm,            # runs GNN-NCM to get causal contributions
    compute_counterfactuals,  # computes counterfactual robustness
    extract_robust_subgraph,  # returns Gamma_r
    plot_graph
)
import networkx as nx

st.set_page_config(page_title="RoCX Causal Graph Explorer", layout="wide")
st.title("RoCX Causal Graph Explorer")

# ----------------------
# Dataset Selection
# ----------------------
dataset_type = st.radio("Select dataset type:", ["Synthetic BA Graph", "Molecular Dataset"])

# ----------------------
# Load / Generate Graph
# ----------------------
G = None
if dataset_type == "Synthetic BA Graph":
    n_nodes = st.number_input("Number of nodes", min_value=5, max_value=100, value=20)
    m_edges = st.number_input("Edges per new node", min_value=1, max_value=5, value=2)
    # generate BA graph
    G = nx.barabasi_albert_graph(n_nodes, m_edges)
    st.success(f"Synthetic BA graph generated with {n_nodes} nodes.")

elif dataset_type == "Molecular Dataset":
    uploaded_file = st.file_uploader("Upload NPZ molecular dataset", type=["npz"])
    if uploaded_file:
        # build molecular graph from SMILES
        G = build_causal_graph(uploaded_file)
        st.success("Molecular graph constructed.")
    else:
        st.info("Please upload a .npz file to generate a molecular causal graph.")

# ----------------------
# Apply RoCX Causal + Counterfactual Reasoning
# ----------------------
if G:
    # 1. Apply GNN Neural Causal Model
    causal_scores = apply_gnn_ncm(G)

    # 2. Compute counterfactual robustness
    robustness_scores = compute_counterfactuals(G, causal_scores)

    # 3. Extract robust explanatory subgraph (Gamma_r)
    Gamma_r = extract_robust_subgraph(G, causal_scores, robustness_scores)

    # 4. Plot the robust causal subgraph
    fig = plot_graph(Gamma_r)
    st.pyplot(fig)
