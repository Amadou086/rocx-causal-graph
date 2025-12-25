import streamlit as st
from rocx_core import build_causal_graph, plot_graph

st.set_page_config(page_title="RoCX Causal Graph", layout="wide")

st.title("RoCX Causal Graph Explorer")

uploaded_file = st.file_uploader("Upload NPZ dataset", type=["npz"])

if uploaded_file:
    G = build_causal_graph(uploaded_file)
    fig = plot_graph(G)
    st.pyplot(fig)
else:
    st.info("Please upload a .npz file to generate a causal graph.")
