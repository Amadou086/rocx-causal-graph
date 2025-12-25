<p align="center">
  <img src="Figures/CRL.png" width="300">
</p>


# RoCX — Integrating Causality and Counterfactual Reasoning for Robust and Explainable Graph Neural Networks.

**RoCX** is a framework that combines **causality, counterfactual reasoning, and structural analysis** to produce robust and interpretable explanations for graph neural network predictions on a wide range of graph types.

This repository accompanies the **RoCX Master’s Thesis** and includes:
- Core causal graph construction logic
- Molecular graph parsing from SMILES
- Interactive visualization via Streamlit
- A deployable causal explanation demo

---

## 🎓 Academic Context

This project is developed as part of a **Master’s thesis in Data Science & Artificial Intelligence (ECE Paris)**.

**Thesis title**  
> *RoCX: Integrating Causality and Counterfactual Reasoning for Robust and Explainable Graph Neural Networks*

The work is inspired by recent advances in:
- Causal inference (Pearl, 2018)
- Graph-based explainability
- Molecular representation learning

---

## 🔍 Motivation

Modern graph-based models (e.g., GNNs) often provide **high predictive accuracy** but **lack causal interpretability**, especially in molecular domains.

RoCX addresses this gap by:

- Moving beyond correlation-based explanations
- Explicitly modeling **cause–effect relations** between molecular components
- Providing **human-interpretable causal subgraphs**
- Enabling **interactive inspection** of causal structures
