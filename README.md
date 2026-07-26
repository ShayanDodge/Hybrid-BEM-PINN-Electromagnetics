# A Novel Hybrid Boundary Element—Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics

[![Project Website](https://img.shields.io/badge/website-BEM_PINN-blue)](https://shayandodge.github.io/)
[![License](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

---

## 🔗 Official Implementation

This repository is the **official implementation** of the paper:

> **“A Novel Hybrid Boundary Element—Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics”**  
> Published in IEEE Access, 2024  
> https://ieeexplore.ieee.org/document/10755077  

It provides a hybrid computational framework that combines the **Boundary Element Method (BEM)** and **Physics-Informed Neural Networks (PINNs)** for solving electromagnetic problems governed by the **Laplace equation**.

The repository reproduces and extends the methodology proposed in the paper, including domain decomposition and coupling via interface boundary conditions.

---

## 🧠 Overview

This work introduces a **hybrid BEM–PINN approach** for numerical solutions in electromagnetics:

- **PINN** is used to solve the PDE in part of the domain  
- **BEM** is applied to the remaining subdomain  
- The two methods are **coupled through interface conditions**  

This approach leverages:
- the flexibility of neural networks  
- the efficiency of boundary-based numerical methods  

to solve partial differential equations in complex domains.

---

## 🔑 Keywords

Boundary Element Method (BEM), Physics-Informed Neural Networks (PINNs), Electromagnetics, Laplace Equation, Scientific Machine Learning, Domain Decomposition

---

## 📌 Versions

<table>

<tr>
<td width="55%" valign="top">

### 🔹 v1.0.0 — Homogeneous L-Shaped Domain  
**Folder:** `Homogeneous L-Shaped Domain (V1.0.0)`

- Baseline implementation of the hybrid BEM–PINN framework  
- Solves the Laplace equation on a homogeneous L-shaped domain  
- Domain decomposition:
  - Left subdomain → solved using PINN  
  - Right-bottom subdomain → solved using BEM  
- Coupling enforced through interface boundary conditions  

👉 This version validates the methodology presented in the paper.

</td>

<td width="45%" align="center">

<img src="https://github.com/user-attachments/assets/0f59587c-7119-4b7e-a790-7432307fa024" width="100%"/>

<sub>BEM–PINN coupling on L-shaped domain</sub>

</td>
</tr>

<tr>
<td width="55%" valign="top">

### 🔹 v2.0.0 — Non-Homogeneous Nested Square Domains *(Coming Soon)*

- Extension to non-homogeneous / multi-material domains  
- Incorporates interface conditions between different materials  
- Improved physical modeling capabilities  

🚧 Currently under development.

</td>

<td width="45%" align="center">

<img width="5805" height="5597" alt="V2 0 0" src="https://github.com/user-attachments/assets/3e56e54f-17b1-43ab-b8e1-fca9e0cb03fa" />

<sub>BEM–PINN coupling on Nested Square Domains</sub>

</td>
</tr>

</table>

---

## 🛠 Requirements

- MATLAB  
- Deep Learning Toolbox  
- Optimization Toolbox  
- Parallel Computing Toolbox  

---

## 📊 Output

Each version produces:

- Scalar potential field  
- Gradient fields  
- BEM reconstruction  
- Combined contour visualization  

---

## 📄 Paper Reference

The scientific foundation of this repository is described in:

> **Barmada, S., Dodge, S., Tucci, M., Formisano, A., Di Barba, P., & Mognaschi, M. E.**  
> *“A Novel Hybrid Boundary Element—Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics.”*  
> IEEE Access, 2024  
> https://ieeexplore.ieee.org/document/10755077  

📌 This repository is the official codebase accompanying the above publication.

---

## 📚 Citation

If you use this work, please cite:

```bibtex
@article{barmada2024hybrid,
  title={A Novel Hybrid Boundary Element—Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics},
  author={Barmada, Sami and Dodge, Shayan and Tucci, Mauro and Formisano, Alessandro and Di Barba, Paolo and Mognaschi, Maria Evelina},
  journal={IEEE Access},
  volume={12},
  pages={171444--171457},
  year={2024},
  publisher={IEEE}
}


