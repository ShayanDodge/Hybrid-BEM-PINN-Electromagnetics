# A Novel Hybrid Boundary Element – Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics

[![Project Website](https://img.shields.io/badge/website-BEM_PINN-blue)](https://shayandodge.github.io/)
[![License](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

This repository contains implementations of a coupled **Boundary Element Method (BEM)** and **Physics-Informed Neural Network (PINN)** framework for solving the **Laplace equation**.

The scientific foundations of this work are presented in:

> **Barmada, S., Dodge, S., Tucci, M., Formisano, A., Di Barba, P., & Mognaschi, M. E.**  
> *“A Novel Hybrid Boundary Element—Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics.”*  
> IEEE Access, 2024, [10.1109/ACCESS.2024.3500039](https://ieeexplore.ieee.org/document/10755077).  
> [Read the article](https://ieeexplore.ieee.org/document/10755077)

This repository implements, reproduces, and extends the methodology described in the above peer-reviewed publication. If you use this code, dataset, or results in your research, please cite the original article.

---

## 📌 Versions

<table>

<tr>
<td width="55%" valign="top">

### 🔹 v1.0.0 — Homogeneous L-Shaped Domain  
**Folder:** `Homogeneous L-Shaped Domain (V1.0.0)`

- Baseline implementation of the BEM-PINN framework  
- Solves the Laplace equation on a homogeneous L-shaped domain  
- Domain decomposition:
  - Left subdomain → solved using PINN  
  - Right-bottom subdomain → solved using BEM  
- Coupling enforced through interface boundary conditions  

👉 This version establishes the core methodology and validates the approach.

</td>

<td width="45%" align="center">

<img src="https://github.com/user-attachments/assets/0f59587c-7119-4b7e-a790-7432307fa024" width="100%"/>

<sub>BEM-PINN coupling on L-shaped domain</sub>

</td>
</tr>

<tr>
<td width="55%" valign="top">

### 🔹 v1.1.0 — Non-Homogeneous L-Shaped Domain *(Coming Soon)*

- Extension to non-homogeneous / multi-material domains  
- Incorporates interface conditions between different materials  
- Improved physical realism and modeling capability  

🚧 This version is currently under development and will be released soon.

</td>

<td width="45%" align="center">

<div style="font-size:80px;">🚧</div>

<sub>Next version in development</sub>

</td>
</tr>

</table>

## 🛠 Requirements

- MATLAB  
- Deep Learning Toolbox  
- Optimization Toolbox  
- Parallel Computing Toolbox  

## 📊 Output

Each version produces:

- Scalar potential field  
- Gradient fields  
- BEM reconstruction  
- Combined contour visualization  

## 📄 Related Paper

The scientific foundations of this work are presented in:

> **Barmada, S., Dodge, S., Tucci, M., Formisano, A., Di Barba, P., & Mognaschi, M. E.**  
> *“A Novel Hybrid Boundary Element—Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics.”*  
> IEEE Access, 2024, [10.1109/ACCESS.2024.3500039](https://ieeexplore.ieee.org/document/10755077).  
> [Read the article](https://ieeexplore.ieee.org/document/10755077)

📌 This repository implements, reproduces, and extends the methodology described in the above peer-reviewed publication. If you use this code, dataset, or results in your research, please cite the original article.

## 📚 Citation

If you use this work, please cite:

🔗 Paper link: https://ieeexplore.ieee.org/document/10755077

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
