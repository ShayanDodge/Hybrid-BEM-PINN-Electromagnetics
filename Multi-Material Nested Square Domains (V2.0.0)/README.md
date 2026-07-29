# Multi-Material Nested Square Domains

This version extends the original **BEM–PINN framework** from a homogeneous domain to **heterogeneous (multi-material) nested square domains**. The computational domain is decomposed into multiple subdomains with different material properties, while continuity conditions are enforced across the shared interfaces.

## Overview

<table>
<tr>
<td valign="top" width="60%">

### Domain decomposition

The computational domain consists of nested square regions with different material properties.

- **Inner domain:** solved using a Physics-Informed Neural Network (PINN)
- **Outer domain:** solved using the Boundary Element Method (BEM)
- **Material interface:** couples both solvers through the exchange of Dirichlet and Neumann boundary data

The hybrid formulation reconstructs the complete solution by iteratively enforcing interface continuity between the BEM and PINN subdomains.

</td>

<td valign="top" width="40%" align="center">

<img src="https://github.com/user-attachments/assets/96a60dab-eda9-4b0c-910c-31f88d6c74a6"
alt="Nested square multi-material domain"
width="320">

</td>

</tr>
</table>

---

# Main Idea

This implementation performs the following steps:

1. Generates a nested square multi-material geometry
2. Defines interface and external boundary conditions
3. Trains the PINN in the interior material region
4. Computes interface potential and flux using automatic differentiation
5. Transfers interface quantities to the BEM solver
6. Solves the surrounding material domain using BEM
7. Enforces continuity of potential and normal flux across the interface
8. Reconstructs and visualizes the coupled multi-material solution

---

# Repository Structure

```text
multi-material-nested-square-bem-pinn/
├── README.md
├── LICENSE
├── .gitignore
├── main/
│   └── bem_pinn_nested_square.mlx
├── src/
│   ├── model.mlx
│   ├── modelLoss.mlx
│   ├── objectiveFunction.mlx
│   ├── bem_interface_solver.m
│   ├── coordin.m
│   ├── intfundamsquare.m
│   ├── intfundamsquare_post.m
│   ├── intnormalsquare.m
│   ├── intnormalsquare_post.m
│   ├── normalder.m
│   ├── initializeHe.mlx
│   ├── initializeZeros.mlx
│   ├── parameterStructToVector.mlx
│   └── parameterVectorToStruct.mlx
├── examples/
│   ├── homogeneous_case/
│   └── multi_material_case/
└── results/
    ├── data.mat
    ├── data_bem.mat
    └── figures/
```

---

# New Features in Version 2

Compared with **v1.0.0**, this version introduces:

- Nested square multi-material geometries
- Multiple material properties
- Interface continuity constraints
- Coupled transmission conditions
- Improved modular BEM–PINN architecture
- Cleaner code organization
- Enhanced visualization and validation tools

---

# Requirements

This code requires:

- MATLAB
- Deep Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox

---

# How to Run

1. Open the repository in MATLAB.
2. Open `bem_pinn_nested_square.mlx`.
3. Run the Live Script.
4. The workflow will:

   - generate the nested square geometry
   - assign material properties
   - train the PINN
   - solve the BEM region
   - couple both solvers across the interface
   - reconstruct the complete solution
   - generate publication-quality figures

---

# Output

The implementation produces:

- scalar potential distribution
- electric (or magnetic) field components
- interface potential
- interface normal flux
- coupled BEM–PINN solution
- comparison with FEM reference solutions
- convergence history
- error analysis

<p align="center">
<img src="YOUR_RESULTS_IMAGE.png" width="900">
</p>

---

# Scientific Contribution

Version 2 extends the original hybrid BEM–PINN methodology by introducing heterogeneous material regions separated by interfaces. The framework enforces the continuity of both scalar potential and normal flux, enabling accurate simulation of transmission problems in multi-material domains while preserving the computational advantages of combining Boundary Element Methods with Physics-Informed Neural Networks.

---

# Related Paper

This repository extends the methodology presented in:

> **Barmada, S., Dodge, S., Tucci, M., Formisano, A., Di Barba, P., & Mognaschi, M. E.**  
> *A Novel Hybrid Boundary Element–Physics-Informed Neural Network Method for Numerical Solutions in Electromagnetics.*  
> IEEE Access, 2024.

Future publications describing the multi-material extension will be added here.

---

# Citation

If you use this repository, please cite the original IEEE Access paper together with the forthcoming publication describing the multi-material extension.

```bibtex
@article{barmada2024hybrid,
  title={A Novel Hybrid Boundary Element--Physics Informed Neural Network Method for Numerical Solutions in Electromagnetics},
  author={Barmada, Sami and Dodge, Shayan and Tucci, Mauro and Formisano, Alessandro and Di Barba, Paolo and Mognaschi, Maria Evelina},
  journal={IEEE Access},
  volume={12},
  pages={171444--171457},
  year={2024},
  publisher={IEEE}
}
```
