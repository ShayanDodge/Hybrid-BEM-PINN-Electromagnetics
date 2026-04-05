
# Homogeneous L-Shaped Domain: Combined BEM-PINN Solver for Laplace Equation

This repository contains the first implementation of a coupled **Boundary Element Method (BEM)** and **Physics-Informed Neural Network (PINN)** framework for solving the **Laplace equation** on a **homogeneous L-shaped domain**.

## Overview

The computational domain is divided into two subdomains:

- **Left subdomain:** solved using a PINN
- **Right-bottom subdomain:** solved using BEM

The two parts are coupled through interface and boundary information to reconstruct the full solution over the L-shaped domain.

## Main Idea

This code performs the following steps:

1. Constructs a 2D L-shaped domain from a rectangular point grid
2. Extracts Dirichlet, Neumann, and interface boundary points
3. Trains a PINN model on one subdomain
4. Uses the PINN-predicted interface information as input to the BEM solver
5. Solves the complementary subdomain with BEM
6. Reconstructs and visualizes the combined scalar field

## Repository Contents

- `BEM_PINN_Laplace_LD.mlx` — main live script
- `model.mlx` — neural network model definition
- `modelLoss.mlx` — loss computation
- `objectiveFunction.mlx` — optimization objective
- `bemforpinn_func_opt.m` — BEM-PINN coupling-related function
- `coordin.m` — boundary node generation for BEM
- `intfundamsquare.m` — fundamental solution integration
- `intfundamsquare_post.m` — post-processing fundamental integration
- `intnormalsquare.m` — normal derivative integration
- `intnormalsquare_post.m` — post-processing normal derivative integration
- `normalder.m` — normal derivative utility
- `initializeHe.mlx` — He initialization
- `initializeZeros.mlx` — zero initialization
- `parameterStructToVector.mlx` — parameter structure to vector conversion
- `parameterVectorToStruct.mlx` — parameter vector to structure conversion

    homogeneous-l-shaped-domain-bem-pinn/
        ├── README.md
        ├── LICENSE
        ├── .gitignore
        ├── main/
        │   └── run_bem_pinn_laplace.m
        ├── src/
        │   ├── model.mlx
        │   ├── modelLoss.mlx
        │   ├── objectiveFunction.mlx
        │   ├── bemforpinn_func_opt.m
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
        ├── results/
        │   ├── data.mat
        │   └── data_bem.mat
        ├── figures/
        │   └── domain_and_boundary_conditions.png
        ├── docs/
        │   └── notes.md
        └── examples/
            └── basic_run.md
  
## Requirements

This code requires:

- MATLAB
- Deep Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox

## How to Run

1. Open this repository in MATLAB.
2. Open `BEM_PINN_Laplace_LD.mlx`.
3. Run the live script.
4. The script will:
   - generate the L-shaped domain
   - train the PINN
   - solve the BEM subdomain
   - save output files such as `data.mat` and `data_bem.mat`
   - plot the combined contour field

## Output

The implementation produces:

- predicted scalar potential field
- predicted gradients
- BEM field reconstruction
- combined contour visualization of the BEM-PINN solution

## Notes

- This repository represents the **first version** of the code for the homogeneous L-shaped domain case.
- The implementation is intended for research and reproducibility purposes.
- File organization and modularization may be improved in future versions.

## Related Paper

This repository is associated with the following paper:

**Solving Laplace Equations with a Combined BEM and PINN Approach**

Please add the full citation here.

## Citation

If you use this code in your research, please cite the related paper.

```bibtex
@article{yourpaper,
  title={Solving Laplace Equations with a Combined BEM and PINN Approach},
  author={Author Name(s)},
  journal={Journal Name},
  year={2024},
  volume={},
  number={},
  pages={}
}
