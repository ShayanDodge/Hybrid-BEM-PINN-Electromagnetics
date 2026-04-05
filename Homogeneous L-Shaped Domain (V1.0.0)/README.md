# Homogeneous L-Shaped Domain: 

This version contains the first implementation of a coupled **Boundary Element Method (BEM)** and **Physics-Informed Neural Network (PINN)** framework for solving the **Laplace equation** on a **homogeneous L-shaped domain**.

## Overview
<table>
  <tr>
    <td valign="top" width="60%">

### Domain decomposition

The computational domain is divided into two subdomains:

- **Left subdomain:** solved using a PINN  
- **Right-bottom subdomain:** solved using BEM  

The two parts are coupled through interface and boundary information to reconstruct the full solution over the L-shaped domain.

</td>
    <td valign="top" width="40%" align="center">
      <img src="https://github.com/user-attachments/assets/ed3e7a41-335a-4161-a946-79d129138d8a" alt="L-shaped domain decomposition" width="320">
    </td>
  </tr>
</table>

## Main Idea

This code performs the following steps:

1. Constructs a 2D L-shaped domain from a rectangular grid  
2. Extracts Dirichlet, Neumann, and interface boundary points  
3. Trains a PINN model on one subdomain  
4. Uses PINN-predicted interface data as input to the BEM solver  
5. Solves the complementary subdomain using BEM  
6. Reconstructs and visualizes the combined scalar field  

## Repository Structure

    homogeneous-l-shaped-domain-bem-pinn/
    ├── README.md
    ├── LICENSE
    ├── .gitignore
    ├── main/
    │   └── bem_pinn_laplace_lshape.mlx
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
    └── results/
        ├── data.mat
        └── data_bem.mat


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
