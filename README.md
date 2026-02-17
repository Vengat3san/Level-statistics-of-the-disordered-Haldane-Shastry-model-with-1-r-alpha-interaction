# Exact Diagonalization (ED) of the generalized Haldane-Shastry model with $1/r^{\alpha}$ interaction.

   

This repository contains the source codes for Exact Diagonalization (ED) of the generalized generalized Haldane-Shastry model with $1/r^{\alpha}$ interaction, the datasets for the mean gap ratio, and the Jupyter notebooks used to reproduce the plots presented in our paper: [Link to Paper Here]().

# Table of Contents

1. [The Model](#the-model)
2. [Repository Structure](#repository-structure)
3. [How to access the data](#how-to-access-the-data)
4. [Running on Google Colab (Non-Julia Users)](#running-on-google-colab-non-julia-users)

# The Model

**Clean limit:** We consider the following variant of the Haldane-Shastry model described by the Hamiltonian

$$
H_0 = \frac{1}{\mathcal{N}}\left( \frac{2\pi }{N}\right)^{\alpha} \sum_{1\leq j{<}k \leq N} \frac{\vec{S}_j \cdot \vec{S}_k}{\left|\eta_j-\eta_k\right|^{\alpha}},
$$



where $\mathcal{N}$ is a scale factor, $\eta_j{=}e^{i\frac{2\pi}{N} j}$ is the position of the $j^{\rm th}$ spin in the complex plane, $\left|\eta_j{-}\eta_k\right|{=}r$ is the chord distance between the sites $j$ and $k$. 

As $\alpha$ increases from $0$, this model interpolates between the all-to-all Heisenberg model ($\alpha{=}0$) and the XXX model ($\alpha{\to}\infty$), with $\alpha{=}2$ being the Haldane-Shastry model. 

**Position disorder:** Introduced by randomly displacing the spins from their original sites in such a way that the neighboring spins do not cross each other.

$$
H_\delta= \frac{1}{\mathcal{N}}\left( \frac{2\pi }{N}\right)^{\alpha} \sum_{1\leq j{<}k \leq N} \frac{\vec{S}_j \cdot \vec{S}_k}{\left|\eta_j^{\delta}-\eta_k^{\delta}\right|^\alpha},
$$

where $\eta_j^{\delta}{=}e^{i\frac{2\pi}{N}(j{+}\zeta^{\delta}_j)}$ is the new displaced position of the $j^{\rm{th}}$ spin with $\zeta^{\delta}_j$ being a uniformly distributed random number in the interval $[{-}\delta/2,\delta/2]$, and $\delta{\in}[0, 1]$ controls the strength of the position disorder. 

**Field disorder:** Introduced by the random longitudinal magnetic field term $H_h{=}(1/\mathcal{N})\sum h_j\, S_j^z$, where $h_j$ is a Gaussian distributed random number with mean $0$ and variance $h^2$. This term breaks the $SU(2)$ symmetry of the Hamiltonian down to $U(1)$. The field strength is measured relative to the interaction energy between nearest-neighbor (NN) spins of the clean Hamiltonian i.e., $J^{\rm NN}_\alpha$.

$$
H=H_\delta+H_h =\frac{1}{\mathcal{N}}\left[\left( \frac{2\pi }{N}\right)^{\alpha} \sum_{1\leq j{<}k \leq N} \frac{\vec{S}j \cdot \vec{S}k}{\left|\eta_j^{\delta}-\eta_k^{\delta}\right|^\alpha} + \sum_{1\leq j \leq N} h j \, S_j^z\right].
$$

# Repository Structure
```
.  
├── data/                   # Mean gap ratio datasets (JLD2 format)
│   ├── (datasets)
│   └── (notebooks for loading data and generating plots)
│
├── scripts/                # Example scripts for ED calculations
│   └── (example scripts)
│
├── src/                    # Julia source codes for ED calculations
│   └── (source files)
│
├── Manifest.toml           # Exact dependency versions
├── Project.toml            # Julia project dependencies
└── README.md               # Project documentation
```
# How to access the data

To run these scripts locally, you need a working Julia installation.
### Clone the repository:

```bash
git clone https://github.com/Vengat3san/Level-statistics-of-the-disordered-Haldane-Shastry-model-with-1-r-alpha-interaction.git
```

### Install dependencies:

The repository includes a Julia environment (Project.toml) to ensure consistent package versions across systems. From the repository root directory, run:

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

After installation, open the Jupyter notebooks located in the data/ folder. These notebooks contains the code for loading the datasets and reproducing the plots used in the manuscript.
# Running on Google Colab (Non-Julia Users)

If you are not familiar with the Julia ecosystem, you can still reproduce our results and plots using Google Colab:

1. Download the data files (.jld2) and the notebooks (.ipynb) from the data/ folder.
2. Go to Google Colab.
3. Upload the notebook and the data files to your session.
4. Set the runtime to Julia.
5. Excute the notebook cells to generate the plots from the manuscript.
