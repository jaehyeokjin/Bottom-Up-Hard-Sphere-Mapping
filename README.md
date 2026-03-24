

# Bottom-Up Hard-Sphere Mapping of Molecular Liquids

<p align="center">
  <img src="https://jinjaehyeok.wordpress.com/wp-content/uploads/2026/03/toc-hs.png" alt="Hard-Sphere Mapping TOC" width="420">
</p>

This repository contains essential data and script files for the manuscript **“Microscopic Hard-Sphere Mapping of Molecular Liquids”** currently under review.

The work introduces a **bottom-up mapping framework** that connects realistic molecular liquids to an effective hard-sphere description using microscopic structural correlations and thermodynamic constraints. This mapping allows the dynamics of complex molecular liquids to be rationalized through a **single-parameter scaling description based on the effective packing fraction**.

The repository contains the **data, scripts, and figure files** used in the manuscript, including:

- Structural correlation functions
- Effective interaction potentials
- Hard-sphere diameter calculations
- Weeks–Chandler–Andersen (WCA) based mapping procedures
- Scripts used to generate figures and scaling analysis

---

# 1. Key Features

- **Bottom-Up Hard-Sphere Mapping**  
  Provides a bottom-up route to map molecular liquids to effective hard-sphere systems.

- **Microscopic Parameterization**  
  Uses pair correlations and effective interactions to determine hard-sphere diameters.

- **Single-Parameter Scaling**  
  Demonstrates that molecular liquids collapse onto universal scaling relations when expressed in terms of the effective packing fraction.

- **Multiple Molecular Liquids**  
  The mapping procedure is demonstrated for several representative molecular liquids:
  
  - Ortho-terphenyl (OTP)
  - Chloroform
  - triphenylchloromethane (TPC)

---

# 2. Repository Structure

```text
.
├── Data
│   ├── Chloroform
│   │   ├── Interaction
│   │   └── RDF
│   │
│   ├── OTP
│   │   ├── Interaction
│   │   └── RDF
│   │
│   └── TPC
│       ├── Interaction
│       └── RDF
│
├── Script
│   ├── Boltzmann
│   │   └── calculate_diameter.py
│   │
│   ├── Scaling-Figure1
│   │   ├── Figure1.ipynb
│   │   └── ortho_terphenyl.csv
│   │
│   └── WCA-Parameter
│       ├── Chloroform
│       ├── OTP
│       ├── TPC
│       └── Model-WCA
│
├── Figure
│   ├── Figure1.pdf
│   ├── Figure2.pdf
│   ├── Figure3.pdf
│   ├── Figure4new.pdf
│   ├── Figure5.pdf
│   ├── Figure6.png
│   ├── Figure7new.pdf
│   ├── Figure8new.pdf
│   └── Figure9new.pdf
│
└── Filelist.tree
```

---

# 3. Data

The `Data` directory contains the structural and interaction data used in the hard-sphere mapping analysis.

## Interaction Potentials

For each molecular liquid, the `Interaction` folders contain the **effective coarse-grained pair interactions** obtained at different temperatures.

Examples include

```
Data/OTP/Interaction/360.table
Data/TPC/Interaction/340.table
Data/Chloroform/Interaction/300.png
```

These potentials are used to determine **effective hard-sphere diameters** via the mapping procedures described in the manuscript.

---

## Radial Distribution Functions

The `RDF` directories contain structural correlation functions used to validate the mapping.

Typical files include:
- `cg` = coarse-grained representation  
- `fg` or `ref` = atomistic reference system

These comparisons confirm that the coarse-grained representation preserves key structural correlations.

---

# 4. Scripts

The `Script` directory contains Python scripts and Jupyter notebooks used to compute hard-sphere diameters and generate the analysis presented in the manuscript.

---

## Boltzmann Diameter

```
Script/Boltzmann/calculate_diameter.py
```

This script computes the **Boltzmann hard-sphere diameter** from the effective interaction potential.

---

## WCA Parameterization

```
Script/WCA-Parameter/
```

This directory contains scripts used to compute **Weeks–Chandler–Andersen (WCA) effective diameters** for different molecular liquids.

Subdirectories include

```
Model-WCA/ # Model system
OTP/ # Main system 1
Chloroform/ # Main system 2
TPC/ # Main system 3
```

Each system contains

- tabulated pair potentials (`*.table`)
- scripts for determining the WCA diameter
- Jupyter notebooks for visualization and analysis.

---

## Single-Parameter Scaling (Figure 1)

```
Script/Scaling-Figure1/
```

This directory contains the notebook used to generate the **single-parameter scaling collapse** shown in Figure 1 of the manuscript. The analysis demonstrates that the dynamics of molecular liquids collapse when expressed in terms of the **effective hard-sphere packing fraction**.

---

# 5. Figures

The `Figure` directory contains the final figures used in the manuscript.

---

# 6. Citation

If you use any data or scripts from this repository, please cite the manuscript: **Microscopic Hard-Sphere Mapping of Molecular Liquids** (arXiv will be updated)

```bibtex
@article{jin2026hardsphere,
  title={Microscopic Hard-Sphere Mapping of Molecular Liquids},
  author={Jin, Jaehyeok and Pedersen, Ulf R. and Dyre, Jeppe C. and Reichman, David R.},
  journal={arXiv preprint arXiv:XXXX.XXXXX},
  year={2026}
}
```

