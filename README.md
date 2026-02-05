# Bow-Geometry
This repository contains the complete, reproducible numerical implementation used to generate all figures and results in the manuscript:"Geometric Signatures of Ancient Bows: A Unified Mechanical
Framework".
The code implements a physics-based model of three bow architectures — self-bow, composite recurve, and asymmetric Japanese Yumi — including bow geometry, force–draw relations, energy storage, efficiency, arrow launch velocity, and sensitivity analyses.

All figures in the manuscript can be regenerated from this repository without modification.
bow-mechanics/
├── src/
│   ├── main.py          # Entry point: reproduces all figures
│   ├── parameters.py    # Physical and geometric parameters (Table 1)
│   ├── geometry.py      # Bow–string kinematics
│   ├── forces.py        # Force–draw laws
│   ├── simulations.py  # Solvers, Yumi coupling, Monte Carlo analysis
│   └── utils.py         # Energy, efficiency, velocity utilities
│
├── figures/
│   ├── Figure-composite.png
│   └── sensitivity.png
│
├── requirements.txt
├── CITATION.cff
├── LICENSE
└── README.md

If you use this code, please cite the associated manuscript and the archived software release.
