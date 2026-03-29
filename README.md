# Bow Geometry and Energy Transfer – Reproducibility Package

This repository accompanies the manuscript:

**"Geometric Signatures of Ancient Bows: A Unified Mechanical Framework"**
Authors: Kishore Dutta


---

## Overview

This repository provides the complete numerical implementation supporting all results and figures reported in the accompanying manuscript. Its purpose is to ensure full transparency, reproducibility, and long-term accessibility of the computational analysis.

The code implements a physics-based model of bow mechanics that captures:

* Limb geometry and curvature
* Rotational stiffness distribution
* String kinematics
* Energy storage and transfer mechanisms

Three bow architectures are analyzed:

* Self-bow
* Composite recurve bow
* Asymmetric Japanese Yumi

All force, draw relationships, stored energy, delivered energy, efficiency metrics, and sensitivity analyses reported in the manuscript are generated directly from this implementation.

The repository is fully self-contained and deterministic, enabling exact reproduction of results across systems.

---

## Repository Structure

* `parameters.py`
  Defines all physical and geometric parameters used in the model

* `geometry.py`
  Implements geometric relationships, including limb angles and draw length calculations

* `forces.py`
  Computes force--draw curves, stored energy, and energy transfer

* `simulations.py`
  Runs full simulations and parameter studies for different bow configurations

* `utils.py`
  Contains helper and utility functions

* `main.py`
  Entry point that executes all simulations and generates manuscript figures

* `requirements.txt`
  Lists required Python packages

* `README.md`
  Documentation

---

## Reproducing the Results

All results in the manuscript can be reproduced by running:

python main.py

This script:

* Runs simulations for all three bow types
* Computes force--draw curves
* Calculates stored and transferred energy
* Performs efficiency and sensitivity analyses
* Generates all figures used in the manuscript

---

## Mapping Between Code and Manuscript

* **Geometric formulation and bow shape modeling**
  → `geometry.py`

* **Force--draw curves (Figures showing draw force vs draw length)**
  → `forces.py`, `simulations.py`

* **Energy storage and transfer calculations**
  → `forces.py`

* **Efficiency analysis and comparative results across bow types**
  → `simulations.py`

* **Sensitivity analysis (parameter variation studies)**
  → `simulations.py`

All figures in the manuscript are generated automatically by `main.py`.

---

## Software Requirements

Python version:
Python 3.10.12 

Required packages:

* numpy
* scipy
* matplotlib

Install dependencies using:

pip install -r requirements.txt

---

## Getting Started

1. Clone or download the repository

2. Navigate to the project directory

3. (Optional) Create a virtual environment:
   python -m venv env

4. Activate the environment:

   * Linux/Mac: source env/bin/activate
   * Windows: env\Scripts\activate

5. Install dependencies:
   pip install -r requirements.txt

6. Run the simulation:
   python main.py

---

## Reproducibility Notes

* The code is fully self-contained and does not require external datasets
* Deterministic execution is ensured (`np.random.seed(0)` is set where applicable)
* Numerical solvers and procedures are explicitly defined
* Results should reproduce exactly given consistent software versions

---

## Data and Code Availability

This repository is archived for long-term access at:

DOI: XXXXX (to be added after Zenodo archival)

Development version (GitHub):
http://github.com/kdkishore77-alt/Bow-Geometry.git

---

## Intended Use

This code is intended for:

* Reviewers verifying reproducibility
* Researchers reproducing or extending the analysis
* Studies of historical, archaeological, or biomechanical bow systems

The implementation is designed to be accessible to users with basic familiarity with Python.
