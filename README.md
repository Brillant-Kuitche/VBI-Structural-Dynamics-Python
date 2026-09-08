# VBI-Structural-Dynamics-Python
Repository containing Python implementations for structural analysis, including longitudinal and transverse influence lines (Guyon-Massonet method), Finite Element Method solvers, and a 1D Vehicle-Bridge Interaction (VBI) model for dynamic behavior analysis under standard traffic loads.

Note: Inline comments within the scripts are written in French to detail the calculation steps, while the mathematical logic follows international structural engineering standards.

## Repository Structure

* `VBI 1D Beam.py` (Root directory)
  Main script modeling the dynamic vehicle-bridge interaction. It evaluates the dynamic response of a bridge deck modeled as a 1D beam under standard moving traffic loads (Fascicule 61).

* `/01_Influence_Lines`
  Scripts for evaluating bridge load distribution and generating influence lines:
  - `Influence line for a 03 span beam.py`: (Computes influence lines for continuous multi-span beams).
  - `LI Longitudinale PSI_DP.py`
  - `LI Transversales PSI_DP.py` (Implementation of the Guyon-Massonnet-Bareš method for transverse load distribution).

* `/02_Finite_Element_Method`
  Numerical solvers for structural mechanics:
  - `FEM for column with variable cross-section.py` (Assembles elementary stiffness matrices and computes the global structural response for column with varying inertia).

## Visual Portfolio

This repository is dedicated strictly to computational scripts and numerical methods. For a visual overview of my practical civil engineering projects-including automated Eurocode Excel sheets, structural modeling (Robot Structural Analysis, SAP2000), and BIM (ArchiCAD)-please refer to my complete portfolio document below.

[Download Visual Portfolio (PDF)](https://github.com/Brillant-Kuitche/Structural-Engineering-Portfolio/raw/main/Brillant_Kuitche_Portfolio.pdf)

---
*Patrice Brillant Kuitche Mbe*
