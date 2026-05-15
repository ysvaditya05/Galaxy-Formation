# Thesis Code Repository

## Author
Yeddula S. V. Aditya

## Acknowledgements

This work was carried out under the guidance of Prof. Kandaswamy Subramanian.

## Description
This repository contains the code developed as part of my thesis work on High-Redshift Galaxy Formation. It includes:
- Star Formation Rate Density (SFRD)
- Cosmological Reionization
- X-ray Background (CXB)
- X-ray Heating of the Intergalactic Medium (IGM)
- Alternate Sources of X-ray Emission

The code implements theoretical models, numerical integrations, and plotting routines used to generate the results presented in the thesis.

---

## Repository Structure

/Data

/Notebooks
    ├── .ipynb_checkpoints/
    ├── Alternat_XRaySources.ipynb
    ├── Cosmological_Reionization.ipynb
    ├── Current_XRay_Background.ipynb
    ├── Theoretical_Setup.ipynb
    ├── XRay_Heating.ipynb

/src
    ├── __pycache__/
    ├── .ipynb_checkpoints/
    └── utils.py

---

## Requirements

- Python
- NumPy
- SciPy
- Matplotlib

---

## How to Run

1. Download the parent folder as a whole
2. Navigate to the /Notebooks directory
3. Run the desired code file; the contents of each code file are specified below

---

## Key Features

- Theoretical_Setup.ipynb contains all code related to the theoretical foundations — PS formalism, abundance calculation, and SFRD calculation
- Cosmological_Reionization.ipynb contains the code related to cosmological reionization and Thomson optical depth
- Current_XRay_Background.ipynb contains the code related to the calculation of the contribution of star-forming galaxies to the CXB at present day (z = 0)
- XRay_Heating.ipynb contains the code related to X-ray heating of the IGM
- Alternate_XRaySources.ipynb contains the code related to alternate sources of X-ray emission, including accretion shocks and emissivities of HMXBs and LMXBs

---

## Notes

- Some scripts may depend on intermediate outputs from other scripts. This is managed by /src/utils.py, where shared functions are defined
- Some parts of the code take a long time to run. To avoid repeated execution, intermediate outputs are stored in the /Data directory
- Files in the /Data directory do not need to be edited manually; updates can be made through the notebooks and saved accordingly
- The code is in its early stages and contains limited comments
- This repository may be updated with additional code and data

---

This code was developed as part of my thesis on High-Redshift Galaxy Formation.

Yeddula, S. V. Aditya (2026)  
Code Repository for High-Redshift Galaxy Formation (UG Thesis)