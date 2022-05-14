# The Neurobiological Basis of Interpersonal Neural Synchronization: Evidence from a Multimodal Data Fusion Approach

This repository accompanies our preprint: 

Leon D. Lotter, Simon Huldreich Kohl, Christian Gerloff, Laura Bell, Alexandra Niephaus, Jana A. Kruppa, Juergen Dukart, Martin Schulte-Rüther, Vanessa Reindl, & Kerstin Konrad (2022). *[The Neurobiological Basis of Interpersonal Neural Synchronization: Evidence from a Multimodal Data Fusion Approach]()*. bioRxiv.

[![DOI: XXX](https://img.shields.io/badge/DOI-XXX-blue)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

---

All data and code in this repository are licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey)](http://creativecommons.org/licenses/by-nc-sa/4.0/)  

---
## Main analyses and results:

- ***[MAsync_analyses.ipynb](MAsync_analyses.ipynb)***: Jupyter notebook to reproduce all analyses except for resting-state data processing. Contains general descriptions of methods and results.
- ***[data/datasets/fmri_coordinates.csv](data/datasets/fmri_coordinates.csv)***: Coordinates and associated study information from fMRI interpersonal neural synchrony studies
- ***[data/datasets/fnirs_coordinates.csv](data/datasets/fnirs_coordinates.csv)***: Coordinates and associated study information from fNIRS interpersonal neural synchrony studies

## All data, code, figures, & literature search result:

#### `/data`
Contains all results, see separate [README](data/README.md) file.

#### `/fig`
Contains code to reproduce figures and the associated files.

#### `/lit` 
Contains results from different stages of the literature search as *.ris files.

#### `/src` 
Contains all functions called in the analysis notebook.

- **[ale.py](src/ale.py)**: Wrapper for Activation Likelihood Estimation workflow in NiMARE
- **[fnirs.py](src/fnirs.py)**: Quantitative analysis of fNIRS data
- **[fsn.py](src/fsn.py)**: Fail-Safe-N validation analysis
- **[loeo.py](src/loeo.py)**: Leave-one-experiment-out validation analysis
- **[neurosynth.py](src/neurosynth.py)**: Get Neurosynth data and generate meta-analytic maps
- **[overlap.py](src/overlap.py)**: Quantify overlap between cluster volumes using relative and absolute distributions
- **[parcellate_pet.py](src/parcellate_pet.py)**: Obtain parcellated PET data from volumes (volumes not in this repo)
- **[rsfc_resample_hcp.py](src/rsfc_resample_hcp.py)**: Resample HCP-data to 3mm voxel-resolution
- **[utils_image.py](src/utils_image.py)**: Image manipulation functions
- **[utils_io.py](src/utils_io.py)**: NiMARE dataset input/output
- **[utils_plot.py](src/utils_plot.py)**: Plotting functions used within analysis notebook
- **[CONNbatch_processHCP.m](src/CONNbatch_processHCP.m)**: Matlab: Batch script to generate resting-state functional connectivity data from FIX-preprocessed HCP-data
- **[juspace_correlations.m](src/juspace_correlations.m)**: Matlab: Wrapper script calling functions from the JuSpace toolbox for spatial correlation analyses  
  
## Contact
If you have questions, comments or suggestions, feel free to [contact me](mailto:leondlotter@gmail.com)! 
