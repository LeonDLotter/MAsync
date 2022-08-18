# Revealing the Neurobiology Underlying Interpersonal Neural Synchronization with Multimodal Data Fusion

This repository accompanies our publication: 

Leon D. Lotter, Simon H. Kohl, Christian Gerloff, Laura Bell, Alexandra Niephaus, Jana A. Kruppa, Juergen Dukart, Martin Schulte-Rüther, Vanessa Reindl, & Kerstin Konrad (2022). *[Revealing the Neurobiology Underlying Interpersonal Neural Synchronization with Multimodal Data Fusion](https://doi.org/10.1101/2022.07.26.501562)*. bioRxiv.

[![DOI](https://img.shields.io/badge/bioRxiv-10.1101/2022.07.26.501562-BD2736)](https://doi.org/10.1101/2022.07.26.501562)  
[![DOI](https://img.shields.io/badge/Zenodo-10.5281/zenodo.7002119-0F80C1)](https://zenodo.org/badge/latestdoi/407815811)  
[![Twitter](https://img.shields.io/badge/Twitter-Thread-1A8CD8)](https://twitter.com/LeonDLotter/status/1559124740450795521)

---

All data and code in this repository are licensed under [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey)](http://creativecommons.org/licenses/by-nc-sa/4.0/)  

Please cite both the preprint (bioRxiv) and the dataset (Zenodo) when using data from this repository.

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
Contains all functions called in the analysis notebook. Also includes a development version of the [JuSpyce](https://github.com/LeonDLotter/JuSpyce) toolbox.

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
  
## Reproduce the analyses

To ease the reproduction of our analyses, we added dependency management via poetry to this repository. To create a virtual python environment with all the package versions we used, perform the steps below. However, you can also install all required packages manually in a python 3.8.8 environment.  

###  Poetry on osx/ linux
Install via curl:  
```
curl -sSL https://install.python-poetry.org | python -
```
###  Poetry on windows
Install via powershell:  
```
(Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | python -
```
### Test if poetry was correctly installed:  
```
poetry --version
```
If the command is not found, add poetry to the `PATH` environment variable via:  
```
export PATH=$PATH:$HOME/.poetry/bin
```
### Run virtual environment
Go to the folder you downloaded the repository to, start a new terminal, and get synchronized :sunglasses:
```
poetry shell  
poetry install
```

## Further resources
[Here](https://leondlotter.github.io/MAsync/citenet), we provide an interactive version of the citation network included in the manuscript (Figure 2B).  
Three toolboxes were developed partly in the context of this paper:  
- [SetYouFree](https://github.com/ChristianGerloff/set-you-free): A tool for automated literature retrieval from API-accessible databases  
- [JuSpyce](https://github.com/LeonDLotter/JuSpyce): A - wait for it - s**py**ced up Python adaptation of the [JuSpace](https://github.com/juryxy/JuSpace) toolbox, incorporating numerous strategies to test for alignment between multimodal neuroimaging data  
- [ABAnnotate](https://github.com/LeonDLotter/ABAnnotate): A MATLAB toolbox for ensemble-based multimodal gene-category enrichment analysis of human neuroimaging data

## Contact
If you have any problems, questions, comments or suggestions, feel free to open an issue or [contact me](mailto:leondlotter@gmail.com)! 
