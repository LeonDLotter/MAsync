---
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
---

# Meta-Analysis and Quantitative Review of Interpersonal Neural Synchrony fMRI and fNIRS Studies

This is the repository accompanying our paper *Title*.   

*Citation*  

### Content

#### `main`
- ***[MAsync_main.ipynb](MAsync_main.ipynb)***: Jupyter notebook to reproduce all analyses that were performed in python. Contains overview descriptions of methods and results.
- **[MAsync_supp.ipynb](MAsync_supp.ipynb)**: Supplementary analyses. @Chris: erstmal nicht beachten :)

#### `src` 
Contains all functions called in the above notebooks.

- **[ma_fun.py](src/ma_fun.py)**: Functions to perform activation likelihood estimations, validation analyses, and functional decoding. 
- **[fnirs_fun.py](src/fnirs_fun.py)**: Quantitative analysis of fNIRS data.
- **[help_fun.py](src/help_fun.py)**: Additional functions to import data, manipulate images, etc.

#### `data`
Contains all results, see separate [README](data/README.md) file.

#### `fig`
Contains code to reproduce figures and the associated files.

### Getting started

## Installation
To repoduce the analysis, you must first install poetry.
Subsequently, poetry will install all other dependencies in a virtual environment.

###  Poetry on osx / linux
Install poetry via curl
```
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
```
###  Poetry windows
1.  Install poetry via powershell
```
(Invoke-WebRequest -Uri https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py -UseBasicParsing).Content | python -
```
2.  Test if poetry was installed
```
poetry --version
```
If the command was not found add poetry to the `PATH` environment variable via
```
export PATH=$PATH:$HOME/.poetry/bin
````

### start virtual environment
Go to the folder of the repo, start a new shell and become in sync :sunglasses:
```
poetry shell
```

## Contact
If you have questions or recommendations, feel free to [contact me](mailto:leondlotter@gmail.com)! 