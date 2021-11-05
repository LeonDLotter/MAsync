---
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
---

# Meta-Analysis and Quantitative Review of Interpersonal Neural Synchrony fMRI and fNIRS StudiesThis is the repository accompanying our paper *Title*.   *Citation*  
### Content

#### `main`
- ***[MAsync_main.ipynb](MAsync_main.ipynb)***: Jupyter notebook to reproduce all analyses that were performed in python. Contains overview descriptions of methods and results.
- **[MAsync_supp.ipynb](MAsync_supp.ipynb)**: Supplementary analyses. @Chris: erstmal nicht beachten :)
#### `src` 
Contains all functions called in the above notebooks.
- **[ma_fun.py](src/ma_fun.py)**: Functions to perform activation likelihood estimations, validation analyses, and functional decoding. - **[fnirs_fun.py](src/fnirs_fun.py)**: Quantitative analysis of fNIRS data.- **[help_fun.py](src/help_fun.py)**: Additional functions to import data, manipulate images, etc.

#### `data`Contains all results, see separate [README](data/README.md) file.
#### `fig`
Contains code to reproduce figures and the associated files.
  
### ContactIf you have questions or recommendations, feel free to [contact me](mailto:leondlotter@gmail.com)! 