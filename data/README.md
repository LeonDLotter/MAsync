# The Neurobiological Basis of Interpersonal Neural Synchronization: Evidence from a Multimodal Data Fusion Approach

Here, you find all data underlying and generated in [MAsync_analyses.ipynb](../MAsync_analyses.ipynb)

DOI: ...

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]  

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

---

## Content

#### `/ale`
Contains results of the ALE analyses, with the following naming convention.  
Analyses:
- **ale_{}**: fMRI data, voxel-level threshold p < 0.001 (default)
- **ale01_{}**: fMRI data, voxel-level threshold p < 0.01
- **aleFNIRS_{}**: fMRI & all fNIRS data
- **aleFNIRSSEL_{}**: fMRI & selected fNIRS data
- **aleNoPseudo_{}**: fMRI data, excluding pseudo-hyperscanning experiments  

Output types:
- **{}_p.nii.gz, {}_stat.nii.gz, {}_z.nii.gz**: unthresholded p-, ale-, p-to-z-values
- **{}_{logp | z}_desc-{mass | size}_level-cluster_corr-FWE_method-montecarlo.nii.gz**: cluster-level-corrected volumes of z or negative log(p) values
- **{}_{logp | z}_level-voxel_corr-FWE_method-montecarlo.nii.gz**: voxel-level-corrected volumes of z or negative log(p) values
- **{}_thresh_4d.nii.gz**: volumes thresholded at cluster-level p < 0.05, 4D-volume with one binarized volume for each cluster
- **{}_thresh_bin.nii.gz**: volumes thresholded at cluster-level p < 0.05, all clusters binarized
- **{}_thresh_idx.nii.gz**: volumes thresholded at cluster-level p < 0.05, clusters labelled from 1 to n_clusters
- **{}\_thresh_{stat | z}.nii.gz**: volumes thresholded at cluster-level p < 0.05, ale or z-values, not binarized
- **ale_thresh_{}.csv**: tables with anatomical cluster and peak information

#### `/atlases`
Used templates, parcellations, and meta-analytic maps.

#### `/context`
Results from contextualization analyses: Distributions of the INS map within meta-analytic maps, conjunction volumes with prior meta-analytic maps, neurosynth topic decoding results, results from spatial correlation with PET and mRNA expression maps.
- **/gcea**: results from gene-category enrichment analyses using [ABAnnotate](https://github.com/LeonDLotter/ABAnnotate)
- **/parcel-maps**: parcellated volumes generated from parcellated PET data for analysis with [JuSpace](https://github.com/juryxy/JuSpace)
- **/topic_maps**: meta-analytic maps generated from Neurosynth LDA topic-decoding selected from the 200-topic LDA dataset

#### `/datasets`
Datasets used for all analyses. **This includes the coordinate data for fNIRS and fMRI studies collected in this work!**
- **/neurosynth**: the downloaded neurosynth database version 7
- **brainmap_NormMapActOnly_{}.pkl.gz**: "normal mapping" & "activation only" experiments from the BrainMap database as a NiMARE dataset
- **fmri_coordinates**: INS fMRI coordinates and associated study information
- **fnirs_coordinates**: INS FNIRS coordinates and associated study information
- **pet_parcellated**: parcellated PET data

#### `/fnirs`
Contains fNIRS INS results.
- **fnirs_atlas_result_{}.csv**: parcel-wise summarized fNIRS data, 'sel': only post-hoc selected studies
- **fnirs_atas_{}_{}.nii.gz**: volumes with parcel-wise results (see notebook)
- **fnirs_coordinates_nearest_{}.csv**: original and nearest within-atlas fNIRS coordinates
- **fnirs_ale_contributions.csv**: contributions of fNIRS studies to clusters obtained in joines fMRI+fNIRS ALE
- **fnirs_ale_randomization.csv**: result from repeated ALE analysis when fNIRS coordinates were slightly randomized

#### `/fsn`
Contains ALE results (see above) of the Fail-Safe-N analysis, only last iteration (original coordinates + the largest number of noise studies that can be added to the original ALE before the rTPJ cluster looses significance).

#### `/loeo`
ALE results (see above) of all fMRI Leave-One-Experiment-Out ALEs. Additionally, the conjunction volume of all LOEO ALEs and csvs with experiment-wise contributions to the fMRI ALE clusters

#### `/macm`
Contains inputs and results for/of MACM analyses (for naming see above).

#### `/rsfc`
Contains RSFC results.
- **fwd.csv**: motion estimates (frame-wise displacement) for exclusion of RSFC subjects 
- **macm_rTPJ_idx_forConn.nii/.csv**: uncompressed MACM volume with associated cluster labels
- **rsfc_cormat{}.csv**: averaged, thresholded and unthresholded semipartial z-transformed correlation matrix between MACM-ROIs
- **rsfc_macm_rTPJ_semipartial.mat**: RSFC result from CONN toolbox
