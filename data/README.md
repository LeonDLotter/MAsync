These are the majority of results published in our manuscript *Title*.   

### Content:
- **[ale_coordinates.csv](ale_coordinates.csv)**: table with fMRI study information and MNI coordinates.
- **[ale_rTPJ_bin.nii.gz](ale_rTPJ_bin.nii.gz)**: binarized cluster volume resulting from main ALE analysis.
- **[fnirs_atlas_result.csv](fnirs_atlas_result.csv)**: table with ROI-wise fNIRS results.
- **[fnirs_coordinates.csv](fnirs_coordinates.csv)**: table with fNIRS study information and reconstructed MNI coordinates.
- **[macm_rTPJ_bin.nii.gz](macm_rTPJ_bin.nii.gz)**: binarized network volume resulting from MACM anylysis.
- **[macm_rTPJ_idx.nii.gz](macm_rTPJ_idx.nii.gz)**: MACM volume with indexed clusters.


- **ale/..**: results of the main ALE analyses ('ale_[...]': voxel-level p < 0.001, 'ale01_[...]': p < 0.01).
  - **ale_logp_[...].nii.gz, ale_z_[...].nii.gz**: cluster/voxel-corrected NiMARE outputs.
  - **ale_p.nii.gz, ale_stat.nii.gz, ale_z.nii.gz**: unthresholded NiMARE output.
  - **ale_thresh_[...].nii.gz**: volumes thresholded at cluster-level p < 0.05.
  - **ale_thresh_[...].csv**: table with cluster and peak information.


- **atlases/..**: atlases/ROIs/networks used in our analyses.


- **distr/..**: spatial overlap between INS results and predefined networks/rois quantified as relative and absolute distributions.


- **fd/..**: functional decoding results.


- **fnirs/..**: fNIRS INS results.
  - **fnirs_atas_[...].nii.gz**: volumes with parcel-wise results (see notebook).
  - **fnirs_coordinates_nearest.csv**: original and nearest within-atlas fNIRS coordinates.


- **fsn/..**: ALE results (see above) from the Fail-Safe-N analysis with the largest number of noise studies that can be added to the original ALE before the rTPJ cluster looses significance.


- **loeo/..**: ALE results (see above) from all Leave-One-Experiment-Out ALEs.


- **macm/..**: MACM input and results.
  - **ale_rTPJ_bin_forSleuth.nii.gz**: ALE rTPJ cluster with 1x1x1mm^3 resolution 
  - **cor_macm_ale.csv**: parcel-wise ALE value estimates for ALE-MACM volume correlations
  - **macm_rTPJ_[...].nii.gz**: ALE results (see above) from MACM analyses
  - **macm_rTPJ_coordinates.txt**: Sleuth text file with MNI coordinates for MACM ALE
  - **macm_rTPJ_project.work**: Sleuth work project  

- **rsfc/..**: RSFC results
  - **ale_rTPJ_bin_forCONN.nii**: ALE rTPJ cluster as uncompressed volume for CONN toolbox
  - **fwd.csv**: motion estimates (frame-wise displacement) for exclusion of RSFC subjects 
  - **macm_rTPJ_idx_forConn.nii/.csv**: uncompressed MACM volume with associated cluster labels
  - **rsfc_cormat.csv**: averaged and thresholded semipartial z-transformed correlation matrix 
  - **rsfc_macm_rTPJ_semipartial.csv**: RSFC result from CONN toolbox
