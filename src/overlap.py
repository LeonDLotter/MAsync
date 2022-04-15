#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:24:59 2021

@author: Leon D. Lotter
"""

from nilearn.image import resample_to_img, math_img
import pandas as pd
import numpy as np
    
import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)


def rel_abs_distributions(roi_vols, target_vol, 
                          roi_labs=None, target_labs=None):   
    """
    Quantifies overlap of a list of ROIs with, e.g., a given target network 
    parcellation. 
    Relative distribution = proportion of 'roi_vols' voxels within a given 
        'target_vol' cluster vs. all 'roi_vols' voxels.
    Absolute distribution  = proportion of 'roi_vols'-voxels within a given 
        'target_vol' cluster vs. all voxels within this 'target_vol' cluster.
    
    Input:
        roi_vols = list of roi volumes 
        roi_labs = list of roi names for output table
        target_vol = target volume with clusters/parcels indexed from 1 to n
        target_lab = list of target cluster/parcel names for output table
    All volumes are resampled to space and voxelsize of the first roi.

    Output:
        pandas dataframe with relative and absolute distributions of each ROI 
        within each cluster of the target volume.     
    """
    
 
    # resample target volume to space of first roi volume 
    t_vol = resample_to_img(target_vol, roi_vols)
    
    t_mat = t_vol.get_fdata().round() # get target data matrix
    t_idc = np.unique(t_mat)[1:].astype(int) # get cluster indices
    
    # target labels
    if target_labs is None:
        target_labs = [f'target{i}' for i in t_idc]
    if len(t_idc) != len(target_labs):
        lgr.error(f'Number of target clusters/parcels ({len(t_idc)}) does not '
                  f'match number of given labels ({len(target_labs)})!')
            
    # roi labels
    if roi_labs is None:
        roi_labs = [f'roi{i}' for i in range(0, len(roi_vols))]
    if len(roi_vols) != len(roi_labs):
        lgr.error(f'Number of ROI volumes ({len(roi_vols)}) does not '
                  f'match number of given labels ({len(roi_labs)})!')

    # Results table with target cluster names and total voxel numbers per cluster
    res = pd.DataFrame({'targetNr': t_idc, 
                        'targetName': target_labs,
                        'nvox_target': [len(t_mat[t_mat == idx]) for idx in t_idc]})  
    
    # Iterate over ROI list
    res_out = []
    for roi, roi_lab in zip(roi_vols, roi_labs):
        r_vol = resample_to_img(roi, roi_vols[0]) # resample to first roi space
        r_bin = math_img('img > 0', img=r_vol) # binarize, just in case
        r_mat = r_bin.get_fdata() # get roi data matrix
        rt_mat = r_mat * t_mat # intersect roi with target matrix
        
        # RESULTS
        res_roi = res.copy()
        # ROI name
        res_roi['roiName'] = roi_lab
        # voxels in ROI
        res_roi['nvox_roi'] = len(r_mat[r_mat == 1])
        # overlapping voxels per target cluster
        res_roi['nvox_overlap'] = [len(rt_mat[rt_mat == nr]) for nr in res_roi.targetNr]
        # distributions
        res_roi['relDistr'] = res_roi.nvox_overlap / res_roi.nvox_roi # relative
        res_roi['absDistr'] = res_roi.nvox_overlap / res_roi.nvox_target # absolute
        
        # save to list
        res_out.append(res_roi)

    res_df = pd.concat(res_out)
    res_df = res_df[['roiName', 'nvox_roi', 'targetNr', 'targetName', 
                     'nvox_target', 'nvox_overlap', 'relDistr', 'absDistr']]
    return(res_df)