#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:24:59 2021

@author: Leon D. Lotter
"""

import nilearn as nl
from nilearn.image import resample_to_img, math_img, load_img, new_img_like
from nilearn.glm import threshold_stats_img
from nimare.dataset import Dataset
from nimare.meta.cbma.ale import ALE
from nimare.stats import null_to_p
from nimare.transforms import p_to_z
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import pandas as pd
import numpy as np

from .utils_image import get_cluster_stats

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

    res_df = pd.concat(res_out).reset_index()
    res_df = res_df[['roiName', 'nvox_roi', 'targetNr', 'targetName', 
                     'nvox_target', 'nvox_overlap', 'relDistr', 'absDistr']]
    return(res_df)


#=============================================================================


def distribution_null_test(ds, target_vol, roi_vol=None, 
                           vox_thresh=0.001, cluster_thresh=25, thresh_type="size", 
                           n_perm=1000, seed=None, n_proc=-1,
                           target_labs=None, template=None):
    """
    Assignes permutation-based p-values to the overlap between meta-analytically derived ROIs with, 
    e.g., a given target network parcellation. 
    Relative distribution = proportion of 'roi_vols' voxels within a given 
        'target_vol' cluster vs. all 'roi_vols' voxels.
    Absolute distribution  = proportion of 'roi_vols'-voxels within a given 
        'target_vol' cluster vs. all voxels within this 'target_vol' cluster.
    If the ROI volume {roi_vol} is None, a true ROI map is calculated from the input dataset {ds}. 
    The 'true' distributions across the target ROIs {target_vol} are calculated. Based on the 
    numbers of subjects and foci in {ds}, {n_perm} null datasets are randomly created. On each of
    these datasets, an ALE ROI map is calculated using the thresholding settings {vox_thresh}, 
    {cluster_thresh} and {thresh_type}. Relative/Absolute distributions are calculated for each
    of these null ROI maps and the resulting null distributions are compared to the true values
    to estimate positive-sided p-values for each ROI in {target_vol}.
    
    Input:
        ds = NiMARE dataset
        target_vol = target volume with clusters/parcels indexed from 1 to n
        roi_vol = input volume. If None, will be generated from NiMARE dataset
        vox_thresh = voxel threshold as p-value
        cluster_thresh = cluster mass or size threshold, depending on {thresh_type}
        thresh_type = 'mass' or 'size'
        n_perm = number of iterations/ null datasets/ maps
        seed = seed for reproducability
        n_proc = number of parallel processes, -1 = number of CPU cores
        target_labs = list of target cluster/parcel names for output table
        template = template (i.e. grey matter) from which to sample random foci coordinates.
            if None, will use nilearn MNI-152 GM TPM thresholded at > 0.2

    Output:
        - pandas dataframe with p-values corresponding to relative and absolute distributions  
            within each cluster of the target volume.  
        - dict with {n_perm} pandas dataframes as resulting from the rel_abs_distributions method
    """
    
    def thresh_img(img_stat, img_z):    
        # voxel threshold 
        img_vthresh = math_img(
            f"stat * (np.abs(z) > {p_to_z(vox_thresh, tail='one')})", 
            stat=img_stat,
            z=img_z)
        # cluster threshold
        img_labels, clust_labels, clust_sizes, clust_masses = get_cluster_stats(img_vthresh)
        if thresh_type=="size":
            clust_labels = clust_labels[clust_sizes > cluster_thresh]
        elif thresh_type=="mass":
            clust_labels = clust_labels[clust_masses > cluster_thresh]
        img_cthresh = np.zeros(img_vthresh.shape)
        for c in clust_labels:
            img_cthresh[img_labels==c] = 1
        # return new image
        return new_img_like(img_vthresh, img_cthresh)
    
    # run real ALE
    lgr.info(f"Estimating 'true' results. Thresholding: voxel p < {vox_thresh}, cluster {thresh_type} > {cluster_thresh}.")
    fit_ale = ALE()
    if roi_vol is None:
        ale_true = fit_ale.fit(ds)
        # apply thresholds
        roi_vol = thresh_img(img_stat=ale_true.get_map("stat"), img_z=ale_true.get_map("z"))
    # get real distributions
    distr_true = rel_abs_distributions(roi_vols=[roi_vol],
                                       target_vol=load_img(target_vol),
                                       target_labs=target_labs)

    # number of experiments, foci, and  subjects for null experiments
    np.random.seed = seed
    n_exp_true = len(ds.ids)
    n_coord_true = [len(ds.coordinates.x[ds.coordinates.id==id]) for id in ds.ids]
    n_sub_true = [n[0] for n in ds.get_metadata(field='sample_sizes')]
    
    # Get coordinates for null experiments
    # If no template given, fetch MNI152 GM mask and resample to target_vol space
    # produces warning: "pixdim[0] (qfac) should be 1 or -1; setting qfac to 1",
    # the result is okay, however.
    if template is None:
        template_1mm = nl.datasets.fetch_icbm152_brain_gm_mask(threshold=0.2)
        template = nl.image.resample_to_img(template_1mm, target_vol, 
                                            interpolation='nearest')
    # Get possible MNI coordinates for null studies = all voxels in template
    dat = template.get_fdata() # get matrix
    x, y, z = np.where(dat == 1.0) # get GM coordinates
    affine = template.affine # get affine matrix for MNI transform
    # get mni coordinates
    mni_all = nl.image.coord_transform(x=x, y=y, z=z, affine=affine) 
    mni_all = np.array(mni_all).T
    
    # Function to generate null datasets
    def null_datasets(i):
        # number of 
        n_coord_null = np.random.choice(n_coord_true, n_exp_true)
        n_sub_null = np.random.choice(n_sub_true, n_exp_true)
        # generate null experiments
        exp_null_dict = {} # empty dict to store null study-wise results
        for ii in range(n_exp_true):  
            study_id = f'null_{ii+1:03}' # study id
            coords_list = mni_all[np.random.choice(mni_all.shape[0], n_coord_null[ii])] # pick random foci
            coords = np.vstack(coords_list) # make coordinates array
            # write into dict as required by NiMARE
            exp_null_dict[study_id] = {
                'contrasts': { '1': {
                    'coords': {
                        'space': 'MNI',
                        'x': list(coords[:, 0]),
                        'y': list(coords[:, 1]),
                        'z': list(coords[:, 2]),
                        },
                    'metadata': { 
                        'sample_sizes': n_sub_null[ii]
                        }}}}
        # convert dict to NiMARE dataset
        ds_null = Dataset(exp_null_dict) 
        ds_null.coordinates['space'] = 'MNI' # avoids error when merging datasets
        # return
        return ds_null
    
    # Generate null data
    lgr.info(f"Generating null datasets (n = {n_perm}).")
    ds_null = Parallel(n_jobs=n_proc)(delayed(null_datasets)(i) for i in tqdm(range(n_perm)))
    ds_null = dict(zip(list(range(n_perm)), ds_null))
            
    # Function to estimate ALE, calculate distributions and return result
    def null_ale(ds_null):
        # perform ALE
        ale_null = fit_ale.fit(ds_null)
        ale_null_thresh = thresh_img(img_stat=ale_null.get_map("stat"), img_z=ale_null.get_map("z"))
        if (ale_null_thresh.get_fdata() > 0).any():
            # get distributions
            distr_null = rel_abs_distributions(roi_vols=[ale_null_thresh],
                                               target_vol=load_img(target_vol))
        else:
            distr_null = None
        # return
        return distr_null
        
    # Get distributions of null volumes across target volumes
    lgr.info(f"Performing null meta-analyses and estimating distributions across target volume.")
    distr_null = Parallel(n_jobs=n_proc)(delayed(null_ale)(ds_null[i]) for i in tqdm(range(n_perm)))
    distr_null = dict(zip(list(range(n_perm)), distr_null))

    # collect distribution data
    lgr.info("Estimating p-values.")
    p_rel, p_abs = dict(), dict()
    for i_target, target in enumerate(distr_true.targetName):
        data_null = np.zeros((n_perm,2))
        for i in range(n_perm):
            if distr_null[i] is not None:
                data_null[i,:] = distr_null[i].loc[i_target, ["relDistr", "absDistr"]]
        p_rel[target] = null_to_p(distr_true["relDistr"][i_target], data_null[:,0], tail="upper")
        p_abs[target] = null_to_p(distr_true["absDistr"][i_target], data_null[:,1], tail="upper")

    # return
    return pd.DataFrame([p_rel, p_abs], index=["relDistr", "absDistr"]).T, distr_null
    
    

