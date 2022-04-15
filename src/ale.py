#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:21:15 2021

@author: Leon D. Lotter
"""

import os
from os.path import join
import numpy as np
import nimare as ni
import nilearn as nl
import math
from atlasreader import get_statmap_info
from .utils_image import index_clusters_bin

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)


def ale(data, work_dir=None, pref='ALE_', 
        vox_thr=0.001, n_perm=1000, n_core=-1, cluster='mass', output=True, 
        save_cl_info=True, atlases=['talairach_gyrus', 'talairach_ba', 'aal'],
        print_glassbrain=True, glassbrain_title=None, 
        sample_size=None):   
    """
    Executes an Activation Likelihood Estimation using NiMARE functions with
    recommended thresholding settings (voxel-level p<.001 uncorrected, 
    cluster-level p<.05 FWE-corrected, based on permutation) (Turkeltaub et 
    al., 2002, 2012; Eickhoff et al., 2012). The MNI-152 2mm template is used 
    (Fonov et al., 2009, 2011), coordinates need to be in MNI space.
    Adopted from (see also for references): 
    https://github.com/neurostuff/nimare/blob/b5082c8/nimare/workflows/ale.py#L16
    
    Input:
        data = NiMARE dataset or Sleuth text file
        work_dir, pref = path and prefix for output, defaults to current dir
        vox_thr, n_perm, n_core = recommended ALE settings
        cluster = feature to determine cluster significance, 'mass' for cluster 
            inference based on the sum of ALE values, 'size' for sum of voxels
        output = if true, save output (NiMARE + custom)
        save_cl_info = if true, save cluster info via AtlasReader
        atlases = atlases to determine cluster centers and extents from
        print_glassbrain = if true, plot two glassbrains with unthresholded
            and thresholded ALE results, title = glassbrain_title.
            
    Output:
        NiMARE ALE results object, thresholded volume, binarized volume
    """
    
    # Working directory
    if work_dir == '' or work_dir == None:
        work_dir = os.getcwd()
        
    # Check if input is sleuth textfile and read into NiMARE data base if it is
    # Else assumes data to be a NiMARE dataset objekt
    if isinstance(data, str):
        path,ext = os.path.splitext(data)
        if ext == '.txt':
            ds = ni.io.convert_sleuth_to_dataset(data, target='mni152_2mm')
    else:
        ds = data
        
    # Override sample_sizes in ds if given
    if sample_size is not None:
        ds.metadata = ds.metadata.assign(sample_sizes=sample_size)
        lgr.info(f'Overriding sample sizes with n = {sample_size}.')

    # cluster mass vs extent
    if cluster == 'mass':
        cluster_map = 'logp_desc-mass_level-cluster_corr-FWE_method-montecarlo'
    elif cluster == 'size':
        cluster_map = 'logp_desc-size_level-cluster_corr-FWE_method-montecarlo'
    else:
        lgr.error(f"Valid input for 'cluster' are 'mass' or 'size' not {cluster}")
    
    # get unique study ids (not contrasts)
    ids = ds.metadata.study_id.unique()
    # sum up subjects
    try:
        n = np.sum([ds.metadata[ds.metadata.study_id==iD].sample_sizes.tolist()[0] 
                    for iD in ids]) # unique subjects
    except:
        lgr.error('No sample sizes in dataset and none provided!')
        
    lgr.info(f'Calculating ALE on {len(ids)} experiments with '
             f'{len(ds.ids)} contrasts, {n} subjects and '
             f'{len(ds.coordinates)} foci.')
    lgr.info(f'Thresholding: voxel-level p < {vox_thr} uncorrected, '
             f'cluster-level p < 0.05 FWE-corrected, based on cluster {cluster} '
             f'with {n_perm} permutations.')
        
    # Estimate ALE map
    ds.update_path(work_dir)
    ale = ni.meta.cbma.ale.ALE() # get ALE class
    ale_res = ale.fit(ds) # fit to dataset
    
    # Plot
    if print_glassbrain == True:
        if glassbrain_title is None:
            glassbrain_title1 = 'Unthresholded ALE z-scores'
        else:
            glassbrain_title1 = glassbrain_title+'; unthresholded'
        nl.plotting.plot_glass_brain(
            ale.results.get_map('z'), title=glassbrain_title1, 
            draw_cross=False, cmap='RdBu_r', display_mode='lyrz')
    
    # FWE correction
    corrector = ni.correct.FWECorrector(method='montecarlo', n_iters=n_perm, 
                                        voxel_thresh=vox_thr, n_cores=n_core)
    ale_cres = corrector.transform(ale_res)
    
    # Resulting cluster-level corrected "logp" volumes have the negative 
    # logarithm to the base 10 of the cluster-level FWE-corrected p-value 
    # stored in all voxels of each cluster and require further thresholding 
    # at -log10(0.05) ~= 1.3.
    t = -math.log10(0.05)        
        
    # Save & plot results
    if output == True:
        
        # Save NiMARE results
        ale_cres.save_maps(output_dir=work_dir, prefix=pref)     
            
        # Save clusters and cluster info if significant clusters were found
        if np.max(ale_cres.get_map(cluster_map).
                  get_fdata()) > t: 
            
            # Threshold logp volume
            img_logp_thr = nl.image.threshold_img(
                ale_cres.get_map(cluster_map), 
                threshold=t)
            
            # Save binarized mask volume
            img_bin = nl.image.math_img('img>0', img=img_logp_thr)
            img_bin.to_filename(join(work_dir, pref + 'thresh_bin.nii.gz'))
            
            # Index clusters
            img_cl_idx, img_cl_4d, cl_sizes = index_clusters_bin(img_bin)
            if len(cl_sizes) > 1: # save volumes if more then one cluster 
                # Save volume with clusters indexed from 1 to n(clusters) 
                # and sorted from large to small
                img_cl_idx.to_filename(join(work_dir, 
                                                    pref + 'thresh_idx.nii.gz'))
                # Save 4D volume with one binarized 3D volume for each cluster
                img_cl_4d.to_filename(join(work_dir, pref + 'thresh_4d.nii.gz'))
            
            # Save thresholded ALE stats volume
            img_stat_thr = nl.image.math_img('img1*img2', img1=img_bin, 
                                             img2=ale_cres.get_map('stat'))
            img_stat_thr.to_filename(join(work_dir, pref + 'thresh_stat.nii.gz'))
            
            # Save thresholded z score volume
            img_z_thr = nl.image.math_img('img1*img2', img1=img_bin, 
                                          img2=ale_cres.get_map('z'))
            img_z_thr.to_filename(join(work_dir, pref + 'thresh_z.nii.gz'))
            
            # Extract cluster details from thresholded stat volume, default to 
            # using talairach atlases as does GingerALE
            if save_cl_info == True:
                cl_info, peak_info = get_statmap_info(
                    join(work_dir, pref + 'thresh_stat.nii.gz'), 
                    cluster_extent=0, voxel_thresh=0, direction='pos', 
                    atlas=atlases)
                cl_info.to_csv(join(work_dir, pref + 'thresh_stat_clusters.csv'), 
                               index=False)
                peak_info.to_csv(join(work_dir, pref + 'thresh_stat_peaks.csv'), 
                                 index=False)
            
            # Plot
            if print_glassbrain == True:
                if glassbrain_title is None:
                    glassbrain_title2 = 'Thresholded ALE z-scores'
                else:
                    glassbrain_title2 = glassbrain_title+'; thresholded'
                nl.plotting.plot_glass_brain(img_bin, title=glassbrain_title2, 
                    draw_cross=False, cmap='RdBu_r', display_mode='lyrz')
        
            # Return if volumes saved
            lgr.info('ALE completed, significant clusters found. '
                     f'Results saved to: {work_dir}\n')
            return(ale_cres, img_stat_thr, img_bin)
        else:
            lgr.info('ALE completed, no significant clusters! '
                     f'NiMARE results saved to: {work_dir}\n')
            return(ale_cres, None, None)
    
    # Return if no volumes saved
    else:
        if np.max(ale_cres.get_map(cluster_map).
                  get_fdata()) > t: 
            lgr.info('ALE completed, significant clusters found. No results saved.\n')
        else:
            lgr.info('ALE completed, no significant clusters! No results saved.\n')
        return(ale_cres, None, None)

