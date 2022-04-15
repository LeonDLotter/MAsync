#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:23:24 2021

@author: Leon D. Lotter
"""

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)

import os
import numpy as np
import pandas as pd
import nimare as ni
import nilearn as nl
from .ale import ale

def fsn(ds, cl_true, template=None, fsn_range=None, max_fsn_it=None, 
        quadrant=True, n_subj_null=None, n_coord_null=None, save_dir=None, 
        pref='fsn_', n_perm=1000, v_thr=0.001, cluster='mass', 
        seed=None, plot_gb=False):   
    """
    Computes the 'Fail-Safe-N' (FSN) for ALE results by determining the maximum 
    number of noise studies that can be added to the analysis before the 
    cluster at hand is no longer significant. Python implementation of the 
    method proposed by Acar et al. (2018). To exclude the possibility that
    noise studies contribute to the 'true' cluster, all noise foci are sampled
    outside of the brain quadrant where the cluster's center of mass is located.
    
    Reference: https://doi.org/10.1371/journal.pone.0208177
    Github - Noise study generation: https://github.com/NeuroStat/GenerateNull
    Github - FSN: https://github.com/NeuroStat/FailSafeN
    
    Input:
        ds = 'true' NiMARE dataset from which the original ALE was computed
        cl_true = list of binarized cluster volumes (nifti object or path)
        template = binarized volume to determine possible random foci coordinates
            from. Defaults to MNI152 grey matter probability map thresholded at 0.2
        fsn_range = minimum and maximum number of noise studies to start the
            FSN-finding process. Defaults to [length(true ds), length(true ds)*10]
        max_fsn_it = maximum number of iterations to find exact FSN. Defaults 
            to 10. Final maximum number of ALE runs = 2 + max_fsn_it
        quadrant = if true, samples noise coordinates only in the brain quadrant
            where the current cluster's center of mass is not located
        n_subj_null, n_coord_null = lists of sample sizes and numbers of foci
            to use in noise studies. By default, samples from true ds
        save_dir, pref, n_perm, v_thr, cluster = ALE settings, see 'ale.py'
        seed = starting seed for sampling procedures
        plot_gb = if True, plots glass brain for each ALE.

    Output:
        list of NiMARE results of final ALEs for each cluster, table with FSN
        and NiMARE noise dataset with n = fsn_range[1] noise studies.
    """

    # Set variables ----------------------------------------------------------
    ds.coordinates['space'] = 'MNI' # no effect but avoids warning when merging datasets
    # set save directory 
    if save_dir is None:
        save_dir = os.getcwd()
    # set randum number generator
    if seed is None:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)
    # set minimum and maximum number noise studies to start FSN finding process
    if fsn_range is None:
        m = len(ds.ids)
        M = len(ds.ids) * 10
    else:
        m = fsn_range[0]
        M = fsn_range[1]
    # set maximum number of iterations to determine FSN
    if max_fsn_it is None:
        max_fsn_it = 10
    # get sample sizes for noise studies from 'true' dataset
    if n_subj_null is None:
        n_subj_true = ds.get_metadata(field='sample_sizes')
        n_subj_null = rng.choice(n_subj_true, M)
    # get numbers of foci per study for noise studies from 'true' dataset
    if n_coord_null is None:
        n_coord_true = [len(ds.coordinates.x[ds.coordinates.id==iD]) 
                        for iD in ds.ids] 
        n_coord_null = rng.choice(n_coord_true, M)
      
    # Template --------------------------------------------------------------- 
    # If no template given, fetch MNI152 GM mask and resample to cl_true[0] space
    # produces warning: "pixdim[0] (qfac) should be 1 or -1; setting qfac to 1",
    # the result is okay, however.
    
    if template is None:
        template_1mm = nl.datasets.fetch_icbm152_brain_gm_mask(threshold=0.2)
        template = nl.image.resample_to_img(template_1mm, cl_true[0], 
                                            interpolation='nearest')
        
    # Get possible MNI coordinates for noise studies = all voxels in template
    dat = template.get_fdata() # get matrix
    x, y, z = np.where(dat == 1.0) # get GM coordinates
    affine = template.affine # get affine matrix for MNI transform
    # get mni coordinates
    mni_all = nl.image.coord_transform(x=x, y=y, z=z, affine=affine) 
    mni_all = np.array(mni_all).T
    
    
    # FSN -------------------------------------------------------------------- 
    
    # define function to run ALE and return whether a significant cluster within
    # the location of the 'true' cluster was found
    def fsn_ale(ds_add, cl, prefix, title):
        # run ale
        res, _, mask = ale(data=ds_add, work_dir=save_dir, pref=prefix, 
                           cluster=cluster, vox_thr=v_thr, n_perm=n_perm,
                           save_cl_info=False, 
                           print_glassbrain=plot_gb, glassbrain_title=title)
        if mask is not None: # if 'run_ale' found any significant cluster
            # intersect significant result with 'true' cluster
            intersect_cl = nl.image.math_img('img1 * img2', img1=cl, img2=mask)
            # check whether any above-zero voxels exist within 'true' cluster
            check_cl = np.sum(intersect_cl.get_fdata()) > 0
        elif mask is None: # if 'run_ale' found no significant cluster at all
            check_cl = False
        return(res, check_cl)
    
    # initialize lists to store results
    fsn_tab = list()
    fsn_ress = list()
            
    # loop over clusters
    for cl_i, cl in enumerate(cl_true, start=1):
        lgr.info(f'Cluster {cl_i}:\n') # print cluster number
        cl_pref = pref+'cl'+str(cl_i)+'_' # prefix to store ALE results
        
        
        # 1. Generate noise studies ------------------------------------------
        
        
        # Make sure that noise coordinates do not overlap with the target 
        # cluster by determining a "nogo" brain quadrant
        if quadrant == True:
            lgr.info(f'Generating noise studies while excluding cluster brain quadrant.')
            # Get center of mass of cluster 
            ctr = nl.plotting.find_parcellation_cut_coords(cl)[0]    
            # Get MNI coordinates we can sample from by removing coordinates 
            # belonging to the 'nogo' brain quadrant where 'ctr' is placed in
            q = [ (1 if c > 0 else -1) for c in ctr ]  
            mni_go = mni_all[(mni_all[:,0] * q[0] < 0) | 
                             (mni_all[:,1] * q[1] < 0) | 
                             (mni_all[:,2] * q[2] < 0), :]  
        else:
            lgr.info(f'Generating noise studies while including all brain quadrants.')
            mni_go = mni_all
                
        # Generate noise dataset
        ds_null_dict = {} # empty dict to store noise study-wise results
        for i in range(0, M):  
            # re-initiate RNG with different seed in every loop, just in case
            rng = np.random.default_rng(seed+i) 
            study_id = 'NullExp_' + str(f'{i+1:03}') # study id
            coords_list = rng.choice(mni_go, n_coord_null[i]) # pick random foci
            coords = np.vstack(coords_list) # make coordinates array
            # write into dict as required by NiMARE
            ds_null_dict[study_id] = {
                'contrasts': { '1': {
                    'coords': {
                        'space': 'MNI',
                        'x': list(coords[:, 0]),
                        'y': list(coords[:, 1]),
                        'z': list(coords[:, 2]),
                        },
                    'metadata': { 
                        'sample_sizes': n_subj_null[i]
                        }}}}
        # convert dict to NiMARE dataset
        ds_null = ni.dataset.Dataset(ds_null_dict) 
        ds_null.coordinates['space'] = 'MNI' # avoids error when merging datasets
        ids_null = ds_null.ids # get all noise study IDs
        
        
        # 2. FSN algorithm ---------------------------------------------------
        
        
        # Step 1: run ale with minimum 'm' noise studies added ---------------
        lgr.info(f'FSN step 1: Add minimum {m} noise studies.\n')
        # add first m noise studies to ds
        ds_add = ds.merge(ds_null.slice(ids_null[:m])) 
        fsn_res, check_cl = fsn_ale(ds_add, cl, cl_pref, 
                                    f'FSN, {m} noise studies') # run
        if check_cl == False:
            lgr.info('FSN finished: Cluster not significant, '
                     'indicative of low robustness.\n')
            m_str = '< '+str(m) # true FSN is below minimum
            fsn_tab.append(['cl' + str(cl_i), m_str, m_str, m_str]) # save 
        elif check_cl == True:
            lgr.info('Cluster significant, proceed.\n')
            
            
            # Step 2: run ale with maximum 'M' noise studies -----------------
            lgr.info(f'FSN step 2: Add maximum {M} noise studies.\n')
            # add M noise studies to ds
            ds_add = ds.merge(ds_null.slice(ids_null[:M])) 
            fsn_res, check_cl = fsn_ale(ds_add, cl, cl_pref, 
                                        f'FSN, {M} noise studies') # run          
            if check_cl == True:
                lgr.info('FSN finished: Cluster still significant, '
                         'indicative of highly influencial studies!\n')
                M_str = '> '+str(M) # true FSN is above maximum
                fsn_tab.append(['cl'+str(cl_i), M_str, M_str, M_str]) # save 
            elif check_cl == False:
                lgr.info('Cluster not significant, proceed.\n')
                
                
                # Step 3 to n: run ale with N noise studies ------------------
                m_ = m
                M_ = M
                N = int((m_ + M_) / 2) # start with average of m and M,
                                       # cut to the lower integer
                for ii in range(0, max_fsn_it):
                    lgr.info(f'FSN step 3.{ii}: Add {N} noise studies.\n')
                    # add N noise studies to ds
                    ds_add = ds.merge(ds_null.slice(ids_null[:N])) 
                    fsn_res, check_cl = fsn_ale(ds_add, cl, cl_pref, 
                                                f'FSN, {N} noise studies') # run   
                    if check_cl == False:
                        lgr.info('Cluster not significant.')
                        m_ = m_     # minimum remains the same 
                        M_ = N      # last N is our new maximum
                        N = int((N + m_) / 2) # new N is mean of min and old N
                    elif check_cl == True:
                        lgr.info('Cluster significant.')
                        m_ = N      # last N is our new minimum
                        M_ = M_     # maximum remains the same
                        N = int((N + M_) / 2) # new N is mean of old N and max
                    lgr.info(f'New min FSN = {m_}, new max FSN = {M_}, '
                             f'new N = {N}.')
                    # loop out if FSN exactly determined:
                    # case 1: N==m_==M__
                    # case 2: m_==n & M_==n+1
                        # m_==n -> significant cluster, 
                        # M_==n+1 -> no significant cluster, thus FSN==m_
                    if (m_ == M_) or (M_ == m_+1):
                        lgr.info(f'Looping out after {ii+1} iterations, '
                                 'FSN exactly determined for the given set of '
                                 'noise studies.\n')
                        M_ = m_ 
                        break 
                        
                    
                # Finished: Show and store results ---------------------------
                # If enough iterations where performed, FSN, min and max are 
                # the same. If not enough iterations, the true FSN lies 
                # between min and max.
                lgr.info(f'Cluster {cl_i}: FSN finished: FSN (~)= {N} '
                         f'(min = {m_}, max = {M_})\n')
                fsn_tab.append(['cl'+str(cl_i), N, m_, M_]) # save FSN values
     
        # store ALE results
        fsn_ress.append(fsn_res) # save last ALE result
    
    # write list to dataframe
    fsn_tab = pd.DataFrame(fsn_tab, columns=['cluster', 'appr_fsn', 
                                             'min_fsn', 'max_fsn'])
    # return
    return(fsn_tab, fsn_ress, ds_null)