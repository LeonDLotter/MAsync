#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 15:02:53 2021

@author: Leon D. Lotter
"""

import pandas as pd
import numpy as np
from atlasreader import get_statmap_info
from scipy.spatial import KDTree
from statsmodels.stats.multitest import fdrcorrection
from nimare.utils import mm2vox, vox2mm
from nimare.stats import null_to_p
from nilearn.image import load_img, get_data, new_img_like, math_img
from tqdm.auto import tqdm

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)

def fnirs_result(data, atlas, labels, publications=None):
    
    """
    Evaluates channel-wise MNI coordinates from fNIRS studies by a given whole-
    brain parcellation. For each coordinate, the corresponding parcel is found.
    For each parcel, a summary of the number of "INS"-coordinates, "non-INS"-
    coordinates and associated study information is created. Additionally, the
    parcel described anatomically  via atlasreader using the AAL atlas.
    If a coordinate is not located within the atlas, it is replaced by the 
    nearest coordinate within the atlas using scipy.spatial.KDTree.
    
    Input:
        data = dataframe with the following columns: 'publication': study name, 
            'n': number of subjects, 'sync': binary indicator of study 
            result, 'x', 'y', 'z': MNI coordinates
        atlas = parcellation volume
        labels = dataframe with atlas labels corresponding to atlas indices
            with following columns: 'idx': parcel indices, 'label': parcel name
        publications = list of publications to include, if None -> all
        
    Output:
        Dataframe with parcel-wise results
        Dict of nifti volumes, 1: number of INS channels per parcel, 2: ratio of 
            INS to all channels per parcel, 3: ratio of INS to all channels per 
            parcel multiplied by all subjects that "contributed" to the parcel as 
            an (arbitrary) index of convergence 
        Dataframe with input study information extended by nearest coordinates
            (real-world, 'i_','j_','k_') and distance ('dist') to original 
            coordinates ('i','j','k'). If distance == 0, nearest == original 
    """

    # fNIRS DATA
    if type(data)==str:
        data = pd.read_csv(data)
    df = data[['publication', 'n', 'sync', 'x', 'y', 'z']]

    # USER-DEFINED publications TO EXCLUDE?
    if publications is not None:
        df = df[df.publication.isin(publications)]

    # EXCLUDE ROWS WITH MISSING COORDINATE DATA
    df = df[~df['x'].isna()]
    lgr.info(f'Loaded fNIRS data, kept {len(df)} rows with coordinate information '
             f'from {len(df.publication.unique())} experiments.')
    # check that relevant cols are numeric
    df[['n', 'sync', 'x', 'y', 'z']] = df[['n', 'sync', 'x', 'y', 'z']].apply(pd.to_numeric, axis = 1)
    
    # ATLAS 
    atlas = load_img(atlas)
    a_dat = get_data(atlas) # data matrix
    
    # TRANSFORM COORDINATES TO WORLD SPACE
    df.loc[:, ['i', 'j', 'k']] = mm2vox(df[['x', 'y', 'z']], atlas.affine)
    
    # FIND NEAREST NEIGHBOR WITHIN ATLAS IF FNIRS COORDINATE NOT IN ATLAS
    ijk = np.argwhere(a_dat > 0) # get indices of non-zero atlas points
    kdtree = KDTree(ijk) # prepare KDTree class
    dist, points = kdtree.query(df[['i', 'j', 'k']], 1) # find nearest neighbors
    df.loc[:, ['i_', 'j_', 'k_']] = ijk[points] # write coordinates to df
    df.loc[:, 'dist'] = dist # write distance to df; if dist == 0, then i,j,k == i_,j_,k_
    lgr.info(f'Finding nearest atlas region for each coordinate. Maximum distance = {df.dist.max():1.2f}.')

    # GET REGION ASSOCIATED WITH COORDINATES
    regions = []
    for _, row in df.iterrows(): # iterate over fNIRS coordinates
        idx = a_dat[(row['i_'], row['j_'], row['k_'])] # get atlas index of coordinate
        label = labels['label'][labels['idx']==idx].values[0] # get atlas label
        regions.append([idx, label]) # save
    df.loc[:, ['atlas_idx', 'atlas_region']] = regions # append to fnirs dataframe

    # GET FNIRS RESULTS
    lgr.info('Extracting parcel-wise results...')
    result = []
    for region_idx in df['atlas_idx'].unique():

        # characterize atlas region anatomically using anatomical atlas, default to AAL
        region_info, _ = get_statmap_info(math_img(f'a=={region_idx}', a=atlas), 1, ['aal', 'talairach_ba'], 1)

        # get all df rows with coordinate in atlas region
        df_region = df[df['atlas_idx']==region_idx]
        
        # get counts
        # channels
        n_ch_all = len(df_region)
        n_ch_sync = df_region['sync'].sum()
        # experiments
        exp_all = df_region['publication'].unique()
        exp_sync = df_region['publication'][df_region['sync']==1].unique()
        n_exp_all = len(exp_all)
        n_exp_sync = len(exp_sync)
        # subjects
        n_sub_all = np.nansum([df_region['n'][df_region['publication']==pub].values[0] for pub in exp_all])
        n_sub_sync = np.nansum([df_region['n'][df_region['publication']==pub].values[0] for pub in exp_sync])
        
        # get ratios
        ch_ratio = n_ch_sync / n_ch_all
        ch_ratio_ch_sync = ch_ratio * n_ch_sync
        ch_ratio_sub_all = ch_ratio * n_sub_all
        ch_ratio_exp_all = ch_ratio * n_exp_all

        # save result
        r = dict(
            region = labels['label'][labels['idx']==region_idx].values[0], 
            region_idx = region_idx,
            n_ch_all = n_ch_all,
            n_ch_sync = n_ch_sync,
            n_exp_all = n_exp_all,
            n_exp_sync = n_exp_sync,
            n_sub_all = n_sub_all,
            n_sub_sync = n_sub_sync,
            ch_ratio = ch_ratio,
            ch_ratio_ch_sync = ch_ratio_ch_sync,
            ch_ratio_sub_all = ch_ratio_sub_all,
            ch_ratio_exp_all = ch_ratio_exp_all,
            exp_all = ', '.join(exp_all),
            exp_sync = ', '.join(exp_sync),
            AAL = region_info['aal'][0],
            BA = region_info['talairach_ba'][0]
        )
        # append to list
        result.append(r) 

    # convert to dataframe, sort, save
    result_df = pd.DataFrame(result)
    result_df.sort_values('ch_ratio_sub_all', ascending=False, inplace=True)
    result_df.reset_index(drop=True)
         
    # NEW VOLUMES
    lgr.info('Creating volumes...')
    # new data matrices
    dat_sync, dat_ratio_sub, dat_ratio_exp = np.zeros(a_dat.shape), np.zeros(a_dat.shape), np.zeros(a_dat.shape) # empty data matrices               
    for _, row in result_df.iterrows():
        s = row['n_ch_sync']
        if s > 0:
            dat_sync[a_dat==row['region_idx']] = s # number of significant channels
            dat_ratio_sub[a_dat==row['region_idx']] = row['ch_ratio_sub_all'] # ratio of sig to all channels 
            dat_ratio_exp[a_dat==row['region_idx']] = row['ch_ratio_exp_all'] # ratio of sig to all channels weighted by total subjects
        else: # fill parcels with -1 to indicate covered areas 
            dat_sync[a_dat==row['region_idx']] = -1
            dat_ratio_sub[a_dat==row['region_idx']] = -1
            dat_ratio_exp[a_dat==row['region_idx']] = -1
    # create dict with volumes    
    vols = {'n_ch_sync': new_img_like(atlas, dat_sync),
            'n_ch_ratio_sub_all': new_img_like(atlas, dat_ratio_sub),
            'n_ch_ratio_exp_all': new_img_like(atlas, dat_ratio_exp)}    
    
    lgr.info('Finished.')
    
    return(result_df, vols, df)


#=============================================================================


def fnirs_permute(fnirs_coord_df, region_indices, n_perm):

    """
    Calculated permutation-based p-values for the fnirs result as estimated in
    the function above.
    
    Input:
        fnirs_coord_df = third output table form fnirs_result function
        region_indices = region_indices to look at, must match {region_idx} in
            fnirs_coord_df
        n_perm = number of permutations
        
    Output:
        Three dictionaries with (1) permutated data, (2) p values, (3) q values,
        in each dict there are three arrays stored for (1) number of sync channels,
        (2) sync/all channel ratio weighted by subject number, and (3) weighted by 
        experiment number.
    """

    # function to extract results ---------------------------------------------------------------------------------------
    def get_fnirs_indices(df, region_idc):

        # empty list to store results
        ch_sync, ch_ratio_sub_all, ch_ratio_exp_all = [], [], []

        # collect parcel-wise data (iterating parcels in the order of the original dataframe)
        for region_idx in region_idc:
            # all df rows with coordinate in atlas region
            df_region = df[df['atlas_idx']==region_idx]

            # get relevant counts
            # channels
            n_ch_all = len(df_region)
            n_ch_sync = df_region['sync'].sum()
            ch_ratio = n_ch_sync / n_ch_all
            # experiments
            exp_all = df_region['publication'].unique()
            #exp_sync = df_region['publication'][df_region['sync']==1].unique()
            n_exp_all = len(exp_all)
            #n_exp_sync = len(exp_sync)
            # subjects
            n_sub_all = np.nansum([df_region['n'][df_region['publication']==pub].values[0] for pub in exp_all])
            #n_sub_sync = np.nansum([df_region['n'][df_region['publication']==pub].values[0] for pub in exp_sync])
            
            # collect results
            ch_sync.append(n_ch_sync)
            ch_ratio_sub_all.append(ch_ratio * n_sub_all)
            ch_ratio_exp_all.append(ch_ratio * n_exp_all)

        return ch_sync, ch_ratio_sub_all, ch_ratio_exp_all
    # ---------------------------------------------------------------------------------------

    # get real result
    real_n_ch_sync, real_ch_ratio_sub_all, real_ch_ratio_exp_all = \
        get_fnirs_indices(fnirs_coord_df, region_indices)

    # get permuted results
    perm_n_ch_sync = np.zeros((len(region_indices), n_perm))
    perm_ch_ratio_sub_all = np.zeros((len(region_indices), n_perm))
    perm_ch_ratio_exp_all = np.zeros((len(region_indices), n_perm))

    for i in tqdm(range(n_perm)):
        # copy dataframe
        df = fnirs_coord_df.copy(deep=True)
        # drop parcel labels (dont fit anymore)
        df.drop('atlas_region', axis=1)
        # randomize coordinate-parcel-assignment
        df['atlas_idx'] = np.random.permutation(df['atlas_idx'].values)
        
        # save results
        perm_n_ch_sync[:,i], perm_ch_ratio_sub_all[:,i], perm_ch_ratio_exp_all[:,i] = \
            get_fnirs_indices(df, region_indices)

    # make dfs
    perm_n_ch_sync = pd.DataFrame(data=perm_n_ch_sync, index=region_indices)
    perm_ch_ratio_sub_all = pd.DataFrame(data=perm_ch_ratio_sub_all, index=region_indices)
    perm_ch_ratio_exp_all = pd.DataFrame(data=perm_ch_ratio_exp_all, index=region_indices)

    # get p values
    p_n_ch_sync, p_ch_ratio_sub_all, p_ch_ratio_exp_all = [], [], []
    for i, region_idx in enumerate(region_indices):
        p_n_ch_sync.append(null_to_p(real_n_ch_sync[i], perm_n_ch_sync.loc[region_idx,:], tail='upper'))
        p_ch_ratio_sub_all.append(null_to_p(real_ch_ratio_sub_all[i], perm_ch_ratio_sub_all.loc[region_idx,:], tail='upper'))
        p_ch_ratio_exp_all.append(null_to_p(real_ch_ratio_exp_all[i], perm_ch_ratio_exp_all.loc[region_idx,:], tail='upper'))

    # fdr correction
    _, q_n_ch_sync = fdrcorrection(p_n_ch_sync)
    _, q_ch_ratio_sub_all = fdrcorrection(p_ch_ratio_sub_all)
    _, q_ch_ratio_exp_all = fdrcorrection(p_ch_ratio_exp_all)

    return(
       dict(
           perm_n_ch_sync = perm_n_ch_sync,
           perm_ch_ratio_sub_all = perm_ch_ratio_sub_all,
           perm_ch_ratio_exp_all = perm_ch_ratio_exp_all
       ),
       dict(
           p_n_ch_sync = p_n_ch_sync, 
           p_ch_ratio_sub_all = p_ch_ratio_sub_all,
           p_ch_ratio_exp_all = p_ch_ratio_exp_all
       ),
       dict(
           q_n_ch_sync = q_n_ch_sync,
           q_ch_ratio_sub_all = q_ch_ratio_sub_all,
           q_ch_ratio_exp_all = q_ch_ratio_exp_all
       ),
    )


#=============================================================================


def rand_nimare_coords(ds, mask, ids=None, radius=10, seed=None, verbose=True):
    """
    Takes a dataset {ds} and randomizes coordinates of studies defined by {ids} 
    within a sphere with a radius of {radius} * voxelsize mm. If ids is None, all coordinates are
    randomized. Dataset and mask must be in the same space!
    """
    
    # copy dataset to leave the input original 
    ds = ds.copy()

    # set seed
    if seed: 
        np.random.seed(seed)

    # get study ids
    if ids is None:
        ids = ds.ids
    
    if verbose: lgr.info(f'Randomizing coordinates within a sphere of {radius} mm.')
    # get possible coordinates to sample from
    m = load_img(mask)
    m_dat = m.get_fdata()
    # all possible coordinates
    ijk_all = np.argwhere(m_dat > 0) # array with all non zero mask coordinates
    kdtree = KDTree(ijk_all) # prepare KDTree class

    # iterate over ids
    for id in ids:
        if verbose: lgr.info(f'Randomizing coordinates from {id}.')
        # get MNI coordinates
        xyz = ds.coordinates[['x', 'y', 'z']][ds.coordinates.id==id]
        # convert to world coordinates
        ijk = mm2vox(xyz, m.affine)
        # get all coordinates in sphere around coordinate
        _, points = kdtree.query(ijk, k=None, distance_upper_bound=radius)
        # loop over coordinates 
        xyz_ = list()
        for p in range(len(points)):
            # get random coordinate
            rand_coord = ijk_all[np.random.choice(points[p])]
            # convert back to mni coordinate
            xyz_.append(vox2mm(rand_coord, m.affine))
        # write back to ds
        ds.coordinates.loc[ds.coordinates.id==id, ['x', 'y', 'z']] = xyz_
    
    return(ds)
