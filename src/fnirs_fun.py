#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 15:02:53 2021

@author: Leon D. Lotter
"""

def fnirs_result(data, atlas, labels, anat_atlas='aal'):
    
    """
    Evaluates channel-wise MNI coordinates from fNIRS studies by a given whole-
    brain parcellation. For each coordinate, the corresponding parcel is found.
    For each parcel, a summary of the number of "INS"-coordinates, "non-INS"-
    coordinates and associated study information is created. Additionally, the
    parcel is anatomically described via atlasreader and a given anatomical atlas.
    If a coordinate is not located within the atlas, it is replaced by the 
    nearest coordinate within the atlas using scipy.spatial.KDTree.
    
    Input:
        data = dataframe with the following columns: 'publication': study name, 
            'n': number of subjects, 'synchrony': binary indicator of study 
            result, 'x', 'y', 'z': MNI coordinates
        atlas = parcellation volume
        labels = dataframe with atlas labels corresponding to atlas indices
            with following columns: 'idx': parcel indices, 'label': parcel name
        anat_atlas = anatomical atlas, see atlasreader documentation
        
    Output:
        Dataframe with parcel-wise results
        List of nifti volumes, 1: number of INS channels per parcel, 2: ratio 
            of INS to all channels per parcel multiplied by all subjects that
            "contributed" to the parcel as an (arbitrary) index of convergence 
        Dataframe with input study information extended by nearest coordinates
            (real-world, 'i_','j_','k_') and distance ('dist') to original 
            coordinates ('i','j','k'). If distance == 0, nearest == original 
    """
    
    import pandas as pd
    import numpy as np
    from atlasreader import get_statmap_info
    from scipy.spatial import KDTree
    from nimare.utils import mm2vox
    from nilearn.image import load_img, get_data, new_img_like, math_img
    
    # fNIRS DATA
    df = data[['publication', 'n', 'synchrony', 'x', 'y', 'z']]
    
    # ATLAS 
    atlas = load_img(atlas)
    a_dat = get_data(atlas) # data matrix
    
    # TRANSFORM COORDINATES TO WORLD SPACE
    df.loc[:, ['i', 'j', 'k']] = mm2vox(df[['x', 'y', 'z']], atlas.affine)
    
    # FIND NEAREST NEIGHBOR WITHIN ATLAS IF FNIRS COORDINATE NOT IN ATLAS
    i, j, k = np.where(a_dat > 0) # get indices of non-zero atlas points
    ijk = np.array((i, j, k)).T # array with all indices
    kdtree = KDTree(ijk) # prepare KDTree class
    dist, points = kdtree.query(df[['i', 'j', 'k']], 1) # find nearest neighbors
    df.loc[:, ['i_', 'j_', 'k_']] = ijk[points] # write coordinates to df
    df.loc[:, 'dist'] = dist # write distance to df; if dist == 0, then i,j,k == i_,j_,k_
    
    # GET REGION ASSOCIATED WITH COORDINATES
    regions = []
    for _, row in df.iterrows(): # iterate over fNIRS coordinates
        idx = a_dat[(row['i_'], row['j_'], row['k_'])] # get atlas index of coordinate
        label = labels['label'][labels['idx']==idx].values[0] # get atlas label
        regions.append([idx, label]) # save
    df.loc[:, ['atlas_idx', 'atlas_region']] = regions # append to fnirs dataframe
    
    # GET FNIRS RESULTS
    result = []
    for region_idx in df['atlas_idx'].unique():
        # characterize atlas region anatomically using anatomical atlas, default to AAL
        region_info, _ = get_statmap_info(math_img('a=={}'.format(region_idx), a=atlas), 1, anat_atlas, 1)
        # all df rows with coordinate in atlas region
        df_region = df[df['atlas_idx']==region_idx]
        
        # synchrony result
        n_sync = sum(df_region['synchrony'])
        n_all = len(df_region)
        ratio = round(n_sync / n_all, 2)
        # experiments
        exps_sync = df_region['publication'][df_region['synchrony']==1].unique()
        exps_all = df_region['publication'].unique()
        # subjects
        n_sub_sync = np.nansum(df[df['publication'].isin(exps_sync)]['n'])
        n_sub_all = np.nansum(df[df['publication'].isin(exps_all)]['n'])
        
        # store in dict
        r = {'region': labels['label'][labels['idx']==region_idx].values[0], 
             'region_idx': region_idx, 
             'n_chan_sync': n_sync, 
             'n_chan_all': n_all, 
             'chan_sync/all': ratio,
             'chan_sync/all*chan_sig': ratio * n_sync,
             'chan_sync/all*sub_all': ratio * n_sub_all,
             'chan_sync/all*exps_all': ratio * len(exps_all),
             'n_exps_sync': len(exps_sync),
             'n_exps_all': len(exps_all),
             'n_sub_sync': n_sub_sync,
             'n_sub_all': n_sub_all,
             'exps_sync':  ', '.join(exps_sync),
             'exps_all':  ', '.join(exps_all),
             'AAL': region_info['aal'][0]
            }
        result.append(r) # append to list
    
    # convert to dataframe, sort, save
    result_df = pd.DataFrame(result)
    result_df.sort_values('chan_sync/all*sub_all', ascending=False, inplace=True)
                       
    # NEW VOLUMES
    # new data matrices
    dat_sync, dat_ratio = np.zeros(a_dat.shape), np.zeros(a_dat.shape) # empty data matrices               
    for _, row in result_df.iterrows():
        s = row['n_chan_sync']
        if s > 0:
            dat_sync[a_dat==row['region_idx']] = s # number of significant channels
            dat_ratio[a_dat==row['region_idx']] = row['chan_sync/all*sub_all'] # ratio of sig to all channels weighted by total subjects
        else: # fill parcels with values > 0 for vizualisation
            dat_sync[a_dat==row['region_idx']] = 0.0001 
            dat_ratio[a_dat==row['region_idx']] = 1
    # create volumes        
    vol_sync = new_img_like(atlas, dat_sync)
    vol_ratio = new_img_like(atlas, dat_ratio)
    
    return(result_df, [vol_sync, vol_ratio], df)
