#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 11:48:09 2021

@author: Leon D. Lotter
"""

from os import getcwd
from os.path import splitext, join
import numpy as np
import pandas as pd
import nibabel as nb
from nilearn.datasets import load_mni152_template
from nilearn.input_data import NiftiLabelsMasker, NiftiMasker
from nilearn.image import load_img, math_img, new_img_like, get_data, resample_to_img, index_img
from nilearn.regions import connected_regions
from nilearn.reporting import get_clusters_table        
from scipy.stats import zscore
from scipy.stats import spearmanr, pearsonr
import matplotlib.pylab as plt    
from seaborn import regplot

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)


def index_clusters_bin(img):
    """
    Extracts clusters from a binarized volume and assigns indexes in the
    order of cluster size to all voxels within each cluster 
    (1=largest to n=smallest)
    
    Input: img=binarized volume, path or volume
    Output: volume with indexes from 1 to n(clusters), 4D volume, 
            pd.df with sizes
    """
        
    # load volume, binarize and remove nan (just in case) 
    img = math_img('np.nan_to_num(img) > 0', img=img)
    # get clusters and store in separate volumes within 4D nifti 
    img4D, _ = connected_regions(img, min_region_size=1, 
                                 extract_type='connected_components')
    # get 4D data
    dat4D = get_data(img4D)
    
    # get region sizes
    sizes = list()
    for i in range(dat4D.shape[3]): # iterate over fourth dimension of data
        dat3D = dat4D[:,:,:,i] # get 3D array
        sizes.append(len(dat3D[dat3D == 1])) # get and store size
    # sort sizes, dataframe index denotes position in 4D volume
    sizes = pd.DataFrame(sizes, columns=['size'])
    sizes.sort_values('size', ascending=False, inplace=True)
    
    # create new 3D array   
    dat = np.zeros(dat4D.shape[:3])
    dat_4d = np.zeros(dat4D.shape)
    for idx_new, idx_old in enumerate(sizes.index, start=1):
        dat = dat + dat4D[:,:,:,idx_old] * idx_new
        dat_4d[:,:,:,idx_new-1] = dat4D[:,:,:,idx_old]
    # write into new volumes and return
    img_idx = new_img_like(img, dat, copy_header=False)
    img_4d = new_img_like(img4D, dat_4d, copy_header=False)
    return(img_idx, img_4d, sizes)


#=============================================================================
    

def correlate_volumes_via_atlas(x_img, y_img, atlas, method='spearman', 
                                colors=None, labels=None, pr=True, pl=True):
    """
    Correlates two volumes using roi wise averaged values defined by a given
    parcellation. Prints correlation coefficients and a scatter plot.
    
    Input: 
        img1/2 = volumes to correlate
        atlas  = parcellation
        method = 'spearman' or 'pearson'
        labels = labels to plot over each marker
        colors = colors of each point, must be 1D array with color for each 
            marker
        pr, pl = print correlation coefficient, plot scatter plot
        
    Output: correlation coefficient, results pd.dataframe
    """

    # resample volumes to atlas space
    i1 = resample_to_img(x_img, atlas)
    i2 = resample_to_img(y_img, atlas)
    
    # get roi-wise data
    masker = NiftiLabelsMasker(atlas)
    i1_dat = masker.fit_transform(i1)[0] 
    i2_dat = masker.fit_transform(i2)[0] 
    
    # correlate
    if method == 'spearman':
        r, p = spearmanr(i1_dat, i2_dat, axis=1)
        if pr is True: print(f'Spearman`s r = {round(r,2)}') 
    if method == 'pearson':
        r, p = pearsonr(i1_dat, i2_dat)
        if pr is True: print(f'Pearson`s r = {round(r,2)}') 
        
    # plot
    if pl is True:
        if colors is None:
            regplot(x=i1_dat, y=i2_dat)
        if colors is not None:
            regplot(x=i1_dat, y=i2_dat, color='black', 
                    scatter_kws={'facecolors':colors}) 
        if labels is not None:
            for i, l in enumerate(labels):
                plt.text(i1_dat[i], i2_dat[i], l)
            
    return(r, pd.DataFrame({'idx': list(range(1,len(i1_dat)+1)), 
                            'dat1': i1_dat, 
                            'dat2': i2_dat, 
                            'label': labels}))
    
    
#=============================================================================


def get_size_of_rois(img):
    """
    Extracts a list of roi sizes from a volume or numpy array with indexed 
    clusters
    
    Input: cluster volume or numpy array
    Output: pandas dataframe with sizes
    """
        
    if type(img) is nb.nifti1.Nifti1Image:
        dat = get_data(img) # get data matrix
    elif type(img) is np.ndarray:
        dat = img # input is data matrix
    idx = np.unique(dat) # get roi indices
    idx = idx[idx != 0] # drop zero from list
    
    # iterate over roi list anr extract roi sizes
    sizes_list = list()
    for i in idx:
        sizes_list.append(len(dat[dat == i]))
        
    # create df and return
    sizes = pd.DataFrame({'idx':idx, 'size':sizes_list})    
    return(sizes)   


#=============================================================================


def drop_atlasregions_by_size(img, threshold):
    """
    Drops regions from a brain parcellation based on their voxel-size. 
    
    Input: img=parcellation volume, threshold=size threshold
    Output: volume without sub-threshold sized regions
    """
    
    # get sizes of rois in volume and rois to keep after thresholding
    sizes = get_size_of_rois(img)
    sizes_nogo = sizes[sizes['size'] <= threshold]
    sizes_go = sizes[sizes['size'] > threshold]
    
    # get data and set all voxels in regions below threshold to zero
    dat = get_data(img)
    for i in sizes_nogo['idx']:
        dat[dat == i] = 0
    
    # write into new volume and return
    img_drop = new_img_like(img, dat, copy_header=False)
    lgr.info(f'Kept {len(sizes_go)} regions with sizes > {threshold} voxels.')
    return(img_drop, sizes_go)
     
    
#=============================================================================


def z_norm_vol(in_file, out_path=None, mask=None):
    """
    Z-normalizes 3D volume to mean and sd of all non-zero voxels (default) 
    or all voxels included in 'mask', if given. 
        
    Input:
        in_file = input volume, nifti object or file path
        out_path = full path to store file, defaults to current path or path
            of input volume (if input volume given as file path)
        mask = binarized mask, if None, determined from background intensity of 
            in_file
        
    Output:
        z normalized volume 
    """
    
    img = math_img('np.nan_to_num(img)', img=in_file) # remove nan, just in case
    
    masker = NiftiMasker(mask_img=mask, standardize=False, standardize_confounds=False)

    img_data = masker.fit_transform(img)
    img_data_z = zscore(img_data, nan_policy='omit')
    img_z = masker.inverse_transform(img_data_z)
    
    if out_path is None and type(in_file) == str: 
        # save as /(in_file)_z.nii.gz                       
        path,_ = splitext(in_file)              
        img_z.to_filename(path+'_z.nii.gz')   

    elif out_path is None and type(in_file) != str: 
        # save as cwd/z.nii.gz
        img_z.to_filename(join(getcwd(), 'z.nii.gz'))

    else:
        # save as out_path
        img_z.to_filename(out_path) 
        
    return(img_z)
        

#=============================================================================


def get_cluster_peaks(vol_file, z=True, tab_save=False, peak_dist=10,
                        vthr=2.5, cthr=None, cthr_prc=0.1):
    """
    Extracts clusters and peak coordinates based on given voxel and cluster 
    thresholds. Cluster threshold either determined directly (e.g., p < 0.05)
    or via percent of non-zero-voxels (e.g., 0.1 %).
    
    Input:
        vol_file = input volume, nifti object or file path
        z = if True, input volume is normalized to mean and sd of all 
            non-zero voxels
        tab_save = if True, stores table with results at path of input volume,
            only possible if input volume is given as file path
        peak_dist = minimum distance between reported subpeaks in mm
        vthr, cthr = voxel- and cluster-level thresholds
            if cthr = None, cthr is set to cthr_prc % of non-zero voxels
            
    Output: results table and cluster-level threshold
    """
    
    vol = math_img('np.nan_to_num(img)', img=vol_file) # remove nan, just in case
    
    if z == True: # z-normalize volume and get in-brain mask
        vol, inbrain = z_norm_vol(vol_file)            
    elif z == False: # get original volume and in-brain mask
        inbrain = math_img('img != 0', img=vol)
        
    if cthr is None:                            
        n_vox = np.sum(inbrain.get_fdata())     # get number of inbrain voxels
        cthr = round(n_vox * cthr_prc / 100)    # cluster threshold based on %
    
    # get cluster table
    tab = get_clusters_table(stat_img=vol, stat_threshold=vthr, 
                             cluster_threshold=cthr, min_distance=peak_dist) 
    if tab_save == True:                        
       path,_ = splitext(vol_file)                   # get path of volume
       if z == True:
           tab.to_csv(path+'_z_clusters.csv')       # save table with suffix
       elif z == False:
           tab.to_csv(path+'_clusters.csv')       # save table with suffix
    else:
       tab.to_csv(tab_save)                     # save table as tab_save
    
    print(f'Found {len(tab)} peaks, thresholding: voxel-value > {vthr}; '
          f'cluster-size > {cthr}')
    print(tab)
    
    return(tab, cthr)
    
    
#=============================================================================


def parcel_data_to_volume(data, atlas, save_path=None, rank=False):
    """
    Writes parcellated data into nifti volume according to shape and labels 
    of 'atlas'. If rank==True, parcel-values will be ranks of values in data.
    """

    # atlas data
    a = load_img(atlas) # load atlas
    a_dat = a.get_fdata() # get atlas data
    idc = np.unique(a_dat) # get atlas indices
    idc = idc[idc!=0] # remove zero from indices

    # check input data
    if len(idc) != len(data):
        lgr.error('Input array not the same length as atlas has indices!')

    # new volume
    dat = np.zeros(np.shape(a_dat)) # make zero 3d array
    # replace data with ranks of data
    if rank:
        data_ranks = [sorted(data).index(x) + 1 for x in data]
        data = data_ranks
    # create volume
    for i, idx in enumerate(idc):
        dat[a_dat==idx] = data[i] # write data in array
    vol = new_img_like(a, dat)

    # save & return
    if save_path is not None:
        vol.to_filename(save_path)
    return(vol, save_path)


#=============================================================================

def combine_atlases(atlases, renumber=True, target_space=None):
    """
    Combines Parcellations, either by adding them up or renumbering parcel 
    indices from 1 to n(rois) starting with the first input atlas. If parcels
    from different atlases overlap, the nth atlas "overwrites" the n-1th.
    
    Input: 
        atlases = list of atlases
        renumber = if True, renumber parcel indices, if False, add atlases up
        target_space = space to resample atlases to, default to mni152_2mm
    
    Output: new atlas volume
    """
    
    # fetch mni152_2mm as default space
    if target_space is None:
        target_space = load_mni152_template()
           
    dat = np.zeros(target_space.shape) # empty matrix 
    n_rois = 0 # starting point for number of rois
    
    # loop over atlases
    for atlas in atlases:
        # load image data and first frame, if 4D volume
        a = load_img(atlas)
        if len(a.shape) == 4:
            a = index_img(a, 0)
        # resample to first volume and get data
        a = resample_to_img(a, target_space, interpolation='nearest')
        a_dat = get_data(a)
        
        # get parcel indices
        idc = np.unique(a_dat) # get indices
        idc = idc[idc != 0] # remove zero
        
        # loop over indices 
        for i_idx, idx in enumerate(idc, start=n_rois+1):
            if renumber == True: # renumber atlas from n_rois+1 to n_rois+n(idc)
                dat[a_dat == idx] = i_idx # new index
            if renumber == False: # use orginal parcel indices
                dat[a_dat == idx] = idx # original index
        n_rois = n_rois + len(idc) # save current largest index

    # write new volume
    combined = new_img_like(target_space, dat.round(0))
    return(combined)

