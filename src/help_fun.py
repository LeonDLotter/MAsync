#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 11:48:09 2021

@author: Leon D. Lotter
"""
#=============================================================================


def csv_to_nimare_ds(file, exp_var, n_var, con_var, spa_var,
                     x_var, y_var, z_var, sheet=None, single_contrasts=False):     
    """
    Imports data from spreadshead (csv or excel) into a NiMARE dataset. 
    Adopted from Daniel Alcalá López: https://git.bcbl.eu/brainhackdonostia/
    material-2019/-/blob/263fa6c9982d1480ea261abf2cc462eab812501e/fMRI/
    neuroimaging_meta-analysis/nimare_meta-analysis.py
    
    Input:
        file = '/Users/leonlotter/MAsync/MAsync_ALE.xlsx'
        sheet = 'Coordinates'       #-> excel sheet if xls(x) file given
        Column names (= first row) to be defined in input variables, examples:
            exp_var = 'publication' #-> any string (experiment/study name)
            n_var = 'n'             #-> number of subjects in experiment
            con_var = 'contrasts'   #-> any string (contrast name/id)
            spa_var = 'space'       #-> 'MNI' or 'TAL'
            x_var = 'x'             #-> any string, x-coordinate
            y_var = 'y'             #-> any string, y-coordinate
            z_var = 'z'             #-> any string, z-coordinate
        single_contrasts = if True and experiments from the spreadsheet have 
            multiple associated contrasts, these contrasts are stored with 
            separate contrast ids but with the same study id. If False, all 
            coordinates from each study are stored using the corresponding 
            first contrast name while ignoring the existence of multiple 
            contrasts. The 'True' version leads to more liberal findings. 
    If space == 'TAL', coordinates will be converted to MNI space using the
    Brainmap 'tal2mni' Lancaster transform as implemented in NiMARE.
    NiMARE datasets default to MNI152 space with 2x2x2mm voxel resolution.
    
    Output:
        NiMARE dataset object
    """
    
    import pandas as pd
    import numpy as np
    pd.options.mode.chained_assignment = None  # default='warn'
    import nimare as ni
    
    # read table
    if file.endswith(('.xls', '.xlsx')):
        df = pd.read_excel(file, sheet_name=sheet)
    elif file.endswith('.csv'):                    
        df = pd.read_csv(file)
    else:
        print('Unsupported file type!')
        
    # get list with experiments and remove empty entries
    exps = df[exp_var].unique()
    exps = exps[~pd.isnull(exps)]
    
    # loop over experiments
    dset_dict = {}
    for exp in exps:
        # this experiment's data
        this_exp_data = df[df[exp_var] == exp]
        
        # convert TAL coordinates to MNI using brainmap's tal2mni function
        if this_exp_data[spa_var].unique()[0] == 'TAL':
            TAL_coord = this_exp_data[[x_var,y_var,z_var]]
            MNI_coord = ni.utils.tal2mni(TAL_coord)
            this_exp_data.loc[:,[x_var, y_var, z_var]] = MNI_coord
            print('Converting TAL coordinates from {} to MNI.'.format(exp))
            
        # Store multiple contrasts independently per study
        if single_contrasts == True:
            # loop over contrasts
            cons = this_exp_data[con_var].unique()
            contrasts = {}
            for con in cons:
                this_con_data = this_exp_data[this_exp_data[con_var] == con]
                
                # write contrast data to dict
                contrasts[con] = {'coords': { 'space': this_con_data[spa_var].unique()[0],
                                              'x': list(this_con_data[x_var]),
                                              'y': list(this_con_data[y_var]),
                                              'z': list(this_con_data[z_var])},
                                  'metadata': { 'sample_sizes': this_exp_data[n_var].unique()[0]}} 
                                  #'metadata': { 'sample_sizes': [this_exp_data[n_var].unique()[0]]}}
                
            # collect data from each experiment
            dset_dict[exp] = {'contrasts': contrasts }
            
        # Store multiple contrasts as one contrast per study            
        elif single_contrasts == False:
            # write data to dict
            contrast = {'coords': { 'space': 'MNI', #this_exp_data[spa_var].unique()[0],
                                    'x': list(this_exp_data[x_var]),
                                    'y': list(this_exp_data[y_var]),
                                    'z': list(this_exp_data[z_var])},
                        'metadata': { 'sample_sizes': [this_exp_data[n_var].unique()[0]]}}
            
            # collect data from each experiment
            dset_dict[exp] = {'contrasts': {this_exp_data[con_var].unique()[0]: contrast }}
   
    # convert dict to NiMARE dataset, target space MNI152
    ds = ni.dataset.Dataset(dset_dict, target='mni152_2mm')
    
    # output
    if single_contrasts == False:
        print('Concatenating coordinates over multiple contrasts per study.')
    elif single_contrasts == True:
        print('Treating multiple contrasts per study independently.')
    ids = ds.metadata.study_id.unique() # unique study ids (not contrasts)
    n = [ds.metadata[ds.metadata.study_id==iD].sample_sizes.tolist() for iD in ids] # unique subjects
    #n = [ds.metadata[ds.metadata.study_id==iD].sample_sizes.tolist()[0] for iD in ids] # unique subjects
    print('Imported data from {} studies, {} participants, and {} foci as {} unique contrasts into NiMARE dataset.'.format(
          len(ids), np.sum(n), len(ds.coordinates), len(ds.ids)))
    
    return(ds)


#=============================================================================


def nimare_ds_to_sleuth(ds, save_path=None):     
    """
    Convert NiMARE dataset to Sleuth/BrainMap foci text file.
    
    Input:
        ds = NiMARE dataset
        save_path = full save path, defaults to 'current dir/sleuth_foci.txt'
    Output:
        save_path   
    """

    import os
    
    # get experiments and number of subjects
    exps = ds.metadata.study_id
    cons = ds.metadata.contrast_id
    ns = ds.metadata.sample_sizes
    
    # output file
    if save_path is None:
        save_path = os.path.join(os.getcwd(), 'sleuth_foci.txt')
        
    # open text file and write reference space
    with open(save_path, 'w') as t:
        t.write('// Reference=MNI\n')
        
        # write experiments with coordinates
        for e, c, n in zip(exps, cons, ns):
            # experiment info
            t.write('// {}: {} \n// Subjects={} \n'.format(e, c, n[0]))
            # coordinates
            coords = ds.coordinates[ds.coordinates.study_id==e][['x', 'y', 'z']]
            t.write(coords.to_csv(header=False, index=False, sep='\t') + '\n')
    
    print('Sleuth foci file saved to: {}'.format(save_path))
    return(save_path)

    
#=============================================================================


def index_clusters_bin(img):
    """
    Extracts clusters from a binarized volume and assigns indexes in the
    order of cluster size to all voxels within each cluster (1=largest to n=smallest)
    
    Input: img=binarized volume, path or volume
    Output: volume with indexes from 1 to n(clusters), 4D volume, pd.df with sizes
    """
    
    from nilearn.image import math_img, new_img_like, get_data
    from nilearn.regions import connected_regions
    import numpy as np
    import pandas as pd
        
    # load volume, binarize and remove nan (just in case) 
    img = math_img('np.nan_to_num(img) > 0', img=img)
    # get clusters and store in separate volumes within 4D nifti 
    img4D, _ = connected_regions(img, min_region_size=1, extract_type='connected_components')
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
    

def correlate_volumes_via_atlas(x_img, y_img, atlas, method='spearman', colors=None, labels=None):
    """
    Correlates two volumes using roi wise averaged values defined by a given
    parcellation. Prints correlation coefficients, p, and a scatter plot.
    
    Input: 
        img1/2 = volumes to correlate
        atlas  = parcellation
        method = 'spearman' or 'pearson'
        labels = labels to plot over each marker
        colors = colors of each point, must be 1D array with color for each marker
        
    Output: results pd.dataframe
    """
    from nilearn.image import resample_to_img
    from nilearn.input_data import NiftiLabelsMasker
    from scipy.stats import spearmanr, pearsonr
    import matplotlib.pylab as plt    
    from seaborn import regplot
    import pandas as pd

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
        print('Spearman`s r = {}, p = {}'.format(round(r,2),round(p,6)))
    if method == 'pearson':
        r, p = pearsonr(i1_dat, i2_dat)
        print('Pearson`s r = {}, p = {}'.format(round(r,2),round(p,6)))
        
    # plot
    if colors is None:
        regplot(x=i1_dat, y=i2_dat)
    if colors is not None:
        regplot(x=i1_dat, y=i2_dat, color='black', scatter_kws={'facecolors':colors}) 
    if labels is not None:
        for i, l in enumerate(labels):
            plt.text(i1_dat[i], i2_dat[i], l)
            
    return(pd.DataFrame({'idx':list(range(1,len(i1_dat)+1)), 'dat1':i1_dat, 'dat2':i2_dat, 'label':labels}))
    
    
#=============================================================================


def get_size_of_rois(img):
    """
    Extracts a list of roi sizes from a volume or numpy array with indexed clusters
    
    Input: cluster volume or numpy array
    Output: pandas dataframe with sizes
    """
    from nilearn.image import get_data
    import numpy as np
    import pandas as pd
    import nibabel as nb
        
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
    from nilearn.image import get_data, new_img_like
    
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
    print('Kept {} regions with sizes > {} voxels.'.format(len(sizes_go), threshold))
    return(img_drop, sizes_go)
     
    
#=============================================================================


def z_norm_vol(in_file, out_path=None, mask=None):
    """
    Z-normalizes 3D volume to mean and sd of all non-zero voxels (default) 
    or all voxels included in 'mask', if given. Utilizes the 'zscore_normalize' 
    function from the intensity_normalizationpackage: 
        https://intensity-normalization.readthedocs.io/
        
    Input:
        in_file = input volume, nifti object or file path
        out_path = full path to store file, defaults to current path or path
            of input volume (if input volume given as file path)
        mask = binarized mask, if None, determined from non-zero voxels in in_file
        
    Output:
        z normalized volume and mask      
    """

    from nilearn.image import math_img
    from intensity_normalization.normalize.zscore import zscore_normalize
    from os.path import splitext, join
    from os import getcwd
    
    img = math_img('np.nan_to_num(img)', img=in_file) # remove nan, just in case
    
    if mask is None:                            
        mask = math_img('img != 0', img=img)    # create non-zero mask
    img_z = zscore_normalize(img, mask)         # z normalize
    
    if out_path is None and type(in_file) == str: ## save as /(in_file)_z.nii.gz                       
        path,_ = splitext(in_file)              
        img_z.to_filename(path+'_z.nii.gz')     
    elif out_path is None and type(in_file) != str: ## save as cwd/z.nii.gz
        img_z.to_filename(join(getcwd(), 'z.nii.gz'))
    else:
        img_z.to_filename(out_path)             ## save as out_path
        
    return(img_z, mask)
        

#=============================================================================


def get_cluster_peaks(vol_file, z=True, tab_save=False, peak_dist=10,
                        vthr=2.5, cthr=None, cthr_prc=0.1):
    """
    Extracts clusters and peak coordinates based on given voxel and cluster 
    thresholds. Cluster threshold either determined directly (e.g., p < 0.05)
    or via percent of non-zero-voxels (e.g., 0.1 %).
    
    Input:
        vol_file = input volume, nifti object or file path
        z = if True, input volume is normalized to mean and sd of all non-0 voxels
        tab_save = if True, stores table with results at path of input volume,
            only possible if input volume is given as file path
        peak_dist = minimum distance between reported subpeaks in mm
        vthr, cthr = voxel- and cluster-level thresholds
            if cthr = None, cthr is set to cthr_prc % of non-zero voxels
            
    Output: results table and cluster-level threshold
    """
    
    import numpy as np
    from nilearn.image import math_img
    from nilearn.reporting import get_clusters_table        
    from os.path import splitext
    
    vol = math_img('np.nan_to_num(img)', img=vol_file) # remove nan, just in case
    
    if z == True: # z-normalize volume and get in-brain mask
        vol, inbrain = z_norm_vol(vol_file)            
    elif z == False: # get original volume and in-brain mask
        inbrain = math_img('img != 0', img=vol)
        
    if cthr is None:                            
        n_vox = np.sum(inbrain.get_fdata())     # get number of inbrain voxels
        cthr = round(n_vox * cthr_prc / 100)    # cluster threshold based on %
    
    tab = get_clusters_table(stat_img=vol, stat_threshold=vthr, 
                             cluster_threshold=cthr, min_distance=peak_dist) # get cluster table
    if tab_save == True:                        
       path,_ = splitext(vol_file)                   # get path of volume
       if z == True:
           tab.to_csv(path+'_z_clusters.csv')       # save table with suffix
       elif z == False:
           tab.to_csv(path+'_clusters.csv')       # save table with suffix
    else:
       tab.to_csv(tab_save)                     # save table as tab_save
    
    print('Found {} peaks, thresholding: voxel-value > {}; cluster-size > {}'.format(
        len(tab), vthr, cthr))
    print(tab)
    
    return(tab, cthr)
    
    
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
    
    from nilearn.image import get_data, resample_to_img, new_img_like, index_img, load_img
    from nilearn.datasets import load_mni152_template
    import numpy as np
    
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


#=============================================================================


def plot_gb(img, title=None, thresh=0, col='RdBu_r'):
    """
    Plots volume on a nilearn glass brain.
    
    Input: img=input volume, title=plot title, thresh=plot threshold, col=color
    """
    
    from nilearn.plotting import plot_glass_brain
    
    gb = plot_glass_brain(img, title=title, threshold=thresh, draw_cross=False, 
                          cmap=col, display_mode='lyrz')


#=============================================================================


def plot_surface_roi(roi, fig=None, ax=None, cmap='red_transparent', 
                     hemi='right', view='lateral', legend=False):  
    """ 
    Plot ROI on right FSL average sulcus surface.
    
    Input: roi=ROI volume
    """
    
    from nilearn.surface import vol_to_surf
    from nilearn.plotting import plot_surf_roi
    from nilearn.datasets import fetch_surf_fsaverage
    
    fsaverage = fetch_surf_fsaverage()
    
    if hemi == 'right':
        s = plot_surf_roi(fsaverage.pial_right, figure=fig, axes=ax, 
                          roi_map=vol_to_surf(roi, fsaverage.pial_right), 
                          hemi='right', view=view, cmap=cmap, colorbar=legend,   
                          bg_map=fsaverage.sulc_right, bg_on_data=False)
    if hemi == 'left':
        s = plot_surf_roi(fsaverage.pial_left, figure=fig, axes=ax, 
                          roi_map=vol_to_surf(roi, fsaverage.pial_left), 
                          hemi='left', view=view, cmap=cmap, colorbar=legend,   
                          bg_map=fsaverage.sulc_left, bg_on_data=False)
    