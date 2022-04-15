#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:28:08 2021

@author: Leon D. Lotter
"""

import os
import pandas as pd
import numpy as np
import nimare as ni
from scipy.spatial import KDTree
from nimare.utils import mm2vox, vox2mm
from nilearn.image import load_img

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)


def csv_to_nimare_ds(file, exp_var, n_var, con_var, spa_var,
                     x_var, y_var, z_var, sheet=None, single_contrasts=False):     
    """
    Imports data from spreadshead (csv or excel, can be loaded) into a NiMARE dataset. 
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

    pd.options.mode.chained_assignment = None  # default='warn'
    
    # read table
    if not isinstance(file, pd.DataFrame):
        if file.endswith(('.xls', '.xlsx')):
            df = pd.read_excel(file, sheet_name=sheet)
        elif file.endswith('.csv'):                    
            df = pd.read_csv(file)
        else:
            lgr.error('Unsupported file type!')
    else:
        df = file
        
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
            lgr.info(f'Converting TAL coordinates from {exp} to MNI.')
            
        # Store multiple contrasts independently per study
        if single_contrasts == True:
            # loop over contrasts
            cons = this_exp_data[con_var].unique()
            contrasts = {}
            for con in cons:
                this_con_data = this_exp_data[this_exp_data[con_var] == con]
                
                # write contrast data to dict
                contrasts[con] = {'coords': { 'space': this_con_data[spa_var].
                                             unique()[0],
                                              'x': list(this_con_data[x_var]),
                                              'y': list(this_con_data[y_var]),
                                              'z': list(this_con_data[z_var])},
                                  'metadata': { 'sample_sizes': this_exp_data[n_var].
                                               unique()[0]}} 
                
            # collect data from each experiment
            dset_dict[exp] = {'contrasts': contrasts }
            
        # Store multiple contrasts as one contrast per study            
        elif single_contrasts == False:
            # write data to dict
            contrast = {'coords': { 'space': 'MNI', 
                                    'x': list(this_exp_data[x_var]),
                                    'y': list(this_exp_data[y_var]),
                                    'z': list(this_exp_data[z_var])},
                        'metadata': { 'sample_sizes': [this_exp_data[n_var].
                                                       unique()[0]]}}
            
            # collect data from each experiment
            dset_dict[exp] = {'contrasts': {
                this_exp_data[con_var].unique()[0]: contrast }}
   
    # convert dict to NiMARE dataset, target space MNI152
    ds = ni.dataset.Dataset(dset_dict, target='mni152_2mm')
    
    # output
    if single_contrasts == False:
        lgr.info('Concatenating coordinates over multiple contrasts per study.')
    elif single_contrasts == True:
        lgr.info('Treating multiple contrasts per study independently.')
    # unique study ids (not contrasts)
    ids = ds.metadata.study_id.unique() 
    # unique subjects
    n = [ds.metadata[ds.metadata.study_id==iD].sample_sizes.tolist() for iD in ids] 
    lgr.info(f'Imported data from {len(ids)} studies, {np.sum(n)} participants, ' 
             f'and {len(ds.coordinates)} foci as {len(ds.ids)} unique contrasts ' 
             'into NiMARE dataset.')
    
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
    
    lgr.info(f'Sleuth foci file saved to: {save_path}')
    return(save_path)


#=============================================================================


def fnirs_to_nimare_dataset(fnirs_file, rand_atlasviewer=False, radius=5, mask_file=None):
    """
    Import fNIRS-coordinates-table into NiMARE dataset. If rand_atlasviewer is 
    True, randomize coordinates with "coordinate_source"=="AtlasViewer" within
    a sphere with {radius} within non-zero voxels of {mask_file} 
    """

    # get data & clean
    if not isinstance(fnirs_file, pd.DataFrame):
        df = pd.read_csv(fnirs_file)
    else:
        df = fnirs_file
    df = df[['publication', 'n', 'contrast', 'coordinate source', 'sync', 'x', 'y', 'z']]
    df = df[~df.x.isna()]
    df[['n', 'sync', 'x', 'y', 'z']] = df[['n', 'sync', 'x', 'y', 'z']].apply(pd.to_numeric, axis = 1)
    df.rename(columns={'coordinate source': 'source'}, inplace=True)
    
    # fill empty study info rows
    df_clean = list()
    for pub in df.publication.unique():
        # get publication rows
        pub_df = df[df.publication==pub].reset_index(drop=True)
        # fill study info rows
        pub_df = pub_df.assign(n=pub_df.n[0], 
                               contrast=pub_df.contrast[0],
                               source=pub_df.source[0])
        df_clean.append(pub_df)
    df_clean = pd.concat(df_clean)
    lgr.info(f'Importing {len(df_clean.publication.unique())} fNIRS experiments.')
 
    # only sync coordinates
    df_sync = df_clean[df_clean.sync==1].reset_index(drop=True)

    # randomize atlasviewer coordinates
    if rand_atlasviewer:
        lgr.info(f'Randomizing AtlasViewer coordinates within a sphere of {radius} mm.')
        # get possible coordinates to sample from
        mask = load_img(mask_file)
        m_dat = mask.get_fdata()
        # all possible coordinates
        ijk_all = np.argwhere(m_dat > 0) # array with all non zero mask coordinates
        kdtree = KDTree(ijk_all) # prepare KDTree class

        # iterate over fNIRS coordinates
        for i, row in df_sync.iterrows(): 
            if row.source == 'AtlasViewer':
                # get world coordinates
                ijk = mm2vox(row[['x', 'y', 'z']], mask.affine)
                # get all coordinates in sphere around coordinate
                _, points = kdtree.query(ijk, k=None, distance_upper_bound=radius) 
                # get random coordinate
                random_coordinate = ijk_all[np.random.choice(points)]
                # get mni coordinate
                x_, y_, z_ = vox2mm(random_coordinate, mask.affine)
                df_sync.loc[i, ['x', 'y', 'z']] = [x_, y_, z_]
    
    # make temporary output table
    df_fnirs = pd.DataFrame(df_sync[['publication', 'n', 'contrast', 'x', 'y', 'z']])
    df_fnirs['space'] = 'MNI'

    # get dataset
    ds = csv_to_nimare_ds(file=df_fnirs, exp_var='publication', n_var='n', 
                          con_var='contrast', spa_var='space', 
                          x_var  ='x', y_var='y', z_var='z',
                          single_contrasts=False) 

    return(ds)