#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:22:25 2021

@author: Leon D. Lotter
"""

import os
import numpy as np
import pandas as pd
from nilearn.input_data import NiftiMasker
from nilearn.masking import intersect_masks
from nilearn.plotting import plot_glass_brain
from .ale import ale

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)


def loeo(ds, v_thr=0.001, n_perm=1000, n_core=-1, 
         save_dir=None, prefix='LOEO_', conj_save=None):    
    """
    Calculates N(experiments) ALEs from NiMARE dataset while excluding 
    one experiment at a time. The intersection of all masks will be saved.
    
    Input:
        ds = NiMARE dataset
        v_thr, n_perm, n_core = ALE settings
        save_dir, prefix = directory and file name prefix for ALE results
        conj_save = path for conjunction volume, defaults to current directory
          
    Output: 
        list of NiMARE ALE results & intersection volume.
    """
    
    if save_dir is None:
        save_dir = os.getcwd()
    if conj_save is None:
        conj_save = os.path.join(os.getcwd(), 'LOEO_conjunction.nii.gz')

    # get experiments and iterate over them
    exps = ds.ids
    ress = list()
    masks = list()
    for i, exp in enumerate(exps):
        lgr.info(f'\nCalculate ALE without experiment: {i+1} {exp}')
        # remove current experiment from ds
        ds_x = ds.slice(np.delete(exps, i))
        # run ALE
        res, _, mask = ale(data=ds_x, work_dir=save_dir, pref=prefix+str(i+1)+'_', 
                           vox_thr=v_thr, n_perm=n_perm, n_core=n_core, 
                           print_glassbrain=False, save_cl_info=False)
        # save NiMARE results and cluster masks
        ress.append(res)
        masks.append(mask)

    # calculate and print conjunction of mask volumes
    conjunction = intersect_masks(masks, 1, connected=False) # calc intersection
    conjunction.to_filename(conj_save) # save volume
    lgr.info(f'Conjunction volume saved to: {conj_save}')
    plot_glass_brain(conjunction, title='Conjunction of thresholded LOEO volumes', 
                     threshold=0, draw_cross=False, cmap='RdBu_r', 
                     display_mode='lyrz')
    
    # return
    return(ress, conjunction)

        
#=============================================================================


def contributions(ale_ALL, ale_EXP, cl_masks, exp_ids):  
    """
    Calculates the average contribution of each experiment's ALE values to 
    the main ALE result within given cluster(s). 
    
    Input:
        ale_ALL = ale value volume when all experiments are included
        ale_EXP = list of ale volumes when excluding one experiment at a time
        cl_masks = list of binary cluster volumes
        exp_ids = list of experiment IDs
    Volumes can be nifti1 objects or file paths.
    
    Output:
        pandas dataframe with results for all clusters.
    """
    
    # Initialize dataframe with experiment IDs
    contributions = pd.DataFrame({'ID': exp_ids})
         
    for i, cl in enumerate(cl_masks, start=1):
        # Set Nilearn masker
        masker = NiftiMasker(mask_img=cl)
        
        # Loop over experiments
        ale_means = list()
        for exp in ale_EXP:
            # get ale values within cluster
            exp_ale_values = masker.fit_transform(exp) 
            # mean of within-cluster values
            exp_ale_mean = np.mean(exp_ale_values) 
            # collect results
            ale_means.append(exp_ale_mean) 
        
        # Main ALE values within cluster 
        main_ale_values = masker.fit_transform(ale_ALL)
        main_ale_mean = np.mean(main_ale_values)
        
        # Save results
        contributions_cl = pd.DataFrame(
            {'cl'+str(i)+'_exp_sum': ale_means,
             'cl'+str(i)+'_all_sum': main_ale_mean,
             'cl'+str(i)+'_contr_prc': (1-ale_means/main_ale_mean)*100,
            })
        contributions = pd.concat([contributions, contributions_cl], axis=1)
    
    lgr.info('Calculated experiment-wise contributions of '
             f'{len(exp_ids)} experiments to {len(cl_masks)} clusters.')
    
    return(contributions)