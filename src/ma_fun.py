#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 16:35:19 2021

@author: Leon D. Lotter
"""
#=============================================================================


def run_ale(data, work_dir=None, pref='ALE_', 
            vox_thr=0.001, n_perm=1000, n_core=-1, output=True, 
            save_cl_info=True, atlases=['talairach_gyrus', 'talairach_ba', 'aal'],
            print_glassbrain=True, glassbrain_title=None):   
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
        output = if true, save output (NiMARE + custom)
        save_cl_info = if true, save cluster info via AtlasReader
        atlases = atlases to determine cluster centers and extents from
        print_glassbrain = if true, plot two glassbrains with unthresholded
            and thresholded ALE results, title = glassbrain_title.
            
    Output:
        NiMARE ALE results object, thresholded volume, binarized volume
    """
    
    import os
    import numpy as np
    import nimare as ni
    import nilearn as nl
    import math
    from atlasreader import get_statmap_info
    from src.help_fun import index_clusters_bin
    
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
        
    ids = ds.metadata.study_id.unique() # unique study ids (not contrasts)
    n = [ds.metadata[ds.metadata.study_id==iD].sample_sizes.tolist()[0] for iD in ids] # unique subjects
    print('Calculating ALE on {} experiments with {} contrasts, {} participants and {} foci.'.format(
          len(ids), len(ds.ids), np.sum(n), len(ds.coordinates)))
    print('Thresholding: voxel-level p < {} uncorrected, cluster-level p < 0.05 FWE-corrected, {} permutations.'.format(
          vox_thr, n_perm)) 
        
    # Estimate ALE map
    ds.update_path(work_dir)
    ale = ni.meta.cbma.ALE() # get ALE method
    ale_res = ale.fit(ds) # fit to data
    
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
    
    # Resulting cluster-level corrected "logp" volumes have the negative logarithm
    # to the base e of the cluster-level FWE-corrected p-value stored in all voxels 
    # of each cluster and require further thresholding at -ln(0.05) ~= 3.
    t = -math.log(0.05)        
        
    # Save & plot results
    if output == True:
        
        # Save NiMARE results
        ale_cres.save_maps(output_dir=work_dir, prefix=pref)     
            
        # Save clusters and cluster info if significant clusters were found
        if np.max(ale_cres.get_map('logp_level-cluster_corr-FWE_method-montecarlo').get_fdata()) > t: 
            
            # Threshold logp volume
            img_logp_thr = nl.image.threshold_img(
                ale_cres.get_map('logp_level-cluster_corr-FWE_method-montecarlo'), threshold=t)
            
            # Save binarized mask volume
            img_bin = nl.image.math_img('img>0', img=img_logp_thr)
            img_bin.to_filename(os.path.join(work_dir, pref + 'thresh_bin.nii.gz'))
            
            # Index clusters
            img_cl_idx, img_cl_4d, cl_sizes = index_clusters_bin(img_bin)
            if len(cl_sizes) > 1: # save volumes if more then one cluster 
                # Save volume with clusters indexed from 1 to n(clusters) and sorted from large to small
                img_cl_idx.to_filename(os.path.join(work_dir, pref + 'thresh_idx.nii.gz'))
                # Save 4D volume with one binarized 3D volume for each cluster
                img_cl_4d.to_filename(os.path.join(work_dir, pref + 'thresh_4d.nii.gz'))
            
            # Save thresholded ALE stats volume
            img_stat_thr = nl.image.math_img('img1*img2', img1=img_bin, img2=ale_cres.get_map('stat'))
            img_stat_thr.to_filename(os.path.join(work_dir, pref + 'thresh_stat.nii.gz'))
            
            # Save thresholded z score volume
            img_z_thr = nl.image.math_img('img1*img2', img1=img_bin, img2=ale_cres.get_map('z'))
            img_z_thr.to_filename(os.path.join(work_dir, pref + 'thresh_z.nii.gz'))
            
            # Extract cluster details from thresholded stat volume, default to using 
            # talairach atlases as does GingerALE
            if save_cl_info == True:
                cl_info, peak_info = get_statmap_info(
                    os.path.join(work_dir, pref + 'thresh_stat.nii.gz'), cluster_extent=0, 
                    voxel_thresh=0, direction='pos', atlas=atlases)
                cl_info.to_csv(os.path.join(work_dir, pref + 'thresh_stat_clusters.csv'), index=False)
                peak_info.to_csv(os.path.join(work_dir, pref + 'thresh_stat_peaks.csv'), index=False)
            
            # Plot
            if print_glassbrain == True:
                if glassbrain_title is None:
                    glassbrain_title2 = 'Thresholded ALE z-scores'
                else:
                    glassbrain_title2 = glassbrain_title+'; thresholded'
                nl.plotting.plot_glass_brain(img_bin, title=glassbrain_title2, 
                    draw_cross=False, cmap='RdBu_r', display_mode='lyrz')
        
            # Return if volumes saved
            print('ALE completed, significant clusters found. Results saved to: {}\n'.format(work_dir))
            return(ale_cres, img_stat_thr, img_bin)
        else:
            print('ALE completed, no significant clusters! NiMARE results saved to: {}\n'.format(work_dir))
            return(ale_cres, None, None)
    
    # Return if no volumes saved
    else:
        if np.max(ale_cres.get_map('logp_level-cluster_corr-FWE_method-montecarlo').get_fdata()) > t: 
            print('ALE completed, significant clusters found. No results saved.\n')
        else:
            print('ALE completed, no significant clusters! No results saved.\n')
        return(ale_cres, None, None)
    
    
#=============================================================================


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
    
    import os
    import numpy as np
    from nilearn.masking import intersect_masks
    from nilearn.plotting import plot_glass_brain
    
    if save_dir is None:
        save_dir = os.getcwd()
    if conj_save is None:
        conj_save = os.path.join(os.getcwd(), 'LOEO_conjunction.nii.gz')

    # get experiments and iterate over them
    exps = ds.ids
    ress = list()
    masks = list()
    for i, exp in enumerate(exps):
        print('\nCalculate ALE without experiment', i+1, exp)
        # remove current experiment from ds
        ds_x = ds.slice(np.delete(exps, i))
        # run ALE
        res, _, mask = run_ale(data=ds_x, work_dir=save_dir, pref=prefix+str(i+1)+'_', 
                                  vox_thr=v_thr, n_perm=n_perm, n_core=n_core, 
                                  print_glassbrain=False, save_cl_info=False)
        # save NiMARE results and cluster masks
        ress.append(res)
        masks.append(mask)

    # calculate and print conjunction of mask volumes
    conjunction = intersect_masks(masks, 1, connected=False) # calculate intersection
    conjunction.to_filename(conj_save) # save volume
    print('Conjunction volume saved to', conj_save)
    plot_glass_brain(conjunction, title='Conjunction of thresholded LOEO volumes', 
                     threshold=0, draw_cross=False, cmap='RdBu_r', display_mode='lyrz')
    
    # return
    return(ress, conjunction)

        
#=============================================================================


def exp_contributions(ale_ALL, ale_EXP, cl_masks, exp_ids):  
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
    
    from nilearn.input_data import NiftiMasker
    import numpy as np
    import pandas as pd
    
    # Initialize dataframe with experiment IDs
    contributions = pd.DataFrame({'ID': exp_ids})
         
    for i, cl in enumerate(cl_masks, start=1):
        # Set Nilearn masker
        masker = NiftiMasker(mask_img=cl)
        
        # Loop over experiments
        ale_means = list()
        for exp in ale_EXP:
            exp_ale_values = masker.fit_transform(exp) # get ale values within cluster
            exp_ale_mean = np.mean(exp_ale_values) # mean of within-cluster values
            ale_means.append(exp_ale_mean) # collect results
        
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
    
    print('Calculated experiment-wise contributions of {} experiments to {} clusters.'.format(
        len(contributions), len(cl_masks)))
    
    return(contributions)


#=============================================================================


def fsn(ds, cl_true, template=None, fsn_range=None, max_fsn_it=None, 
        quadrant=True, n_subj_null=None, n_coord_null=None, save_dir=None, 
        pref='fsn_', n_perm=1000, v_thr=0.001, seed=None, plot_gb=False):   
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
        save_dir, pref, n_perm, v_thr = ALE settings, see above
        seed = starting seed for sampling procedures
        plot_gb = if True, plots glass brain for each ALE.

    Output:
        list of NiMARE results of final ALEs for each cluster, table with FSN
        and NiMARE noise dataset with n = fsn_range[1] noise studies.
    """
    
    import os
    import numpy as np
    import pandas as pd
    import nimare as ni
    import nilearn as nl
    #import random as rd

    # Set variables ----------------------------------------------------------
    ds.coordinates['space'] = 'MNI' # has no influence but avoids warning when merging datasets
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
        n_subj_true = ds.get_metadata(field="sample_sizes")
        n_subj_null = rng.choice(n_subj_true, M)
    # get numbers of foci per study for noise studies from 'true' dataset
    if n_coord_null is None:
        n_coord_true = [len(ds.coordinates.x[ds.coordinates.id==iD]) for iD in ds.ids] 
        n_coord_null = rng.choice(n_coord_true, M)
      
    # Template --------------------------------------------------------------- 
    # If no template given, fetch MNI152 GM mask and resample to cl_true[0] space
    # produces warning: "pixdim[0] (qfac) should be 1 or -1; setting qfac to 1"
    # the result is okay, however.
    if template is None:
        template_1mm = nl.datasets.fetch_icbm152_brain_gm_mask(threshold=0.2)
        template = nl.image.resample_to_img(template_1mm, cl_true[0], interpolation='nearest')
        
    # Get possible MNI coordinates for noise studies = all voxels in template
    dat = template.get_fdata() # get matrix
    x, y, z = np.where(dat == 1.0) # get GM coordinates
    affine = template.affine # get affine matrix for MNI transform
    mni_all = nl.image.coord_transform(x=x, y=y, z=z, affine=affine) # mni coordinates
    mni_all = np.array(mni_all).T
    
    # FSN -------------------------------------------------------------------- 
    
    # define function to run ALE and return whether a significant cluster within
    # the location of the 'true' cluster was found
    def fsn_ale(ds_add, cl, prefix, title):
        # run ale
        res, _, mask = run_ale(data=ds_add, work_dir=save_dir, pref=prefix, 
                               save_cl_info=False, vox_thr=v_thr, n_perm=n_perm, 
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
        print('\nCluster {}:'.format(cl_i)) # print cluster number
        cl_pref = pref+'cl'+str(cl_i)+'_' # prefix to store ALE results
        
        # 1. Generate noise studies ------------------------------------------
        
        # Make sure that noise coordinates do not overlap with the target 
        # cluster by determining a "nogo" brain quadrant
        if quadrant == True:
            # Get center of mass of cluster 
            ctr = nl.plotting.find_parcellation_cut_coords(cl)[0]    
            # Get MNI coordinates we can sample from by removing coordinates 
            # belonging to the 'nogo' brain quadrant where 'ctr' is placed in
            q = [ (1 if c > 0 else -1) for c in ctr ]  
            mni_go = mni_all[(mni_all[:,0] * q[0] < 0) | 
                             (mni_all[:,1] * q[1] < 0) | 
                             (mni_all[:,2] * q[2] < 0), :]  
        else:
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
        print('\nFSN step 1: Add minimum {} noise studies.\n'.format(m))
        ds_add = ds.merge(ds_null.slice(ids_null[:m])) # add first m noise studies to ds
        fsn_res, check_cl = fsn_ale(ds_add, cl, cl_pref, 'FSN, {} noise studies'.format(m)) # run
        if check_cl == False:
            print('FSN finished: Cluster not significant, indicative of low robustness.')
            m_str = '< '+str(m) # true FSN is below minimum
            fsn_tab.append(['cl'+str(cl_i), m_str, m_str, m_str]) # save FSN values
        elif check_cl == True:
            print('Cluster significant, proceed.')
            
            # Step 2: run ale with maximum 'M' noise studies -----------------
            print('\nFSN step 2: Add maximum {} noise studies.\n'.format(M))
            ds_add = ds.merge(ds_null.slice(ids_null[:M])) # add M noise studies to ds
            fsn_res, check_cl = fsn_ale(ds_add, cl, cl_pref, 'FSN, {} noise studies'.format(M)) # run          
            if check_cl == True:
                print('FSN finished: Cluster still significant, indicative of highly influencial studies!')
                M_str = '> '+str(M) # true FSN is above maximum
                fsn_tab.append(['cl'+str(cl_i), M_str, M_str, M_str]) # save FSN values
            elif check_cl == False:
                print('Cluster not significant, proceed.')
                
                # Step 3 to n: run ale with N noise studies ------------------
                m_ = m
                M_ = M
                N = int((m_ + M_) / 2) # start with average of m and M, cut to the lower integer
                for ii in range(0, max_fsn_it):
                    print('\nFSN step 3+{}: Add {} noise studies.\n'.format(ii, N))
                    ds_add = ds.merge(ds_null.slice(ids_null[:N])) # add N noise studies to ds
                    fsn_res, check_cl = fsn_ale(ds_add, cl, cl_pref, 'FSN, {} noise studies'.format(N)) # run   
                    if check_cl == False:
                        print('Cluster not significant.')
                        m_ = m_     # minimum remains the same 
                        M_ = N      # last N is our new maximum
                        N = int((N + m_) / 2) # new N is average of mimum and old N
                    elif check_cl == True:
                        print('Cluster significant.')
                        m_ = N      # last N is our new minimum
                        M_ = M_     # maximum remains the same
                        N = int((N + M_) / 2) # new N is average of old N and maximum
                    print('New min FSN = {}, new max FSN = {}, new N = {}.'.format(m_, M_, N))
                    # loop out if FSN exactly determined, case 1: N==m_==M__;
                    # case 2: m_==n & M_==n+1, explanation: m_==n -> significant cluster, 
                    # M_==n+1 -> no significant cluster, thus FSN==m_
                    if (m_ == M_) or (M_ == m_+1):
                        print('Looping out after {} iterations, FSN exactly determined for the given set of noise studies.'.format(ii+1))
                        M_ = m_ 
                        break 
                        
                # Finished: Show and store results ---------------------------
                # If enough iterations where performed, FSN, min and max are the same.
                # If not enough iterations, the true FSN lies between min and max.
                print('\nCluster {}: FSN finished: FSN (~)= {} (min = {}, max = {})'.format(
                      cl_i, N, m_, M_))
                fsn_tab.append(['cl'+str(cl_i), N, m_, M_]) # save FSN values
     
        # store ALE results
        fsn_ress.append(fsn_res) # save last ALE result
    
    # write list to dataframe
    fsn_tab = pd.DataFrame(fsn_tab, columns=['cluster', 'appr_fsn', 'min_fsn', 'max_fsn'])
    # return
    return(fsn_tab, fsn_ress, ds_null)


#=============================================================================


def rel_abs_distributions(roi_vols, target_vol, target_descr):   
    """
    Quantifies overlap of a list of ROIs with, e.g., a given target network 
    parcellation. 
    Relative distribution = proportion of 'roi_vols' voxels within a given 
        'target_vol' cluster vs. all 'roi_vols' voxels.
    Absolute distribution  = proportion of 'roi_vols'-voxels within a given 
        'target_vol' cluster vs. all voxels within this 'target_vol' cluster.
    
    Input:
        roi_vols = list of roi volumes 
        roi_descr = list of roi names for output table
        target_vol = target volume with clusters/parcels indexed from 1 to n
        target_descr = list of target culster/parcel names for output table
    All volumes are resampled to space and voxelsize of the first roi.

    Output:
        pandas dataframe with relative and absolute distributions of each ROI 
        within each cluster of the target volume.     
    """
    
    from nilearn.image import resample_to_img, math_img
    import pandas as pd
   
    # resample target volume to space of first roi volume 
    t_vol = resample_to_img(target_vol, roi_vols)
    
    t_mat = t_vol.get_fdata().round() # get target data matrix
    target_nrs = range(1,len(target_descr)+1) # target cluster indices

    # Results table with target cluster names and total voxel numbers per cluster
    res = pd.DataFrame({'targetNr': target_nrs, 
                        'targetName': target_descr,
                        'nvox_target': [len(t_mat[t_mat == nr]) for nr in target_nrs]})  
    
    # Iterate over ROI list
    res_out = []
    for roi in roi_vols:
        r_vol = resample_to_img(roi, roi_vols[0]) # resample to first roi space
        r_bin = math_img('img > 0', img=r_vol) # binarize, just in case
        r_mat = r_bin.get_fdata() # get roi data matrix
        rt_mat = r_mat * t_mat # intersect roi with target matrix
        
        # RESULTS
        res_roi = res.copy()
        # voxels in ROI
        res_roi['nvox_roi'] = len(r_mat[r_mat == 1])
        # overlapping voxels per target cluster
        res_roi['nvox_overlap'] = [len(rt_mat[rt_mat == nr]) for nr in res_roi.targetNr]
        # distributions
        res_roi['relDistr'] = res_roi.nvox_overlap / res_roi.nvox_roi # relative
        res_roi['absDistr'] = res_roi.nvox_overlap / res_roi.nvox_target # absolute
        # save to list
        res_out.append(res_roi)

    return(res_out)
        
        
#=============================================================================


