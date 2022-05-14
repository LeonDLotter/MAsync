# ========================================================================
# Figure S5 - additional fNIRS
# ========================================================================

#%% 
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from nilearn.image import math_img
from nilearn.datasets import fetch_surf_fsaverage
from neuromaps.transforms import mni152_to_fsaverage
from surfplot import Plot

# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values
c_ale001 = ListedColormap(cmap[0,:])

# atlas
atlas = join(wd, 'atlases', 'Schaefer100-7_2mm.nii.gz')

# results tables
fnirs_res = pd.read_csv(join(wd, 'fnirs', 'fnirs_atlas_result.csv'))
fnirs_res_sel = pd.read_csv(join(wd, 'fnirs', 'fnirs_atlas_result_sel.csv'))
fnirs_res_perm = pd.read_csv(join(wd, 'fnirs', 'fnirs_atlas_perm_ch_ratio_sub_all.csv'), index_col=0)

# Surface plot function
def plot_surf(fnirs_vol, outline_vol, title, save_path):
    # get template
    fsaverage = fetch_surf_fsaverage()
    # make plot
    p = Plot(surf_lh=fsaverage['pial_left'], surf_rh=fsaverage['pial_right'], layout='row', 
            size=(1000, 250), zoom=1.25, views=['lateral', 'medial'], flip=False)
    # add fnirs map
    map_lh, map_rh = mni152_to_fsaverage(fnirs_vol, method='nearest', fsavg_density='10k')
    p.add_layer({'left': map_lh, 'right': map_rh}, cmap='viridis', cbar=True)
    # outline significant parcels
    outline_lh, outline_rh = mni152_to_fsaverage(outline_vol, method='nearest', fsavg_density='10k')
    p.add_layer({'left': outline_lh, 'right': outline_rh}, cmap=c_ale001, as_outline=True, cbar=False)
    # legend
    kws = dict(location='bottom', draw_border=True, aspect=15, shrink=.2, decimals=0, pad=-0.08)
    # plot
    fig = p.build(cbar_kws=kws) 
    plt.title(title, loc="left", y=0.9, size=15, c='w', bbox=dict(facecolor='k', edgecolor='k'))
    # save
    plt.savefig(save_path, bbox_inches='tight', transparent='true')


#%% figS5a - histograms for main result

fnirs_res.sort_values('ch_ratio_sub_all', ascending=False, inplace=True)
fnirs_res.reset_index(drop=True, inplace=True)

# make plot
fig, axes = plt.subplots(1,10, figsize=(21,2.2))
for i, region_idx in enumerate(fnirs_res.region_idx.reset_index(drop=True)[:10]):
    # plot histogram
    sns.histplot(fnirs_res_perm.loc[region_idx,:], kde=True, ax=axes[i])
    # index & p value
    ind = fnirs_res[fnirs_res.region_idx==region_idx]['ch_ratio_sub_all'].values
    p = fnirs_res[fnirs_res.region_idx==region_idx]['p_ch_ratio_sub_all'].values
    c = 'red' if p < .05 else 'blue'
    # plot line
    axes[i].axvline(x=ind, c=c)
    # plot p value
    axes[i].annotate(f'{p[0]:.3f}', xy=(ind+20,250), ha='left', c=c)
    # plot title
    axes[i].set(xlabel=None, ylabel=None)
    axes[i].set_title(fnirs_res.region[i], size=10)
plt.tight_layout()
plt.savefig(join(sdir, 'figS5a.pdf'), bbox_inches='tight', transparent='true')


#%% figS5b - all indices, all experiments

# number of channels
fnirs = join(wd, 'fnirs', 'fnirs_atlas_syncChannels.nii.gz')
test = ' | '.join([f'(a=={i})' for i in fnirs_res.query("p_n_ch_sync < 0.05").region_idx.to_list()])
fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
plot_surf(fnirs, fnirs_sig, 'All: Number of channels', join(sdir, 'figS5b1.pdf'))

# ratio, weighted by subjects
fnirs = join(wd, 'fnirs', 'fnirs_atlas_ratioSyncAllChannelsSubj.nii.gz')
test = ' | '.join([f'(a=={i})' for i in fnirs_res.query("p_ch_ratio_sub_all < 0.05").region_idx.to_list()])
fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
plot_surf(fnirs, fnirs_sig, 'All: Ratio of INS to all channels weighted by number of subjects', join(sdir, 'figS5b2.pdf'))

# ratio weighted by experiments
fnirs = join(wd, 'fnirs', 'fnirs_atlas_ratioSyncAllChannelsExps.nii.gz')
test = ' | '.join([f'(a=={i})' for i in fnirs_res.query("p_ch_ratio_exp_all < 0.05").region_idx.to_list()])
fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
plot_surf(fnirs, fnirs_sig, 'All: Ratio of INS to all channels weighted by number of experiments', join(sdir, 'figS5b3.pdf'))


#%% figS5c - all indices, selected experiments

# number of channels
fnirs = join(wd, 'fnirs', 'fnirs_atlas_syncChannels_sel.nii.gz')
test = ' | '.join([f'(a=={i})' for i in fnirs_res_sel.query("p_n_ch_sync < 0.05").region_idx.to_list()])
fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
plot_surf(fnirs, fnirs_sig, 'Restricted: Number of channels', join(sdir, 'figS5c1.pdf'))

# ratio, weighted by subjects
fnirs = join(wd, 'fnirs', 'fnirs_atlas_ratioSyncAllChannelsSubj_sel.nii.gz')
test = ' | '.join([f'(a=={i})' for i in fnirs_res_sel.query("p_ch_ratio_sub_all < 0.05").region_idx.to_list()])
fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
plot_surf(fnirs, fnirs_sig, 'Restricted: Ratio of INS to all channels weighted by number of subjects', join(sdir, 'figS5c2.pdf'))

# ratio weighted by experiments
fnirs = join(wd, 'fnirs', 'fnirs_atlas_ratioSyncAllChannelsExps_sel.nii.gz')
test = ' | '.join([f'(a=={i})' for i in fnirs_res_sel.query("p_ch_ratio_exp_all < 0.05").region_idx.to_list()])
fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
plot_surf(fnirs, fnirs_sig, 'Restricted: Ratio of INS to all channels weighted by number of experiments', join(sdir, 'figS5c3.pdf'))



# %%
