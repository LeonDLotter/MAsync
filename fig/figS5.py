# ========================================================================
# Figure S5 - additional fNIRS
# ========================================================================

#%% 
wd = '/Users/llotter/MAsync/project/data'
sdir = '/Users/llotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from nilearn.image import load_img, math_img, new_img_like
from nilearn.datasets import fetch_surf_fsaverage
from neuromaps.transforms import mni152_to_fsaverage
from surfplot import Plot

# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values
c_ale001 = ListedColormap(cmap[0,:])

# atlas
atlas = load_img(join(wd, 'atlases', 'Schaefer100-7_2mm.nii.gz'))

# results tables
fnirs_res = pd.read_csv(join(wd, 'fnirs', 'fnirs_atlas_result_perm.csv'))
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
    plt.savefig(save_path, bbox_inches='tight', transparent='true', dpi=200)


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


#%% figS5b & c

for dataset, dataset_title, abc in zip(["", "_sel"],
                                       ["All", "Selected"],
                                       ["b", "c"]):
    for j, file, index, title in zip([1, 2, 3],
                                     ["fnirs_atlas_syncChannels", 
                                      "fnirs_atlas_ratioSyncAllChannelsSubj", 
                                      "fnirs_atlas_ratioSyncAllChannelsExps"],
                                     ["p_n_ch_sync", 
                                      "p_ch_ratio_sub_all", 
                                      "p_ch_ratio_exp_all"],
                                     ['Number of channels', 
                                      'Ratio of INS to all channels weighted by subject number', 
                                      'Ratio of INS to all channels weighted by experiment number']):
        # plot
        fnirs = join(wd, 'fnirs', f'{file}{dataset}.nii.gz')
        test = ' | '.join([f'(a=={i})' for i in fnirs_res.query(f"{index} < 0.05").region_idx.to_list()])
        fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
        plot_surf(fnirs, fnirs_sig, f'{dataset_title}: {title}', join(sdir, f'figS5{abc}{j}.pdf'))



# %% figS5d & e

for var, var_title, abc in zip(["rand_median", "rand_perc"],
                               ["Randomized (M)", "Randomized (%)"],
                               ["d", "e"]):
    for j, file, index, title in zip([1, 2, 3],
                                     ["fnirs_atlas_syncChannels", 
                                      "fnirs_atlas_ratioSyncAllChannelsSubj", 
                                      "fnirs_atlas_ratioSyncAllChannelsExps"],
                                     ["p_n_ch_sync", 
                                      "p_ch_ratio_sub_all", 
                                      "p_ch_ratio_exp_all"],
                                     ['Number of channels', 
                                      'Ratio of INS to all channels weighted by subject number', 
                                      'Ratio of INS to all channels weighted by experiment number']):
        # get p volume
        p_vol = new_img_like(atlas, data=np.zeros(atlas.shape))
        for region_idx in fnirs_res.region_idx.to_list():
            p = fnirs_res.query("region_idx==@region_idx")[f"{index}_{var}"].values
            if var=="rand_median":
                p = -np.log10(p)
            else:
                p = p * 100
            p_vol = math_img(f'np.where(a=={region_idx}, {p}, p_vol)', a=atlas, p_vol=p_vol)
        # get p <  0.05 indicator
        fnirs = join(wd, 'fnirs', f'{file}.nii.gz')
        test = ' | '.join([f'(a=={i})' for i in fnirs_res.query(f"{index} < 0.05").region_idx.to_list()])
        fnirs_sig = math_img(f'np.where({test}, 1, 0)', a=atlas)
        # plot
        plot_surf(p_vol, fnirs_sig, f'{var_title}: {title}', join(sdir, f'figS5{abc}{j}.pdf'))

# %%
