# ========================================================================
# Figure S3 - ALE, MACM, sensitivity
# ========================================================================

#%% 
pp = '/Users/leonlotter/MAsync/project'
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

from os.path import join
import sys
sys.path.append(join(pp, 'src'))
from utils_image import get_size_of_rois
import pandas as pd
import numpy as np
from nilearn.image import math_img, load_img
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import to_hex
from nilearn.plotting import plot_glass_brain
from seaborn import regplot

# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values
c_ale001 = ListedColormap(cmap[0,:])


#%% Fig S3AB ALE, SCALE =========================================================================

ale_ps = join(wd, 'ale', 'aleNoPseudo_thresh_bin.nii.gz')
scale_z = join(wd, 'macm', 'scale_rTPJ_z.nii.gz')
scale_p = join(wd, 'macm', 'scale_rTPJ_logp.nii.gz')

# plot
figS3ab = plt.figure(figsize=(11,6), constrained_layout=False)
gs = figS3ab.add_gridspec(ncols=1, nrows=2, hspace=0)
# Row 1: ale no pseudo
ax1 = figS3ab.add_subplot(gs[0,0])
plot_ale_z = plot_glass_brain(ale_ps, figure=figS3ab, axes=ax1, display_mode='lyrz',
                            title='ALE (no pseudo-hyperscanning, p < .001)', cmap=c_ale001)
# Row 2: scale
ax2 = figS3ab.add_subplot(gs[1,0])
plot_ale01 = plot_glass_brain(math_img('z*(p>3)', z=scale_z, p=scale_p), 
                              figure=figS3ab, axes=ax2, display_mode='lyrz',
                              title='SCALE (p < .05)', cmap=c_ale001)

plt.savefig(join(sdir, 'figS3ab.pdf'), transparent=True, bbox_inches='tight')


# %% S3C macm-ale cor

# load correlation, atlas and macm result
cor_macm_ale = pd.read_csv(join(wd, 'macm', 'cor_macm_ale.csv'))
cor_macm_ale.rename(columns={'dat1': 'ale', 'dat2': 'macm'}, inplace=True)
atlas = load_img(join(wd, 'atlases', 'Schaefer100-7_TianS1_2mm.nii.gz'))
macm_idx = load_img(join(wd, 'macm', 'macm_rTPJ_thresh_idx.nii.gz'))
macm_idx = math_img('macm * (macm<10)', macm=macm_idx)
macm_labels = pd.read_csv(join(wd, 'macm', 'macm_rTPJ_thresh_labels.csv'), header=None)
macm_labels = macm_labels[0].to_list()[0:10]

# MARKER COLORS FOR SCATTER PLOT 
cl_info = get_size_of_rois(macm_idx) # get macm cluster indices and sizes
colors = [to_hex(cmap[i,:]) for i in range(len(cmap))]
cor_macm_ale['color'] = 'lightgrey' # non-cluster marker colors
cor_macm_ale['cluster'] = ''
# iterate over macm clusters
for i_cl in cl_info['idx']:
    cl_dat = macm_idx.get_fdata() # get cluster data
    cl_parcel = (cl_dat==i_cl) * atlas.get_fdata() # intersect current cluster with atlas
    cl_size = cl_info['size'][cl_info['idx'] == i_cl].values[0] # size of current cluster 
    cl_parcel_sizes = get_size_of_rois(cl_parcel) # sizes of parcel-cluster intersections
    # keep parcel-cluster intersections that make at least 10% of cluster size
    cl_parcel_thr = cl_parcel_sizes[cl_parcel_sizes['size'] >= cl_size * 0.10] 
    # set color and save cluster name for all parcels overlapping with the current cluster
    cor_macm_ale['color'][cor_macm_ale['idx'].isin(cl_parcel_thr['idx'])] = colors[int(i_cl)-1] 
    cor_macm_ale['cluster'][cor_macm_ale['idx'].isin(cl_parcel_thr['idx'])] = macm_labels[int(i_cl)-1]

# get data and plot
figS3c = plt.figure(figsize=(5,4), constrained_layout=False)
corplot = regplot(data=cor_macm_ale, x='ale', y='macm', color='black', 
                  scatter_kws={'facecolors':cor_macm_ale.color}) 
# legend
for i, grp in cor_macm_ale.groupby(['color', 'cluster']):
    grp.plot(kind='scatter', x='ale', y='macm', c=i[0], ax=corplot, label=i[1], zorder=0)       
corplot.legend(bbox_to_anchor=(1.02, 0.88))
# labels
corplot.set_xlabel('INS ALE Z-scores', fontsize=12)
corplot.set_ylabel('INS MACM Z-scores', fontsize=12)

plt.savefig(join(sdir, 'figS3c.pdf'), transparent=True, bbox_inches='tight')

## %%

# %%
