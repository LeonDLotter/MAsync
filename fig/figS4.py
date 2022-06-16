# ========================================================================
# Figure S4 - prior meta-analyses
# ========================================================================

#%% 
wd = '/Users/llotter/MAsync/project/data'
sdir = '/Users/llotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
from nilearn.image import math_img
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import cm
from nilearn.plotting import plot_glass_brain
import seaborn as sns

# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values

c = ListedColormap(cmap[4,:])
c_nw = ListedColormap(cmap[1,:])
c_cl = ListedColormap(cmap[0,:])

# function to plot bars
def plot_bars(data, ax, x, roiName=['rTPJ', 'MACM'], xlim=(0,1),
              palette={'rTPJ':cmap[0,:], 'MACM':cmap[1,:]}):
    data.rename(columns=dict(relDistr='Relative distribution', absDistr='Absolute distribution'), inplace=True)
    data.roiName = roiName
    x = 'Relative distribution' if x=='rel' else 'Absolute distribution'
    sns.barplot(y=data.roiName, x=data[x], orient='h', ax=ax, palette=palette)
    ax.set(ylabel='', xlim=xlim)
    ax.set_xlabel(x, fontweight='bold', size=11)
    ax.set_yticklabels(roiName, size=11)
    return(ax)


# %%

feng = join(wd, 'atlases', 'Feng_socialNonsocial.nii.gz')
feng_cl = join(wd, 'context', 'feng_ale_conjunction.nii.gz')
feng_nw = join(wd, 'context', 'feng_macm_conjunction.nii.gz')
feng_dist = pd.read_csv(join(wd, 'context', 'feng_distributions.csv'))

schurz_a = join(wd, 'atlases', 'Schurz_ToM_affective.nii.gz')
schurz_a_cl = join(wd, 'context', 'schurz_affective_ale_conjunction.nii.gz')
schurz_a_nw = join(wd, 'context', 'schurz_affective_macm_conjunction.nii.gz')
schurz_a_dist = pd.read_csv(join(wd, 'context', 'schurz_affective_distributions.csv'))

schurz_i = join(wd, 'atlases', 'Schurz_ToM_intermediate.nii.gz')
schurz_i_cl = join(wd, 'context', 'schurz_intermediate_ale_conjunction.nii.gz')
schurz_i_nw = join(wd, 'context', 'schurz_intermediate_macm_conjunction.nii.gz')
schurz_i_dist = pd.read_csv(join(wd, 'context', 'schurz_intermediate_distributions.csv'))

schurz_c = join(wd, 'atlases', 'Schurz_ToM_cognitive.nii.gz')
schurz_c_cl = join(wd, 'context', 'schurz_cognitive_ale_conjunction.nii.gz')
schurz_c_nw = join(wd, 'context', 'schurz_cognitive_macm_conjunction.nii.gz')
schurz_c_dist = pd.read_csv(join(wd, 'context', 'schurz_cognitive_distributions.csv'))

ficco = join(wd, 'atlases', 'Ficco_PredEncodeViolate.nii.gz')
ficco_cl = join(wd, 'context', 'ficco_ale_conjunction.nii.gz')
ficco_nw = join(wd, 'context', 'ficco_macm_conjunction.nii.gz')
ficco_dist = pd.read_csv(join(wd, 'context', 'ficco_distributions.csv'))

bzdok = join(wd, 'atlases', 'BzdokTPJ_combined_2mm.nii.gz')
bzdok_cl = math_img('b*i', b=bzdok, i=join(wd, 'ale', 'ale_thresh_bin.nii.gz')) 
bzdok_dist = pd.read_csv(join(wd, 'context', 'bzdok_distributions.csv'))


# plot
figS4 = plt.figure(figsize=(15,15.5), constrained_layout=False)
h = 0.4
gs = figS4.add_gridspec(ncols=3, nrows=18, hspace=0, width_ratios=[1,0.3,0.3],
                           height_ratios=[h,1,h,h,1,h,h,1,h,h,1,h,h,1,h,h,1,h])

# Row 1: feng
ax1 = figS4.add_subplot(gs[0:3,0])
pl_feng = plot_glass_brain(feng, figure=figS4, axes=ax1, display_mode='lyrz',
                           title='Feng et al.: social interaction > no interaction', cmap=c)
pl_feng.add_overlay(feng_nw, cmap=c_nw)
pl_feng.add_overlay(feng_cl, cmap=c_cl)
ax11 = figS4.add_subplot(gs[1,1])
ax11 = plot_bars(feng_dist, ax11, x='rel')
ax12 = figS4.add_subplot(gs[1,2])
ax12 = plot_bars(feng_dist, ax12, x='abs')

# Row 2: schurz affective
ax2 = figS4.add_subplot(gs[3:6,0])
pl_schurz_a = plot_glass_brain(schurz_a, figure=figS4, axes=ax2, display_mode='lyrz',
                              title='Schurz et al.: affective ToM', cmap=c)
pl_schurz_a.add_overlay(schurz_a_nw, cmap=c_nw)
pl_schurz_a.add_overlay(schurz_a_cl, cmap=c_cl)
ax21 = figS4.add_subplot(gs[4,1])
plot_bars(schurz_a_dist, ax21, 'rel')
ax22 = figS4.add_subplot(gs[4,2])
plot_bars(schurz_a_dist, ax22, 'abs')

# Row 3: schurz intermediate
ax3 = figS4.add_subplot(gs[6:9,0])
pl_schurz_i = plot_glass_brain(schurz_i, figure=figS4, axes=ax3, display_mode='lyrz',
                              title='Schurz et al.: intermediate ToM', cmap=c)
pl_schurz_i.add_overlay(schurz_i_nw, cmap=c_nw)
pl_schurz_i.add_overlay(schurz_i_cl, cmap=c_cl)
ax31 = figS4.add_subplot(gs[7,1])
plot_bars(schurz_i_dist, ax31, 'rel')
ax32 = figS4.add_subplot(gs[7,2])
plot_bars(schurz_i_dist, ax32, 'abs')

# Row 4: schurz cognitive
ax4 = figS4.add_subplot(gs[9:12,0])
pl_schurz_c = plot_glass_brain(schurz_c, figure=figS4, axes=ax4, display_mode='lyrz',
                              title='Schurz et al.: cognitive ToM', cmap=c)
pl_schurz_c.add_overlay(schurz_c_nw, cmap=c_nw)                                       
pl_schurz_c.add_overlay(schurz_c_cl, cmap=c_cl) 
ax41 = figS4.add_subplot(gs[10,1])
plot_bars(schurz_c_dist, ax41, 'rel')
ax42 = figS4.add_subplot(gs[10,2])
plot_bars(schurz_c_dist, ax42, 'abs')

# Row 5: ficco 
ax5 = figS4.add_subplot(gs[12:15,0])
pl_ficco = plot_glass_brain(ficco, figure=figS4, axes=ax5, display_mode='lyrz',
                           title='Ficco et al.: prediction encoding/violation', cmap=c)
pl_ficco.add_overlay(ficco_nw, cmap=c_nw)
pl_ficco.add_overlay(ficco_cl, cmap=c_cl)
ax11 = figS4.add_subplot(gs[13,1])
ax11 = plot_bars(ficco_dist, ax11, x='rel')
ax12 = figS4.add_subplot(gs[13,2])
ax12 = plot_bars(ficco_dist, ax12, x='abs')

# Row 6: bzdok rTPJ
ax6 = figS4.add_subplot(gs[15:,0])
pl_bzdok = plot_glass_brain(bzdok, figure=figS4, axes=ax6, display_mode='lyrz',
                              title='Bzdok et al.: rTPJ subunits', cmap=c)
pl_bzdok.add_overlay(bzdok_cl, cmap=c_cl) 
ax51 = figS4.add_subplot(gs[16,1])
plot_bars(bzdok_dist, ax51, 'rel', roiName=['aTPJ', 'pTPJ'],
          palette={'aTPJ':cmap[0,:], 'pTPJ':cmap[0,:]}, xlim=None)
ax52 = figS4.add_subplot(gs[16,2])
plot_bars(bzdok_dist, ax52, 'abs', roiName=['aTPJ', 'pTPJ'],
          palette={'aTPJ':cmap[0,:], 'pTPJ':cmap[0,:]}, xlim=None)

plt.savefig(join(sdir, 'figS4.pdf'), transparent=True, bbox_inches='tight')



# %%
