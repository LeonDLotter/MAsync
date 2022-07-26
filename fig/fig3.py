# ========================================================================
# Figure 3 - Brains & Functional Context
# ========================================================================

#%% 
wd = '/Users/llotter/MAsync/project/data'
sdir = '/Users/llotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
from nilearn.image import math_img
from nilearn.plotting import plot_glass_brain, plot_connectome
from nilearn.datasets import fetch_surf_fsaverage
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from seaborn import heatmap
from neuromaps.transforms import mni152_to_fsaverage
from surfplot import Plot
from utils_radar import *
import seaborn as sns
from adjustText import adjust_text
from nimare.transforms import p_to_z


# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values
c_continous = None
c_ale001 = ListedColormap(cmap[0,:])
c_ale01 = ListedColormap(cmap[[0,0,0,0,0,0,0,0,3],:])
c_alefnirs = ListedColormap(cmap)
c_overlay = ListedColormap('#000000')
c_macm =  ListedColormap(np.concatenate((cmap, cmap)))


#%% Fig 2AC ALE =========================================================================

ale_z = join(wd, 'ale', 'ale_z.nii.gz')
ale_cl001 = join(wd, 'ale', 'ale_thresh_idx.nii.gz')
ale_cl01 = math_img('np.where(i>3, 0, i)', i=join(wd, 'ale', 'ale01_thresh_idx.nii.gz'))
loeo_cl001 = join(wd, 'loeo', 'loeo_conjunction.nii.gz')
loeo_cl01 = join(wd, 'loeo', 'loeo01_conjunction.nii.gz')
macm_idx = join(wd, 'macm', 'macm_rTPJ_thresh_idx.nii.gz')
macm_bin = join(wd, 'macm', 'macm_rTPJ_thresh_bin.nii.gz')
macm_labels = pd.read_csv(join(wd, 'macm', 'macm_rTPJ_thresh_labels.csv'), header=None)
macm_labels = macm_labels[0].to_list()[0:10]
macm_coord = pd.read_csv(join(wd, 'macm', 'macm_rTPJ_thresh_stat_clusters.csv'))
macm_coord = macm_coord.loc[:9, ['peak_x', 'peak_y', 'peak_z']]

# plot
fig2ac = plt.figure(figsize=(11,12.5), constrained_layout=False)
gs = fig2ac.add_gridspec(ncols=1, nrows=4, hspace=0)
# Row 1: ale z
ax1 = fig2ac.add_subplot(gs[0,0])
plot_ale_z = plot_glass_brain(ale_z, figure=fig2ac, axes=ax1, display_mode='lyrz',
                            title='ALE (no threshold)', cmap=c_continous)
# Row 2: ale cl, p < 0.001
ax2 = fig2ac.add_subplot(gs[1,0])
plot_ale01 = plot_glass_brain(ale_cl001, figure=fig2ac, axes=ax2, display_mode='lyrz',
                              title='ALE (p < .001)', cmap=c_ale001)
plot_ale01.add_overlay(loeo_cl001, cmap=c_overlay)
# Row 3: ale cl, p < 0.01
ax3 = fig2ac.add_subplot(gs[2,0])
plot_ale02 = plot_glass_brain(ale_cl01, figure=fig2ac, axes=ax3, display_mode='lyrz',
                              title='ALE (p < .01)', cmap=c_ale01)
plot_ale02.add_overlay(loeo_cl01, cmap=c_overlay)
# Row 4: macm
ax1 = fig2ac.add_subplot(gs[3,0])
plot_macm = plot_glass_brain(macm_idx, figure=fig2ac, axes=ax1, display_mode='lyrz',
                             title='MACM', cmap=c_macm, vmin=0)
plot_macm.add_contours(math_img('i*(i<2)',i=ale_cl001), colors='k', linewidths=0.7)                            

plt.savefig(join(sdir, 'fig3ac.pdf'), transparent=True, bbox_inches='tight')


#%% Fig 2B1 fnirs surf =======================================================================

fnirs = join(wd, 'fnirs', 'fnirs_atlas_ratioSyncAllChannelsSubj.nii.gz')
fnirs_sig = math_img('np.where((a==43) | (a==74) | (a==53) | (a==21), 1, 0)', 
                     a=join(wd, 'atlases', 'Schaefer100-7_2mm.nii.gz'))

# template
fsaverage = fetch_surf_fsaverage()
map_lh, map_rh = mni152_to_fsaverage(fnirs, method='nearest', fsavg_density='10k')

# make plot
p = Plot(surf_lh=fsaverage['pial_left'], surf_rh=fsaverage['pial_right'], layout='row', 
         size=(1000, 250), zoom=1.25, views=['lateral', 'medial'], flip=False)

# add fnirs map
p.add_layer({'left': map_lh, 'right': map_rh}, cmap='viridis', cbar=True)

# outline significant parcels
outline_lh, outline_rh = mni152_to_fsaverage(fnirs_sig, method='nearest', fsavg_density='10k')
p.add_layer({'left': outline_lh, 'right': outline_rh}, cmap=c_ale001, as_outline=True, cbar=False)

kws = dict(location='bottom', draw_border=True, aspect=15, shrink=.2, decimals=0, pad=-0.08)
fig2b1 = p.build(cbar_kws=kws) 
plt.title('fNIRS (Ratio of INS to all channels weighted by number of subjects)', 
          loc="left", y=0.9, size=15, c='w', bbox=dict(facecolor='k', edgecolor='k'))

plt.savefig(join(sdir, 'fig3b1.pdf'), dpi=400, bbox_inches='tight', transparent='true')


#%% Fig 2B2 fnirs gb =========================================================================

fnirs_ale_cl = join(wd, 'ale', 'aleFNIRS_thresh_idx.nii.gz')
c_alefnirs = get_cmap("viridis")
c_alefnirs = ListedColormap(np.concatenate([
    cmap[[0,0,0,0,0,0],:],
    np.array([[0.477504, 0.821444, 0.318195, 1.0],
              [0.134692, 0.658636, 0.517649, 1.0],
              [0.266941, 0.748751, 0.440573, 1.0]])]))

# plot
fig2b2, ax = plt.subplots(1, figsize=(11,3), constrained_layout=False)
# ale fnirs
plot_fnirs = plot_glass_brain(fnirs_ale_cl, figure=fig2b2, axes=ax, display_mode='lyrz',
                              title='ALE (fMRI & fNIRS, p < .001)', cmap=c_alefnirs)
plot_fnirs.add_contours(fnirs_sig, colors='k', linewidths=0.5)  

plt.savefig(join(sdir, 'fig3b2.pdf'), bbox_inches='tight', transparent='true')


# %%% =========================================================================
# Fig 2D1 rsfc

rsfc_cor = pd.read_csv(join(wd, 'rsfc','rsfc_cormat.csv'))
rsfc_cor.columns, rsfc_cor.index = macm_labels, macm_labels

fig2c = plt.figure(figsize=(5.7,6), constrained_layout=False)
gs = fig2c.add_gridspec(ncols=2, nrows=2, hspace=0, wspace=-0.03,
                        width_ratios=[1.2,1], height_ratios=[1,1.2])
ax1 = fig2c.add_subplot(gs[0,0])
plot_rsfc1 = plot_connectome(np.array(rsfc_cor), np.array(macm_coord), figure=fig2c, 
                            title='RSFC', axes=ax1, display_mode='l', node_color=cmap)
ax2 = fig2c.add_subplot(gs[0,1])
plot_rsfc2 = plot_connectome(np.array(rsfc_cor), np.array(macm_coord), figure=fig2c, 
                            axes=ax2, display_mode='y', node_color=cmap)
ax3 = fig2c.add_subplot(gs[1,0])
plot_rsfc1 = plot_connectome(np.array(rsfc_cor), np.array(macm_coord), figure=fig2c, 
                            axes=ax3, display_mode='r', node_color=cmap)
ax4 = fig2c.add_subplot(gs[1,1])
plot_rsfc2 = plot_connectome(np.array(rsfc_cor), np.array(macm_coord), figure=fig2c, 
                            axes=ax4, display_mode='z', node_color=cmap)

plt.savefig(join(sdir, 'fig3d1.pdf'), transparent=True, bbox_inches='tight')


#%% Fig 2D2 heatmap =========================================================================

# unthresholded correlation data
rsfc_cor_nothresh = pd.read_csv(join(wd, 'rsfc','rsfc_cormat_nothresh.csv'))
rsfc_cor_nothresh.columns, rsfc_cor_nothresh.index = macm_labels, macm_labels

# annotation array
rsfc_cor_annot = np.where(rsfc_cor>0, '✱', '')

fig2c = heatmap(rsfc_cor_nothresh, cmap='bwr', center=0,
                square=True, annot=rsfc_cor_annot, fmt="",   
                annot_kws={'fontsize':15, 'color':'white'})
for _, spine in fig2c.spines.items():
    spine.set_visible(True)

fig2c = fig2c.get_figure()

plt.savefig(join(sdir, 'fig3d2.pdf'), transparent=True, bbox_inches='tight') 


# %% Fig 2E RSN =====================================================================

#sns.set_style("white")

rsn = pd.read_csv(join(wd, 'context', 'rsn_distributions.csv'))
colors = [cmap[0,:], cmap[1,:], cmap[0,:], cmap[1,:]]
labs = ['rTPJ cluster', 'MACM network']
ticks = [np.arange(0),
         np.arange(0),
         np.arange(0),
         np.arange(0)]

# prepare data
data = pd.concat([rsn[rsn.roiName=='ale_rTPJ_cluster'].relDistr.reset_index(drop=True),
                  rsn[rsn.roiName=='macm_rTPJ_network'].relDistr.reset_index(drop=True),
                  rsn[rsn.roiName=='ale_rTPJ_cluster'].absDistr.reset_index(drop=True),
                  rsn[rsn.roiName=='macm_rTPJ_network'].absDistr.reset_index(drop=True)], 
                  axis=1)
data.columns = ["ale_rel", "macm_rel", "ale_abs", "macm_abs"]
data_p = pd.concat([rsn[rsn.roiName=='ale_rTPJ_cluster'].relDistr_p.reset_index(drop=True),
                    rsn[rsn.roiName=='macm_rTPJ_network'].relDistr_p.reset_index(drop=True),
                    rsn[rsn.roiName=='ale_rTPJ_cluster'].absDistr_p.reset_index(drop=True),
                    rsn[rsn.roiName=='macm_rTPJ_network'].absDistr_p.reset_index(drop=True)], 
                  axis=1)
data_p.columns = data.columns
data_q = pd.concat([rsn[rsn.roiName=='ale_rTPJ_cluster'].relDistr_q.reset_index(drop=True),
                    rsn[rsn.roiName=='macm_rTPJ_network'].relDistr_q.reset_index(drop=True),
                    rsn[rsn.roiName=='ale_rTPJ_cluster'].absDistr_q.reset_index(drop=True),
                    rsn[rsn.roiName=='macm_rTPJ_network'].absDistr_q.reset_index(drop=True)], 
                  axis=1)
data_q.columns = data.columns

# label rotation
N = 7
theta = np.array([n / float(N) * 2 * np.pi for n in range(N)])

## PLOT
fig_rsn = plt.figure(figsize=(10, 4.8))
#fig_rsn = plt.figure(figsize=(6,7))

ax1 = plt.subplot(141, polar=True, zorder=0)
ax2 = plt.subplot(142, polar=True, zorder=0)
ax3 = plt.subplot(143, polar=True, zorder=0)
ax4 = plt.subplot(144, polar=True, zorder=0)
#ax1 = plt.subplot(221, polar=True, zorder=0)
#ax2 = plt.subplot(222, polar=True, zorder=0)
#ax3 = plt.subplot(223, polar=True, zorder=0)
#ax4 = plt.subplot(224, polar=True, zorder=0)
line = []
for i, (ax, analysis) in enumerate(zip([ax1,ax2,ax3,ax4], data.columns)):
    # values
    line.append(ax.plot(theta, data[analysis], linewidth=1.5, color=colors[i], alpha=0.75, zorder=2))
    #ax.scatter(theta[data_q[analysis] < 0.05], data[analysis][data_q[analysis] < 0.05],
    #           marker="$☆$", color="k", s=150, alpha=0.5, zorder=10)
    close_line(line[i][0])
    ax.fill(theta, data[analysis], facecolor=colors[i], alpha=0.3, zorder=2)
    # y axis
    #ax.set_rticks(ticks[i])
    ax.set_rlabel_position(80)
    # adjust size
    labels = rsn.targetName[0:7].copy()
    for ii, (p,q) in enumerate(zip(data_p[analysis],data_q[analysis])):
        if p < 0.05:
            labels[ii] = "$\\bf{" + labels[ii] + "}$"
        if q < 0.05:
            labels[ii] = labels[ii] + "*"
    set_rotated_labels(labels, fig_rsn, ax, 0.18)
    
# legend
leg_loc = 'lower right'
leg_pos = (1.6,-0.65)
#leg_loc = 'center'
#leg_pos = (1.2,-0.27)
leg_abs = ax1.legend([l[0] for l in line[0:2]], labs, title="$\\bf{Relative}$", 
                    loc=leg_loc, bbox_to_anchor=leg_pos)
leg_rel = ax3.legend([l[0] for l in line[2:4]], labs, title="$\\bf{Absolute}$", 
                    loc=leg_loc, bbox_to_anchor=leg_pos)
# increase legend linewidth
for l1,l2 in leg_abs.get_lines(),leg_rel.get_lines():
    l1.set_linewidth(3)  
    l2.set_linewidth(3)

# adjust layout
plt.subplots_adjust(wspace = 0.25)
#plt.subplots_adjust(wspace = 0.4, hspace = 0.4)

# save and close figure
plt.savefig('fig3e.pdf', transparent=True, bbox_inches='tight')


# %% Fig 2F neurosynth =====================================================================

ns_rTPJ = pd.read_csv(join(wd, 'context', 'topics_roi_ale_rTPJ.csv'))
ns_rTPJ.Term = [t.split('__')[1] for t in ns_rTPJ.Term]
ns_rTPJ.sort_values(by='Term', inplace=True)
ns_rTPJ.reset_index(drop=True, inplace=True)

ns_whole = pd.read_csv(join(wd, 'context', 'topics_wholebrain.csv'), index_col=0)
ns_whole['topic'] = ['_'.join(t.split('_')[2:6]) for t in ns_whole.index]
ns_whole.sort_values(by='topic', inplace=True)
ns_whole.reset_index(drop=True, inplace=True)

ns = pd.concat([ns_rTPJ, ns_whole], axis=1)

x = np.abs(ns.zReverse) #-np.log(ns.pReverse)
y = np.abs(ns.zForward)#-np.log(ns.pForward)
labs, labsx, labsy = [],[],[]
for i,l in enumerate(ns.topic):
    if ns.pReverse[i] < 0.05 or ns.pForward[i] < 0.05 or ns.q[i] < 0.05:
        lab = '-'.join(l.split('_')[1:])
        if ns.q[i] < 0.05:
            lab = lab+'*'
        labs.append(lab)
        labsx.append(x[i])
        labsy.append(y[i])

sns.set(font_scale=1.5)
sns.set_style("ticks")
fig = plt.figure(figsize=(12,8))
scat = sns.scatterplot(x=x, 
                        y=y,
                        hue=-np.log10(ns.p),
                        size=-np.log10(ns.p),
                        sizes=(20,10000),
                        palette='inferno',
                        alpha=0.4, linewidth=1, edgecolor='k')
leg = plt.legend(bbox_to_anchor=(1.02, 0.95), loc=2, borderaxespad=0, borderpad=0.8, 
           handletextpad=2.2, labelspacing=2.5, frameon=False,
           title="Spearman's rho,\n-log10(p)\n(whole-brain)")
for lh in leg.legendHandles: 
    lh.set_alpha(0.4)
    lh.set_edgecolor('k')
scat.axhline(p_to_z(0.05), c='k', linewidth=1, linestyle='--', alpha=0.8)
scat.axvline(p_to_z(0.05), c='k', linewidth=1, linestyle='--', alpha=0.8)
scat.set_xlabel('Z reverse inference (rTPJ)')
scat.set_ylabel('Z forward likelihood (rTPJ)')

texts=[]
for l,x,y in zip(labs, labsx, labsy):
    texts.append(scat.text(x=x, y=y, s=l, 
                           size=17, horizontalalignment='center', verticalalignment='center'))
adjust_text(texts, expand_text=(1.5, 1.1), arrowprops=dict(arrowstyle="-", color='k', alpha=0.8))

plt.savefig('fig3f.pdf', transparent=True, bbox_inches='tight')




# %%
