# ========================================================================
# Figure 3 - Molecular Context
# ========================================================================

#%% 
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import seaborn as sns

def colors_from_values(values, palette_name):
    # normalize the values to range [0, 1]
    normalized = (values - min(values)) / (max(values) - min(values))
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    return np.array(palette).take(indices, axis=0)


# %% Fig 3A scatters ========================================================================

# pet data
pet_corr = pd.read_csv(join(wd, 'context', 'PETmRNA_ale_z_fdr.csv'))
pet_sig = list(pet_corr.query('(zr>0) & (sig==True)').PETmRNA)
pet_zr = list(pet_corr.query('(zr>0) & (sig==True)').zr)
pet_p = list(pet_corr.query('(zr>0) & (sig==True)').p)

# parcellated data
ale_parc = pd.read_csv(join(wd, 'macm', 'cor_macm_ale.csv'))['dat1']
ale_parc.name = 'ALE Z'
pet_parc = pd.read_csv(join(wd, 'datasets', 'pet_parcellated_data.csv'))
fnirs_data = pd.read_csv(join(wd, 'fnirs', 'fnirs_atlas_result.csv'))
fnirs_vals = np.zeros((116))
for i in range(116):
    if i+1 in list(fnirs_data.region_idx):
        fnirs_vals[i] = fnirs_data[fnirs_data.region_idx==i+1].p_ch_ratio_sub_all.values
    else:
        fnirs_vals[i] = np.nan
fnirs_nan_mask = np.isnan(fnirs_vals)

# make plot
fig3a, axes = plt.subplots(len(pet_sig),1, figsize=(3.2,11), sharex=True)
# iterate axes
for i, pet in enumerate(pet_sig):
    sns.scatterplot(x=ale_parc, y=pet_parc[pet], ax=axes[i], 
                    hue=-np.log10(fnirs_vals), size=-np.log10(fnirs_vals), palette='viridis',
                    edgecolor='k', alpha=0.6, sizes=(25,100))
    sns.scatterplot(x=ale_parc[fnirs_nan_mask], y=pet_parc[pet][fnirs_nan_mask], ax=axes[i],
                    edgecolor='k', alpha=0.9, color='w', size=3, legend=None)
    sns.regplot(x=ale_parc, y=pet_parc[pet], ax=axes[i], scatter=None)
    axes[i].set_ylabel(f'{pet} density (Z)', fontsize=12)
    axes[i].set_xlabel(None)
    axes[i].legend().set_visible(False)
    p = f'= {pet_p[i]:.03f}' if round(pet_p[i],3) > 0 else '< 0.001'
    axes[i].annotate(f'Z(rho) = {pet_zr[i]:.2f}, p {p}', 
                     xy=(0.97,0.05), xycoords='axes fraction', ha='right')
# set y label for bottom plot
axes[i].set_xlabel('INS ALE (Z)', fontsize=12)
# join legend
handles, labels = axes[i].get_legend_handles_labels()
fig3a.legend(handles, labels, title='fNIRS\n-log10(p)', title_fontsize=12, 
             labelspacing=1,
             loc='upper right',  bbox_to_anchor=(1.3, 0.99))

plt.tight_layout()
plt.savefig(join(sdir, 'fig3a.pdf'), transparent=False, bbox_inches='tight')


# %% Fig 3B GCEA cells & disease =====================================================================

# data
cell = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_ale_z_PsychEncodeTPM.csv'))
cell.sort_values(by='cScorePheno', ascending=False, inplace=True)
cell.reset_index(drop=True, inplace=True)
disease = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_ale_z_DisGeNET.csv'))
disease = disease[disease.pValPermCorr < 0.05]
disease.cDesc1[disease.cDesc1=='MAJOR AFFECTIVE DISORDER 2'] = 'Major Affective Disorder 2'
disease.sort_values(by='cScorePheno', ascending=False, inplace=True)
disease.reset_index(drop=True, inplace=True)

# plot function
def plot_bars(data, ax, y, xlim, color, mark_labels=True):
    # plot
    sns.barplot(data=data, x='cScorePheno', y=y, ax=ax,
                palette=colors_from_values(-np.log(data.pValZ), color))
    # make labels nice
    for i, p in enumerate(ax.patches):
        # category label
        x = p.get_x() + p.get_width() if data.cScorePheno[i] > 0 else p.get_x()
        xy = (5, -12)
        label = data[y][i]+'*' if (data.pValPermCorr[i]<0.05) & mark_labels else data[y][i]
        weight = 'bold' if (data.pValPermCorr[i]<0.05) & mark_labels else 'normal'
        ax.annotate(label, (x, p.get_y()), xytext=xy, textcoords='offset points', size=11, weight=weight)
        # category size
        #x = p.get_x() if data.cScorePheno[i] > 0 else p.get_x() + p.get_width() 
        #xy = (-5, -12)
        #label = data['cSize'][i]
        #ax.annotate(label, (x, p.get_y()), xytext=xy, textcoords='offset points', size=11, ha='right')
    # add details
    ax.axvline(0, color='k', linewidth=1)
    ax.set_xlabel("Average Z (Spearman's rho)", fontsize=11)
    ax.set(yticklabels=[], ylabel=None, xlim=xlim)
    ax.tick_params(left=False)
    # add colorbar
    cb = plt.colorbar(
        mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(
                vmin=-np.log10(data.pValZ.max()), 
                vmax=-np.log10(data.pValZ.min())), 
            cmap=color),
        orientation='horizontal', 
        ax=ax,
        location='bottom',
        pad=0.06) 
    cb.set_label('- log10(p)', size=11)

# make figure
fig3b, axes = plt.subplots(1,2, figsize=(8,12), gridspec_kw=dict(wspace=0.7))
# cell plot
plot_bars(data=cell, ax=axes[0], y='cLabel', xlim=(-0.3,0.5), color='cividis')
# disease plot
plot_bars(data=disease, ax=axes[1], y='cDesc1', xlim=(0,0.5), color='magma', mark_labels=False)
# remove y axes
sns.despine(fig3b, None, True, True, True, False)
plt.savefig(join(sdir, 'fig3b.pdf'), transparent=True, bbox_inches='tight')


# %% Fig3C brainspan =======================================

from heatmap import heatmap
import seaborn as sns
sns.set(color_codes=True, font_scale=1.1)

devel = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_ale_z_BrainSpan.csv'))
devel['age'] = [l.split('_')[0] for l in devel.cLabel]
devel['region'] = [l.split('_')[1] for l in devel.cLabel]

age_order = ['adult', 'adolescent', 'child', 'infant', 'prenatal']
region_order = ['OFC', 'MFC', 'DFC', 'VFC', 'M1C', 'S1C', 'STC', 'ITC', 'IPC', 
                'A1C', 'V1C', 'HIP', 'MD', 'STR', 'AMY', 'CBC']

plt.figure(figsize=(8, 2.5))
heatmap(x=devel['region'], y=devel['age'], y_order=age_order, x_order=region_order,
        color=-np.log10(devel['pValZ']), size=devel['cScorePheno'], marker='o',
        palette=sns.color_palette('magma_r', 500)[::-1], 
        color_range=(0, round((-np.log10(devel['pValZ'])).max(),1)))
plt.annotate('test', (2.2,0.5))
plt.savefig(join(sdir, 'fig3c.pdf'), transparent=True, bbox_inches='tight')



