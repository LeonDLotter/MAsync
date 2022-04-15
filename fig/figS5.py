# ========================================================================
# Figure S5 - Pet and cell types
# ========================================================================

#%% 
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

from os.path import join
from re import A
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import seaborn as sns


# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values
c_continous = None
c_ale001 = ListedColormap(cmap[0,:])
c_ale01 = ListedColormap(cmap[[0,0,0,0,0,0,0,0,3],:])
c_overlay = ListedColormap('#000000')
c_macm =  ListedColormap(np.concatenate((cmap, cmap)))


def colors_from_values(values, palette_name):
    # normalize the values to range [0, 1]
    normalized = (values - min(values)) / (max(values) - min(values))
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    return np.array(palette).take(indices, axis=0)



# %% FigS5 pet and cell ========================================================================

fsize = 12

# prepare data
# pet
pet = pd.read_csv(join(wd, 'context', 'PETmRNA_ale_z_fdr.csv'))
pet.sort_values(by='zr', ascending=False, inplace=True)
pet.reset_index(drop=True, inplace=True)
# cell
cell = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_ale_z_PsychEncodeTPM.csv'))
cell.sort_values(by='cScorePheno', ascending=False, inplace=True)
cell.reset_index(drop=True, inplace=True)

# PLOT
figS3, axes = plt.subplots(1,2, figsize=(7.5,12))


# PET ========================================================================
pet_plot = sns.barplot(data=pet, x='zr', y='PETmRNA', ax=axes[0],
                       palette=colors_from_values(-np.log(pet.p), "viridis"))
# labels
for i, p in enumerate(pet_plot.patches):
    x = p.get_x() + p.get_width() if pet.zr[i] > 0 else p.get_x()
    xy = (5, -15)
    label = pet.PETmRNA[i]+'*' if pet.sig[i]==True else pet.PETmRNA[i]
    weight = 'bold' if pet.sig[i]==True else 'normal'
    pet_plot.annotate(label, (x, p.get_y()), xytext=xy, textcoords='offset points', size=fsize, weight=weight)
# style
pet_plot.axvline(0, color='k', linewidth=1)
pet_plot.set_xlabel("Z (Spearman's rho)", fontsize=fsize)
pet_plot.set(yticklabels=[], ylabel=None, xlim=(-0.55,0.55))
pet_plot.tick_params(left=False)
# cbar
cb = plt.colorbar(mpl.cm.ScalarMappable(
    norm=mpl.colors.Normalize(vmin=-np.log(pet.p.max()), vmax=-np.log(pet.p.min())), cmap='viridis'),
    orientation='vertical', cax = figS3.add_axes([0.18, 0.6, 0.03, 0.28]))
cb.set_label('- log(p)', size=fsize)


# CELL =======================================================================
cell_plot = sns.barplot(data=cell, x='cScorePheno', y='cLabel', ax=axes[1],
                        palette=colors_from_values(-np.log(cell.pValZ), "cividis"))
# labels
for i, p in enumerate(cell_plot.patches):
    x = p.get_x() + p.get_width() if cell.cScorePheno[i] > 0 else p.get_x()
    xy = (5, -15)
    label = cell.cLabel[i]+'*' if cell.pValPermCorr[i]<0.05 else cell.cLabel[i]
    weight = 'bold' if cell.pValPermCorr[i]<0.05 else 'normal'
    cell_plot.annotate(label, (x, p.get_y()), xytext=xy, textcoords='offset points', size=fsize, weight=weight)
# style
cell_plot.axvline(0, color='k', linewidth=1)
cell_plot.set_xlabel("Average Z (Spearman's rho)", fontsize=fsize)
cell_plot.set(yticklabels=[], ylabel=None, xlim=(-0.4,0.4))
cell_plot.tick_params(left=False)
# cbar
cb = plt.colorbar(mpl.cm.ScalarMappable(
    norm=mpl.colors.Normalize(vmin=-np.log(cell.pValZ.max()), vmax=-np.log(cell.pValZ.min())), cmap='cividis'),
    orientation='vertical', cax = figS3.add_axes([0.59, 0.6, 0.03, 0.28]))
cb.set_label('- log(p)', size=fsize)


sns.despine(figS3, None, True, True, True, False)
plt.savefig(join(sdir, 'figS5.pdf'), transparent=True, bbox_inches='tight')


# %% prelim plot brains? =============================================================

import sys
sys.path.append(join('/Users/leonlotter/MAsync/project', 'src'))
from utils_image import parcel_data_to_volume
from nilearn.input_data import NiftiLabelsMasker
from nilearn.plotting import plot_glass_brain
from nilearn.image import load_img
from matplotlib.cm import get_cmap

c = get_cmap('viridis').colors
atlas = join(wd, 'atlases', 'Schaefer100-7_TianS1_2mm.nii.gz')

# ale volume
ale = pd.read_csv(join(wd, 'macm', 'cor_macm_ale.csv'))['dat1']
ale_vol, _ = parcel_data_to_volume(ale, atlas, rank=True)
ale_cl = load_img(join(wd, 'ale', 'ale_thresh_bin.nii.gz'))

# pet 
masker = NiftiLabelsMasker(atlas)
gaba = masker.fit_transform(join(wd, 'context', 'parcel_maps', 'GABAa.nii.gz'))[0] 
ser = masker.fit_transform(join(wd, 'context', 'parcel_maps', '5HT2a.nii.gz'))[0] 
gaba_vol, _ = parcel_data_to_volume(gaba, atlas, rank=True)
ser_vol, _ = parcel_data_to_volume(ser, atlas, rank=True)

# cells
ex3 = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_cell_genes_Adult-Ex3.csv')).mean(axis=1)
ex3_vol, _ = parcel_data_to_volume(ex3, atlas, rank=True)
in6 = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_cell_genes_Adult-In6.csv')).mean(axis=1)
in6_vol, _ = parcel_data_to_volume(in6, atlas, rank=True)
in5 = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_cell_genes_Adult-In5.csv')).mean(axis=1)
in5_vol, _ = parcel_data_to_volume(in5, atlas, rank=True)

fig, axes = plt.subplots(6,1, figsize=(9.5,16))
plt.subplots_adjust(hspace=0)
for i,(v,l) in enumerate(zip([ale_vol, gaba_vol, ser_vol, ex3_vol, in6_vol, in5_vol],
                         ['INS', 'GABAa', '5HT2a', 'Excitatory 3', 'Inhibitory 6', 'Inhibitory 5'])):
    gb = plot_glass_brain(v, title=l, colorbar=False, cmap=ListedColormap(np.concatenate((c,c))), 
                     display_mode='lyrz', figure=fig, axes=axes[i])
    gb.add_contours(ale_cl, colors='k', linewidths=0.7)    

plt.savefig(join(sdir, 'pet_cell_brains.pdf'), transparent=True, bbox_inches='tight')


# %%
