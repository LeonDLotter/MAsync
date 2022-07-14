# ========================================================================
# Figure S6 - Pet and cell types
# ========================================================================

#%% 
wd = '/Users/llotter/MAsync/project/data'
sdir = '/Users/llotter/MAsync/project/fig'

from os.path import join
from re import A
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import seaborn as sns
import sys
sys.path.append(join('/Users/llotter/MAsync/project', 'src'))
from utils_image import parcel_data_to_volume
from nilearn.input_data import NiftiLabelsMasker
from nilearn.plotting import plot_glass_brain
from nilearn.image import load_img
from matplotlib.cm import get_cmap

# get atlas
atlas = join(wd, 'atlases', 'Schaefer100-7_TianS1_2mm.nii.gz')

# prepare data
pet = pd.read_csv(join(wd, 'context', 'PETmRNA_ale_z.csv'), index_col=0)
pet.sort_values(by='zr', ascending=False, inplace=True)
pet.reset_index(drop=False, inplace=True)
pet_parc = pd.read_csv(join(wd, 'datasets', 'pet_parcellated_data.csv'))

# brain color
c = get_cmap('viridis').colors

def colors_from_values(values, palette_name):
    # normalize the values to range [0, 1]
    normalized = (values - min(values)) / (max(values) - min(values))
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    return np.array(palette).take(indices, axis=0)


# %% FigS6A pet ========================================================================

fsize = 12

# PLOT
figS6, ax = plt.subplots(1,1, figsize=(4,12))

pet_plot = sns.barplot(data=pet, x='zr', y='index', ax=ax,
                       palette=colors_from_values(-np.log10(pet.p), "viridis"))
# labels
for i, p in enumerate(pet_plot.patches):
    x = p.get_x() + p.get_width() if pet.zr[i] > 0 else p.get_x()
    xy = (5, -15)
    label = pet["index"][i]+'*' if pet.q[i]<0.05 else pet["index"][i]
    weight = 'bold' if pet.p[i]<0.05 else 'normal'
    pet_plot.annotate(label, (x, p.get_y()), xytext=xy, textcoords='offset points', size=fsize, weight=weight)
# style
pet_plot.axvline(0, color='k', linewidth=1)
pet_plot.set_xlabel("Z (Spearman's rho)", fontsize=fsize)
pet_plot.set(yticklabels=[], ylabel=None, xlim=(-0.55,0.55))
pet_plot.tick_params(left=False)
# cbar
cb = plt.colorbar(mpl.cm.ScalarMappable(
    norm=mpl.colors.Normalize(vmin=-np.log10(pet.p.max()), vmax=-np.log10(pet.p.min())), cmap='viridis'),
    orientation='vertical', cax = figS6.add_axes([0.18, 0.6, 0.03, 0.28]))
cb.set_label('- log10(p)', size=fsize)

sns.despine(figS6, None, True, True, True, False)
plt.savefig(join(sdir, 'figS6a.pdf'), transparent=True, bbox_inches='tight')


# %% FigS6B PET brains =============================================================

# ale volume
ale_data = pd.read_csv(join(wd, 'ale', 'ale_z_parc.csv'))
ale_vol = parcel_data_to_volume(ale_data["ALE z"], atlas, rank=False)[0]
ale_vol_rank = parcel_data_to_volume(ale_data["ALE z"], atlas, rank=True)[0]
# ale cluster
ale_cl = load_img(join(wd, 'ale', 'ale_thresh_bin.nii.gz'))

# pet 
pet_labels = pet.query('p < 0.05')["index"].to_list()
pet_maps = list()
pet_maps_rank = list()
for map in pet_labels:
    pet_vol = parcel_data_to_volume(pet_parc[map], atlas, rank=False)[0]
    pet_vol_rank = parcel_data_to_volume(pet_parc[map], atlas, rank=True)[0]
    pet_maps.append(pet_vol)
    pet_maps_rank.append(pet_vol_rank)

# combine
maps_all = [ale_vol] + pet_maps
maps_rank_all = [ale_vol_rank] + pet_maps_rank
labels_all = ['INS'] + pet_labels

fig, axes = plt.subplots(6,2, figsize=(20,14))
axes = axes.ravel()
plt.subplots_adjust(hspace=0)
ax=0
for i, map in enumerate(maps_all):
    # real map
    gb = plot_glass_brain(map, title=f'{labels_all[i]}: parcellated', colorbar=False, cmap=ListedColormap(np.concatenate((c,c))), 
                          display_mode='lyrz', figure=fig, axes=axes[ax])
    gb.add_contours(ale_cl, colors='k', linewidths=0.7)    
    # ranked map
    gb = plot_glass_brain(maps_rank_all[i], title=f'{labels_all[i]}: parcellated & ranked', colorbar=False, cmap=ListedColormap(np.concatenate((c,c))), 
                          display_mode='lyrz', figure=fig, axes=axes[ax+1])
    gb.add_contours(ale_cl, colors='k', linewidths=0.7)   
    ax += 2 

plt.savefig(join(sdir, 'figS6b.pdf'), transparent=True, bbox_inches='tight')


# %%
