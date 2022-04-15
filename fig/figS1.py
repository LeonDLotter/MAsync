# ========================================================================
# Figure S1 - single fMRI studies 
# ========================================================================

# %%
pd = '/Users/leonlotter/MAsync/project'
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

import sys
from os.path import join
sys.path.append(join(pd, 'src'))
from utils_io import csv_to_nimare_ds
from nimare.meta.cbma.ale import ALE
import matplotlib.pyplot as plt
from nilearn.plotting import plot_glass_brain


ds = csv_to_nimare_ds(file   =join(wd, 'datasets', 'fMRI_coordinates.csv'), # spreadsheed with study info
                      exp_var='publication', # column indicating experiment names
                      n_var  ='n', # column indicating sample sizes
                      con_var='contrasts', # column indicating contrast names
                      spa_var='space', # column indicating coordinate space
                      x_var  ='x', y_var='y', z_var='z', # columns indicating coordinates
                      single_contrasts=False) # concatenate contrasts per experiment
ids = ds.ids
exps = ds.metadata.study_id

# %%

# ALE
ale = ALE() # get ALE class

# FIGURE
figS1, axes = plt.subplots(11,2, figsize=(20,30), constrained_layout=False)
axes = axes.ravel()
plt.subplots_adjust(wspace=0.05, hspace=0)

for i, id in enumerate(ids):
    # current exp
    ds_exp = ds.slice([id])
    print(id)
    # get map
    ds_ale = ale.fit(ds_exp) # fit to dataset
    z_exp = ds_ale.get_map('z')
    # plot
    plot_glass_brain(z_exp, figure=figS1, axes=axes[i], display_mode='lyrz',
                     title=f'{exps[i]}: {len(ds_exp.coordinates)} foci, {ds_exp.metadata.sample_sizes.values[0][0]} subjects')

plt.savefig(join(sdir, 'figS1.pdf'), transparent=False, bbox_inches='tight')
plt.savefig(join(sdir, 'figS1.png'), transparent=False, bbox_inches='tight')


# %%
