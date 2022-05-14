# ========================================================================
# Figure S2 - single fNIRS studies 
# ========================================================================

# %%
pdir = '/Users/leonlotter/MAsync/project'
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

import sys
from os.path import join
import pandas as pd
sys.path.append(pdir)
from src.utils_io import fnirs_to_nimare_dataset
from nimare.meta.cbma.ale import ALE
import matplotlib.pyplot as plt
from nilearn.plotting import plot_glass_brain


# spreadsheet with fnirs results
fnirs_tab = pd.read_csv(join(wd, 'datasets', 'fnirs_coordinates.csv'))
ds = fnirs_to_nimare_dataset(fnirs_tab, False)
ids = ds.ids
exps = ds.metadata.study_id

# %% get plot

# ALE
ale = ALE() # get ALE class

# FIGURE
figS1, axes = plt.subplots(20,3, figsize=(25,45), constrained_layout=False)
axes = axes.ravel()
plt.subplots_adjust(wspace=0.05, hspace=0)

for i, id in enumerate(ids):
    # current exp
    ds_exp = ds.slice([id])
    print(f'{i+1}/{len(ids)}: {id}')
    # get map
    ds_ale = ale.fit(ds_exp) # fit to dataset
    z_exp = ds_ale.get_map('z')
    # plot
    plot_glass_brain(z_exp, figure=figS1, axes=axes[i], display_mode='lyrz',
                     title=f'{exps[i]}: {len(ds_exp.coordinates)} channel(s), {ds_exp.metadata.sample_sizes.values[0][0]:2.0f} subjects')

plt.savefig(join(sdir, 'figS2.pdf'), transparent=False, bbox_inches='tight')
plt.savefig(join(sdir, 'figS2.png'), transparent=False, bbox_inches='tight')


# %%
