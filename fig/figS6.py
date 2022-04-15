# ========================================================================
# Figure S6 - GABA validation
# ========================================================================

#%% 
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import seaborn as sns

# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values


#%% cluster-heatmap

gaba_all = pd.read_csv(join(wd, 'context', 'GABAall_parcellated_data.csv'))
gaba_all_cor = gaba_all.corr(method='spearman')

figS6a = sns.clustermap(gaba_all_cor, figsize=(10,10))

figS6a.savefig(join(sdir, 'figS6a.pdf'), transparent=True, bbox_inches='tight') 


# %% scatters

gaba_cl = pd.read_csv(join(wd, 'context', 'GABAcl_parcellated_data.csv'))
ale_parc = pd.read_csv(join(wd, 'macm', 'cor_macm_ale.csv'))['dat1']
ale_parc.name = 'INS ALE Z-scores'

figS6b = plt.figure(figsize=(3,12), constrained_layout=False)
gs = figS6b.add_gridspec(ncols=1, nrows=4, hspace=0.3, wspace=0, width_ratios=[1])
for i, cl in enumerate(gaba_cl.columns):
    ax = figS6b.add_subplot(gs[i,0])
    sns.regplot(x=ale_parc, y=gaba_cl[cl], ax=ax, 
                scatter_kws={'edgecolor':'k', 'alpha':0.5, 'color':to_hex(cmap[1,:])})
    ax.set_ylabel(f'GABA cluster {i+1}')

figS6b.savefig(join(sdir, 'figS6b.pdf'), transparent=True, bbox_inches='tight') 


# %%
