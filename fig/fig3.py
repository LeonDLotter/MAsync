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
from matplotlib.colors import ListedColormap, to_hex
import seaborn as sns
from wordcloud import WordCloud 

# colors
cmap = pd.read_csv(join(sdir, 'colors.csv'), header=None).values


# %% Fig 3A scatters ========================================================================

pet_parc = pd.read_csv(join(wd, 'datasets', 'pet_parcellated_data.csv'))
ale_parc = pd.read_csv(join(wd, 'macm', 'cor_macm_ale.csv'))['dat1']
ale_parc.name = 'ALE Z'

reglkws = {'color':'#3E6282'}
scatkws = {'edgecolor':'k'}
fsize=12

fig3a, axes = plt.subplots(1,2, figsize=(10,3.5))
axes = axes.ravel()
# GABA
sns.regplot(x=ale_parc, y=pet_parc['GABAa'], ax=axes[0], line_kws={'color':cmap[1,:]}, 
            scatter_kws={'edgecolor':'k', 'alpha':0.5, 'color':to_hex(cmap[1,:])})
axes[0].set_ylabel('GABAa receptor availability', fontsize=fsize)
axes[0].set_xlabel('INS ALE Z-scores', fontsize=fsize)
# 5HT2A
sns.regplot(x=ale_parc, y=pet_parc['5HT2a'], ax=axes[1], line_kws={'color':cmap[0,:]}, 
            scatter_kws={'edgecolor':'k', 'alpha':0.5, 'color':to_hex(cmap[0,:])})
axes[1].set_ylabel('5HT2a receptor availability', fontsize=fsize)
axes[1].set_xlabel('INS ALE Z-scores', fontsize=fsize)

plt.savefig(join(sdir, 'fig3a.pdf'), transparent=False, bbox_inches='tight')


# %% Fig 3B scatters ========================================================================

ale_parc = pd.read_csv(join(wd, 'macm', 'cor_macm_ale.csv'))['dat1']
ale_parc.name = 'ALE Z'

cell_parc_ex3 = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_cell_genes_Adult-Ex3.csv'))
cell_parc_ex3.columns = [f'Ex3: {c}' for c in cell_parc_ex3.columns]
df_ex3 = pd.melt(pd.concat([cell_parc_ex3, ale_parc], axis=1), id_vars=['ALE Z'], var_name='Gene', value_name='Gene expression')
df_ex3['Cell type'] = 'Ex3'

cell_parc_in6 = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_cell_genes_Adult-In6.csv'))
cell_parc_in6.columns = [f'In6: {c}' for c in cell_parc_in6.columns]
df_in6 = pd.melt(pd.concat([cell_parc_in6, ale_parc], axis=1), id_vars=['ALE Z'], var_name='Gene', value_name='Gene expression')
df_in6['Cell type'] = 'In6'

cell_parc_in5 = pd.read_csv(join(wd, 'context', 'gcea', 'GCEA_cell_genes_Adult-In5.csv'))
cell_parc_in5.columns = [f'In5: {c}' for c in cell_parc_in5.columns]
df_in5 = pd.melt(pd.concat([cell_parc_in5, ale_parc], axis=1), id_vars=['ALE Z'], var_name='Gene', value_name='Gene expression')
df_in5['Cell type'] = 'In5'

fackws = {'sharey':True, 'sharex':True, 'despine':False}

cells = sns.lmplot(data=pd.concat([df_ex3, df_in6, df_in5], axis=0), col_wrap=7, height=1.3, aspect=0.95, 
                   x='ALE Z', y='Gene expression', col='Gene', hue='Cell type', legend=False,
                   palette=dict(Ex3=cmap[0,:], In6=cmap[1,:], In5=cmap[3,:]), 
                   scatter_kws={'edgecolor':'k', 'alpha':0.1, 's':15}, facet_kws=fackws)
cells.set_titles("{col_name}", fontweight='bold')
cells.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
cells.fig.text(0.52,0.07, 'INS ALE Z-scores ', fontsize=13, rotation='horizontal', ha='center')
cells.fig.text(0.015,0.54, 'Gene expression', fontsize=13, rotation='vertical', va='center')
cells.savefig(join(sdir, 'fig3b.pdf'), transparent=False, bbox_inches='tight')


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
        color=-np.log(devel['pValZ']), size=devel['cScorePheno'], marker='o',
        palette=sns.color_palette('magma_r', 500)[::-1], color_range=(0,12))
plt.annotate('test', (2.2,0.5))
plt.savefig(join(sdir, 'fig3c.pdf'), transparent=True, bbox_inches='tight')


# %% Fig3D2 wordclouds ========================================================================

gofigure_terms = pd.read_csv(join(sdir, 'gofigure', 'gofigure_res_descriptions_5.csv'))

fig3b2, axes = plt.subplots(5,1, figsize=(4,4), constrained_layout=True)
axes = axes.ravel()

cloud = WordCloud(width=1500, height=250, margin=1, background_color='white', colormap='magma', 
                    stopwords=['regulation', 'of', 'positive', 'negative', 'in', 'for', 'to'])
for i, terms in enumerate(gofigure_terms.descs_str):
    cloud.generate(terms)
    axes[i].imshow(cloud, interpolation="bilinear")
    axes[i].axis("off")
    axes[i].set_title(str(i+1), x=-0.05, y=0.3, fontname='Arial', fontweight='bold')

plt.savefig(join(sdir, 'fig3d2.png'), dpi=300, transparent=False, bbox_inches='tight')


# %%

# %%
