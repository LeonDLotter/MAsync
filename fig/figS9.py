# ========================================================================
# Figure S9 - GO GCEA
# ========================================================================

#%% 
wd = '/Users/leonlotter/MAsync/project/data'
sdir = '/Users/leonlotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from wordcloud import WordCloud 

# %% FigS8b wordclouds ========================================================================

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

plt.savefig(join(sdir, 'figS9b.png'), dpi=300, transparent=False, bbox_inches='tight')


# %%
