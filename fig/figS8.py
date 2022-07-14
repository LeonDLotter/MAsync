
# ========================================================================
# Figure S8 - Dominance Analysis
# ========================================================================

#%% 
wd = '/Users/llotter/MAsync/project/data'
sdir = '/Users/llotter/MAsync/project/fig'

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

def colors_from_values(values, palette_name):
    # normalize the values to range [0, 1]
    normalized = (values - min(values)) / (max(values) - min(values))
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    return np.array(palette).take(indices, axis=0)

cmap = "viridis"

#%% get data & plot
da_data = pd.read_csv(join(wd, "context", "PETmRNAcells_dominance_analysis.csv"), index_col=0)
da_totalr2 = da_data["total dominance"].sum()
predictors = list(da_data.index)
#colors = sns.color_palette("tab10", len(da_data))
#colors.reverse()

fig, axes = plt.subplots(3,1, figsize=(7,9), 
                       gridspec_kw=dict(height_ratios=(0.13,1,1)),
                       sharex=False)
axes = axes.ravel()    

# 1: full model
predictors.reverse()
colors = np.flip(colors_from_values(da_data["total dominance"], cmap),0)
perc = da_totalr2
for i, predictor in enumerate(predictors):
    perc_df = pd.DataFrame(
        data={"Total Dominance": perc},
        index=["Full Model"]
    )
    sns.barplot(data=perc_df,
                y=perc_df.index,
                x="Total Dominance",
                ax=axes[0],
                color=colors[i],
                linewidth=1,
                edgecolor="w")
    perc = perc - da_data["total dominance"][predictor]
axes[0].set_xlim(0, 0.4)
axes[0].set_xlabel("Explained Variance", size=12)
axes[0].set_yticklabels(["Full Model"], size=12)
p = axes[0].patches[0]
axes[0].annotate("100.0%", (p.get_x() + p.get_width(), p.get_y()), xytext=(5, -15), textcoords='offset points', size=11)

# 2: total dominance    
sns.barplot(data=da_data,
            y=da_data.index,
            x="total dominance",
            ax=axes[1],
            palette=colors_from_values(da_data["total dominance"], cmap),
            linewidth=1,
            edgecolor="w")
# details
for i, p in enumerate(axes[1].patches):
    x = p.get_x() + p.get_width()
    perc = da_data["total dominance"][i] / da_totalr2 * 100
    label = f"{perc:.01f}%" 
    axes[1].annotate(label, (x, p.get_y()), xytext=(5, -15), textcoords='offset points', size=11)
axes[1].set_xlim(0, 0.4)
axes[1].set_yticklabels(axes[1].get_yticklabels(), size=12)
axes[1].set_xlabel("Explained Variance (Total Dominance)", size=12)

# 3: individual dominance
sns.barplot(data=da_data,
            y=da_data.index,
            x="individual dominance",
            ax=axes[2],
            palette=colors_from_values(da_data["individual dominance"], cmap),
            linewidth=1,
            edgecolor="w")
# details
axes[2].set_xlim(0, 0.4)
axes[2].set_yticklabels(axes[1].get_yticklabels(), size=12)
axes[2].set_xlabel("Explained Variance (Individual Dominance)", size=12)

# save
fig.tight_layout()
fig.savefig(join(sdir, "figS8.pdf"), transparent=True, bbox_inches='tight')
# %%
