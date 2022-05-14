#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 12:30:29 2021

@author: Leon D. Lotter
"""

from nilearn.plotting import plot_glass_brain
from nilearn.datasets import fetch_surf_fsaverage
from neuromaps.transforms import mni152_to_fsaverage
from surfplot import Plot


def plot_gb(img, title=None, thresh=0, col='RdBu_r'):
    """
    Plots volume on a nilearn glass brain.
    
    Input: img=input volume, title=plot title, thresh=plot threshold, col=color
    """
        
    gb = plot_glass_brain(img, title=title, threshold=thresh, draw_cross=False, 
                          cmap=col, display_mode='lyrz')


#=============================================================================


def plot_surf(vol, title='', interp='linear', legend=True, cmap='viridis',
              overlay=None, views=['lateral', 'medial'], size=(1000, 250),
              flip=False, zoom=1.2):

    fsaverage = fetch_surf_fsaverage()
    map_lh, map_rh = mni152_to_fsaverage(vol, method=interp, fsavg_density='10k')

    p = Plot(surf_lh=fsaverage['pial_left'], surf_rh=fsaverage['pial_right'], 
             layout='row', size=size, zoom=zoom, views=views, flip=flip)
    p.add_layer({'left': map_lh, 'right': map_rh}, cmap=cmap, cbar=legend)
    
    if overlay is not None:
        outline_lh, outline_rh = mni152_to_fsaverage(overlay, method='nearest', fsavg_density='10k')
        p.add_layer({'left': outline_lh, 'right': outline_rh}, 
                    cmap='red_transparent', as_outline=True, cbar=False)

    if legend == True:
        kws = dict(location='bottom', draw_border=True, aspect=15, shrink=.2,
                   decimals=0, pad=-0.08)
        fig = p.build(cbar_kws=kws) 
    else:
        fig = p.build()

    fig.axes[0].set_title(title, pad=-3)
    #fig.show()
