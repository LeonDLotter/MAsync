#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 18:28:38 2021

@author: Leon D. Lotter
"""
#=============================================================================


def close_line(line):
    # helper function for radar plots
    
    import numpy as np
    
    x,y = line.get_data()
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    line.set_data(x, y)


def set_rotated_labels(labels, fig, ax, y_dist=0, fdr=None, unc=None):
    # helper function for radar plots
    
    import numpy as np
    
    # if no significance given
    if fdr is None:
        fdr = [0 for i in range(len(labels))]
    if unc is None:
        unc = [0 for i in range(len(labels))]
    
    # set ticks
    ticks = np.linspace(0, 360, len(labels)+1)[:-1] 
    ax.set_xticks(np.deg2rad(ticks))
    ax.set_xticklabels(labels, fontsize=10)
    # get rotation angles
    angles = np.linspace(0,2*np.pi, len(ax.get_xticklabels())+1)
    angles[np.cos(angles) < 0] = angles[np.cos(angles) < 0] + np.pi
    angles = np.rad2deg(angles)
    # rotate labels & set colors according to significance
    fig.canvas.draw
    labels = []
    for label, theta, f, u in zip(ax.get_xticklabels(), angles, fdr, unc):
        x, y = label.get_position()        
        ha = "left" if x <= np.pi/2 or x > np.pi*1.5 else "right"
        color = "darkred" if f == 1 or u == 1 else "black"
        weight = "bold" if f == 1 else "normal"
        lab = ax.text(x, y+y_dist, label.get_text(), transform=label.get_transform(),
                      ha=ha, va=label.get_va(), rotation_mode="anchor", 
                      c=color, fontweight=weight)
        lab.set_rotation(theta)
        labels.append(lab)
    ax.set_xticklabels([])
    
    
def make_radarplot(labels, vals, fowRev, labs, colors, alpha, 
                   tickLoc=45, unc=None, fdr=None, leg_inf=True, leg_dom=True):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    sns.set_style("white")
    
    # label positions
    N = len(labels)
    theta = [n / float(N) * 2 * np.pi for n in range(N)]
    
    # create figure
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(polar=True, zorder=1)
    
    # ri / fi estimates
    val_plots = []
    val_fills = []
    for i,val in enumerate(vals):
        plot = ax.plot(theta, val, linewidth=1.5, color=colors[i], alpha=0.75, zorder=3+i)
        close_line(plot[0])
        fill = ax.fill(theta, val, facecolor=colors[i], alpha=alpha[i], zorder=3+i, )
        val_plots.append(plot)
        val_fills.append(fill)
        
    # y axis
    ax.set_rlabel_position(tickLoc)
    
    # x axis rotation
    set_rotated_labels(labels, fig, ax, fdr=fdr, unc=unc)
    
    # add category pie
    categ = [8,17,15,11,9]
    categ_lab = ["Action", "Cognition", "Emotion", "Interoception", "Perception"]
    ax2 = fig.add_subplot(111, label="pie axes", zorder=0)
    categ_pie,_,_ = ax2.pie(categ, radius=1.32, autopct='%1.1f%%', startangle=-3)
    
    # add ri / fi legend
    if leg_inf is True:
        if fowRev == "fi":
            lab_title = "$\\bf{Forward\ Inference}$" "\n" "${P(Activation\ |\ Domain)}$"
        elif fowRev == "ri":
            lab_title = "$\\bf{Reverse\ Inference}$" "\n" "${P(Domain\ |\ Activation)}$"
        leg1 = ax.legend([v[0] for v in val_plots], labs, loc="lower left", bbox_to_anchor=(-0.51, -0.5),
                   title=lab_title)
        # increase legend linewidth
        for line in leg1.get_lines():
            line.set_linewidth(3.0)
        
    # add category legend
    if leg_dom is True:
        ax2.legend(categ_pie, categ_lab, loc="lower right", bbox_to_anchor=(1.38, -0.5),
                   title="$\\bf{BrainMap\ Domain}$")
    
    # adjust size
    top = 0.72
    right = 0.95
    plt.subplots_adjust(top=top, bottom=1-top, left=1-right, right=right)

    return (fig)
    
