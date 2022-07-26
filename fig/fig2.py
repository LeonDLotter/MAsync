# ========================================================================
# Figure 2 - (Prisma) + citation netowork
# ========================================================================

# %%
pp = '/Users/llotter/MAsync/project'
wd = '/Users/llotter/MAsync/project/data'
sdir = '/Users/llotter/MAsync/project/fig'

import sys
from os.path import join
sys.path.append(join(pp, 'src'))
from utils_io import csv_to_nimare_ds, fnirs_to_nimare_dataset
from citenet import citeNet
import pandas as pd
import numpy as np
import json
import networkx
from pyvis.network import Network
import matplotlib.pyplot as plt
from nilearn.plotting import plot_markers, plot_connectome
from matplotlib import cm
from matplotlib.colors import to_hex

# colors
cmap = list(cm.get_cmap('tab10').colors)
col_fMRI = [to_hex(cmap[3]), to_hex(cmap[3])+'cc']
col_fNIRS = [to_hex(cmap[0]), to_hex(cmap[0])+'cc']


# %% get data

# nimare ds
nim_fMRI = csv_to_nimare_ds(join(wd, 'datasets', 'fMRI_coordinates.csv'), # spreadsheed with study info
                      exp_var='publication', # column indicating experiment names
                      n_var  ='n', # column indicating sample sizes
                      con_var='contrasts', # column indicating contrast names
                      spa_var='space', # column indicating coordinate space
                      x_var  ='x', y_var='y', z_var='z', # columns indicating coordinates
                      single_contrasts=False) # concatenate contrasts per experiment
nim_fNIRS = fnirs_to_nimare_dataset(join(wd, 'datasets', 'fNIRS_coordinates.csv'), False)
nim_ds = nim_fMRI.merge(nim_fNIRS)

# metadata
ds_fMRI = pd.read_csv(join(wd, 'datasets', 'fMRI_coordinates.csv'))
ds_fMRI['type'] = 'fMRI'
ds_fMRI.reset_index(drop=True, inplace=True)
ds_fNIRS = pd.read_csv(join(wd, 'datasets', 'fNIRS_coordinates.csv'))
ds_fNIRS.rename(columns={'contrast': 'contrasts'}, inplace=True)
ds_fNIRS['type'] = 'fNIRS'
ds_fNIRS.reset_index(drop=True, inplace=True)
ds_merged = pd.concat([ds_fMRI, ds_fNIRS], ) 

# get dois as identifier
dois = ds_merged.doi.unique().tolist()
dois = [d for d in dois if str(d)!='nan']

# general data
study_info = pd.read_csv(join(sdir, 'citenet', 'included_studies_data.csv'))

# create nw data dict
nw_data = {}
for n_doi, doi in enumerate(dois):
    print(doi)
    # get doi
    ds_doi = ds_merged[ds_merged.doi==doi]
    # references
    references = citeNet.get_references_from_doi(doi, which_list="references"),
    citations = citeNet.get_references_from_doi(doi, which_list="citations")
    # coordinates

    # save data
    nw_data[doi] = {
        'title': study_info.Title[study_info.DOI==doi].values[0],
        'author': study_info.Author[study_info.DOI==doi].values[0],
        'year': study_info.Year[study_info.DOI==doi].values[0],
        'type': ds_doi.type.unique()[0],
        'reference': [r for r in references[0]],
        'citation': citations, 
        'citation_count': len(citations),
        'exps': list(),
        'foci': list()
    }
    # get & save experiment(s)
    exps_doi = ds_doi.publication.unique()
    for exp in exps_doi:
        ds_exp = ds_doi[ds_doi.publication==exp]
        fem = ds_exp.propFem.unique()[0] if ds_exp.type.unique()[0]=='fMRI' else ds_exp.nFem.unique()[0]/ds_exp.n.unique()[0]
        nw_data[doi]['exps'].append({
            'exp': exp,
            'area': ds_exp.area.unique()[0],
            'n': ds_exp.n.unique()[0],
            'nDyads': ds_exp.nDyads.unique()[0],
            'propFem': fem,
            'ageMean': ds_exp.ageMean.unique()[0],
            'task': ds_exp.task.unique()[0],
            'contrast': ds_exp.contrasts.unique().tolist()
        })
        nw_data[doi]['foci'].append(
            nim_ds.coordinates[nim_ds.coordinates.study_id==exp][['x','y','z']].values)
    # create glassbrain
    # coords
    foci = nw_data[doi]['foci']
    foci = np.concatenate(foci)
    marker_color = col_fMRI[1] if nw_data[doi]['type']=='fMRI' else col_fNIRS[1]
    fig = plt.figure(figsize=(2,2))
    gb = plot_connectome(np.full([len(foci),len(foci),], np.nan), foci, figure=fig, colorbar=None, 
                         display_mode='z', node_size=130, node_color=marker_color, annotate=False) #output_file=brain_path, 
    # fNIRS: all channel coords
    if nw_data[doi]['type']=='fNIRS':
        all_foci = ds_fNIRS.query("doi==@doi")[["x","y","z"]].values
        gb.add_markers(all_foci, alpha=0.4, marker_color="black")
    # save
    plt.savefig(join(sdir, 'citenet_brains', f'doi_{n_doi}.svg'), transparent=True, bbox_inches='tight')
    plt.close()
    nw_data[doi]['brain'] = f'citenet_brains/doi_{n_doi}.svg'


# save
#with open(join(sdir, 'citenet', 'included_studies_fMRI.json'), 'w') as f:
#    json.dump(nw_data, f)

# %% citation network ========================================================================

# load oc data
#with open(join(sdir, 'citenet', 'included_studies_fMRI.json')) as f:
#    nw_data = json.load(f)

# get list of dois
dois = list(nw_data.keys())

# get areas for coloring
#areas = [nw_data[record]['exps'][0]['area'] for record in nw_data]
#communities = list()
#for area in list(set(areas)):
#    communities.append([dois[i] for i,a in enumerate(areas) if a==area])
    
## create Networkx network -----------------------------------------------    
nx = networkx.DiGraph() # initialize network
nx.add_nodes_from(dois) # add nodes, nodes are referred to by doi
    
## draw edges
for record in nw_data:   
    ## get dois of references within our network and add edges       
    refs = [r for r in dois if r in nw_data[record]["reference"]]
    [nx.add_edge(record, ref) for ref in refs]   
    ## do the same but for citations - redundant, necessary?
    cits = [c for c in dois if c in nw_data[record]["citation"]]
    [nx.add_edge(cit, record) for cit in cits]

## create pyvis network --------------------------------------------------  
    
## layout
# network layout is passed to visjs as a string in json format
nw_layout = """
    var options = {
        "nodes": {
            "scaling": {
                "min": 30,
                "max": 50
            },
            "borderWidth": 1.5,
            "borderWidthSelected": 3.5,
            "font": {
                "size": 20,
                "background": "rgba(255,255,255,0.5)"
                },
            "shadow": {
                "enabled": true
                }
            },
        "edges": {
            "font": {
                "color": "transparent",
                "strokeColor": "transparent",
                "background": "transparent"
                },
            "color": {
                "inherit": true,
                "opacity": 0.6
                },
            "smooth": {
                "type": "continuous",
                "forceDirection": "none"
                }
            },
        "interaction": {
            "keyboard": {
                "enabled": true
                },
            "navigationButtons": true,
            "zoomView": true
            },
        "layout": {
            "improvedLayout": false,
            "randomSeed": 4
        },
        "physics": {
            "enabled": true,
            "forceAtlas2Based": {
                "gravitationalConstant": -100,
                "springLength": 220,
                "springConstant": 0.01,
                "avoidOverlap": 0.5
                },
            "maxVelocity": 50,
            "solver": "forceAtlas2Based",
            "timestep": 0.5,
            "stabilization": {
                "enabled": true,
                "fit": true
                }
            }
        }
    """
    
    ## initialize pyvis network and populate with nx data 
nw = Network(directed=True,   
             height='100%', width='100%', bgcolor='#ffffff', font_color=False, 
             heading="")
nw.from_nx(nx) # load nx data
nw.set_options(nw_layout)
    
## iterate over records, find respective node, add metadata & set details
for i, record in enumerate(nw_data):
    # author & year
    authoryear = nw_data[record]['author']+', '+nw_data[record]['year']
    # experiments
    exp_info = ''
    for i, exp in enumerate(nw_data[record]['exps']):
        exp_info += """
            <b>Experiment {i}: {name}</b><br>
            <u>Area</u>: {area}<br> 
            <u>Task</u>: {task}<br>
            <u>Contrast(s)</u>: {contrast}<br>
            <u>Sample</u>: {nsub:.0f} subj; {ndy:.0f} dyads; {age:.2f} y; {fem:.0%} ♀ <br>
            """.format(i=i+1, name=exp['exp'], area=exp['area'], task=exp['task'], 
                       contrast='; '.join(exp['contrast']), nsub=exp['n'], ndy=exp['nDyads'],
                       age=exp['ageMean'], fem=exp['propFem'])
    ## add metadata to respective node
    nw.get_node(record)["label"] = authoryear
    nw.get_node(record)["title"] = """
        <div style="background-image:url('{img}'); background-repeat:no-repeat; background-position:center;
                    background-color:white; font-family:arial; padding:0; margin:0;
                    border-style:solid; border-width:10px; border-color:white"> 
        <div style="background:rgba(255,255,255,0.6); padding:0; margin:0">
            <p style="font-size:15; padding:0; margin:0">
                <b>{a} ({type})</b>: <br> 
                <a target="_blank" href=https://doi.org/{d}><i>{t}</i></a> <br>
            </p>
            <p style="font-size:13; margin:0; padding:0">
                <br> {exps} <br>
                <u>Citation-count</u>: NW: {c_nw}, OC: {c_all} <br>
                <u>Reference-count</u>: NW: {r_nw}, OC: {r_all}
            </p>
        </div>
        </div>
        """.format(a=authoryear, type=nw_data[record]["type"], t=nw_data[record]["title"], d=record, exps=exp_info,
                    c_nw=nx.in_degree[record], r_nw=nx.out_degree[record], 
                    c_all=nw_data[record]["citation_count"], r_all=len(nw_data[record]["reference"]),
                    img=nw_data[record]['brain'])
    
    ## set node size
    nw.get_node(record)["value"] = nx.in_degree[record]
    ## node shape
    nw.get_node(record)["shape"] = "circularImage"
    nw.get_node(record)["image"] = nw_data[record]['brain']
    #nw.get_node(record)["shape"] = "dot" if nw_data[record]["type"] == 'fMRI' else 'diamond'
    ## node color
    nw.get_node(record)["color"] = {
        "border": col_fMRI[0] if nw_data[record]["type"]=='fMRI' else col_fNIRS[0],
        "background": 'white',
        "highlight": {
            "border": col_fMRI[0] if nw_data[record]["type"]=='fMRI' else col_fNIRS[0],
            "background": 'white'
        }
    }

    ## color-code communities
    # iterate over nodes in communites and set colors
    #for i, c in enumerate(communities):
    #    for node in c:
    #        nw.get_node(node)["color"] = {
    #            "border": colors[i],
    #            "background": 'white',
    #            "highlight": {
    #                "border": colors[i],
    #                "background": 'white'
    #            }
    #        }

# legend
nw.add_node('fMRI', label='fMRI', shape="circularImage", value=40, borderWidth=2,
            image='citenet_brains/legend.svg', physics=False)
nw.get_node('fMRI')["color"] = {"border": col_fMRI[0], "background": 'white',
                                "highlight": {"border": col_fMRI[0], "background": 'white'}}
nw.get_node('fMRI')["x"] = -750
nw.get_node('fMRI')["y"] = 800
nw.get_node('fMRI')['font'] = {"size": 30}
nw.add_node('fNIRS', label='fNIRS', shape="circularImage", value=40, borderWidth=2,
            image='citenet_brains/legend.svg', physics=False)
nw.get_node('fNIRS')["color"] = {"border": col_fNIRS[0], "background": 'white',
                                "highlight": {"border": col_fNIRS[0], "background": 'white'}}
nw.get_node('fNIRS')["x"] = -640
nw.get_node('fNIRS')["y"] = 800
nw.get_node('fNIRS')['font'] = {"size": 30}
## save
nw.show('fig2B.html')

# %%
