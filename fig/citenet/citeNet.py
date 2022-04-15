#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 06:24:20 2021

@author: Leon D. Lotter
"""

# ========================================================================================


def get_records_from_dois(doi_list, in_chunks=10, max=None):
    
    import sys
    import re    
    import requests

    print("""Warning: This will request all metadata from OpenCitations! 
    It is likely that some records will be missing. Consider exporting your 
    data to bibtex and using 'get_records_from_bibtex' which will query 
    OpenCitations only for reference and citation lists.
    """)

    ## check input type
    if not isinstance(doi_list, list):
        sys.exit("Wrong input type, has to be list!") 

    # cut list to max
    if max:
        doi_list = doi_list[:max]

    ## query function
    # OpenCitations API
    base_url = "https://opencitations.net/index/coci/api/v1/metadata/"
    # query
    def query_oc(dois):
        response = requests.get(base_url+"__".join(dois))
        response_list = response.json()
        oc_dat = {}
        for response in response_list:
            oc_dat[response.get("doi")] = {
                "author": response.get("author").split("; "),
                "year": response.get("year"),
                "title": response.get("title"),
                "reference": response.get("reference").split("; "),
                "citation": response.get("citation").split("; "),
                "citation_count": response.get("citation_count")
            }
        return(oc_dat)

    ## clean doi list
    dois = clean_doi_list(doi_list)
    
    ## get OpenCitations metadata 
    print(f"Downloading OpenCitations data for {len(dois)} records, could take a while...") 
    if in_chunks is not None:
        chunks = [dois[i:i+in_chunks] for i in range(0, len(dois), in_chunks)]
        oc_dat = {}
        for i, d in enumerate(chunks):
            print("Downloading chunk {}/{}: {} doi(s) ({})".format(
                i+1, len(chunks), len(d), "; ".join(d)))
            oc_dat.update(query_oc(d))
    else:
        oc_dat = query_oc(dois)
        
    oc_dois = list(oc_dat.keys()) # get OpenCitation dois
    
    # if a doi is not found, it will be ignored
    if len(dois) != len(oc_dois):
        diff = list(set(dois) - set(oc_dois))
        print(f"Did not find the following {len(diff)} record(s):", *diff, sep = "\n")
    else:
        print("All records found.")
    
    ## return
    print("Return OpenCitations data.\n")
    return(oc_dat)


# ========================================================================================


def clean_doi_list(doi_list):
    import re

    ## clean doi list
    # use regex to remove whitespaces, "doi.org", etc
    dois = [re.sub("\ |https://|www.|doi.org/", "", d) for d in doi_list]
    # remove potential duplicates
    dois = list(set(dois))
    # make lower-case
    dois = [d.lower() for d in dois]

    return(dois)


# ========================================================================================


def get_references_from_doi(doi, which_list="references"):

    import requests

    # OpenCitations API
    base_url = "https://opencitations.net/index/coci/api/v1/"
    # query
    which_field = "cited" if which_list=="references" else "citing"
    response = requests.get(base_url+which_list+'/'+doi)
    try:
        response_list = response.json()
        ref_list = list()
        for response in response_list:
            ref_list.append(response.get(which_field))
    except:
        print(f"Problem with doi: {doi}. Proceed.")
        ref_list = []
        pass

    return(ref_list)


# ========================================================================================


def get_dois_from_bibtex(bib_file):
    
    from pybtex.database import parse_file
    from crossref.restful import Works
    works = Works()
    
    ## load bibtex
    print("Importing from bibtex.")
    bib = parse_file(bib_file)
    
    ## iterate over entries and get dois
    doi_list = []
    for entry in bib.entries:
        record = bib.entries[entry]
        
        ## get doi
        # search crossref if doi not in bibtex
        if "doi" not in record.fields or record.fields["doi"] == "":
            # get title, year and first author
            titleyear = record.fields["title"] + ", " + record.fields["year"]
            author = str(record.persons["author"][0])
            print(f"""No DOI for {author}, {record.fields["year"]}. Retrieving from Crossref.""")
            # search crossref
            cf_data = works.query(bibliographic=titleyear, author=author)
            # assumes first result to be correct!
            doi = [r["DOI"] for r in cf_data][0]
            print(f"Adding DOI: {doi}")
        # if doi in bibtex, use it
        else:
            doi = record.fields["doi"]
        # create list
        doi_list.append(doi)

    doi_list = clean_doi_list(doi_list)

    print(f"Extracted {len(doi_list)} dois from Bibtex.\n")

    return(doi_list)


# ========================================================================================


def get_records_from_bibtex(bib_file):

    from pybtex.database import parse_file
    from crossref.restful import Works
    import re
    works = Works()
    
    ## load bibtex
    print("Importing from bibtex.")
    bib = parse_file(bib_file)
    
    ## iterate over entries and get dois
    oc_data = {}
    for i, entry in enumerate(bib.entries):
        record = bib.entries[entry]
        print(f'Looking into entry {i+1}/{len(bib.entries)}: {record.fields["title"]}')

        # ------------------------------------------------------------------
        ## get doi
        # search crossref if doi not in bibtex
        if "doi" not in record.fields or record.fields["doi"] == "":
            # get title, year and first author
            titleyear = record.fields["title"] + ", " + record.fields["year"]
            author = str(record.persons["author"][0])
            print(f"""No DOI for {author}, {record.fields["year"]}. Retrieving from Crossref.""")
            # search crossref
            cf_data = works.query(bibliographic=titleyear, author=author)
            # assumes first result to be correct!
            doi = [r["DOI"] for r in cf_data][0]
            print(f"Adding DOI: {doi}")
        # if doi in bibtex, use it
        else:
            doi = record.fields["doi"]
        # clean
        doi = clean_doi_list([doi])
        doi = doi[0]

        # ------------------------------------------------------------------
        # get data
        references = get_references_from_doi(doi, which_list="references"),
        citations = get_references_from_doi(doi, which_list="citations")
        rm = r'[{}.]'
        author = [re.sub(rm, '', str(a)) for a in record.persons["author"]]
        title = re.sub(rm, '', record.fields["title"])
        oc_data[doi] = {
            "author": author,
            "year": record.fields["year"],
            "title": title,
            "reference": [r for r in references[0]],
            "citation": citations, 
            "citation_count": len(citations)
        }
        
    print(f"Bibtex with {len(oc_data)} records imported.\n")

    return(oc_data)


# ========================================================================================


def citenet(oc_dat, shapes=None, h="95%", w="100%", node_colors=None, border_colors=None, 
            bg_color="#ffffff", font_color=False, 
            font_size=17, head="", physics=True, gravitation=None, spring_len=100, load_stable=True, 
            communities=True,seed=None, nav_buttons=True, zoom_mousewheel=True): 
 
    import networkx
    from pyvis.network import Network
    from scipy.stats import percentileofscore

    # function to calculate percentiles
    def perc(dictionary, doi):
        score = dictionary.pop(doi)
        p = percentileofscore(list(dictionary.values()), score, "rank")
        return(round(p, 1))

    # function to get abbreviated author names
    def get_authoryear(author_data, year):
        if len(author_data) == 1:
            authoryear = author_data[0].split(",")[0]+", "+year
        if len(author_data) == 2:
            authoryear = author_data[0].split(",")[0]+" & "+author_data[1].split(",")[0]+", "+year
        if len(author_data) > 2:
            authoryear = author_data[0].split(",")[0]+" et al., "+year
        return(authoryear)

    # get list of dois
    oc_dois = list(oc_dat.keys())

    # check shapes
    if shapes:
        if (len(shapes) != len(oc_dois)):
            print("""Shapes must be a list of shape strings with the same length as oc_dat! 
            Won't use it now. Currently available: 
            'diamond', 'dot', 'star', 'triangle', 'triangleDown', 'hexagon', 'square'""")
            shapes = None
    
    ## create Networkx network -----------------------------------------------
    print(f"Creating network from {len(oc_dois)} records.")
    
    nx = networkx.DiGraph() # initialize network
    nx.add_nodes_from(oc_dois) # add nodes, nodes are referred to by doi
    
    ## draw edges
    for record in oc_dat:   
        ## get record doi & author(s)
        #authoryear = get_authoryear(oc_dat[record]["author"], oc_dat[record]["year"])
        
        ## get dois of references within our network and add edges       
        refs = [r for r in oc_dois if r in oc_dat[record]["reference"]]
        [nx.add_edge(record, ref) for ref in refs]
        
        ## do the same but for citations - redundant, necessary?
        #cits = [c for c in oc_dois if c in record.citation]
        #[nw.add_edge(cit, doi) for cit in cits]
    
    ## compute centrality estimate       
    nx_undir = nx.to_undirected() # convert to undirected nw
    cent = networkx.degree_centrality(nx_undir)
    cent_eig = networkx.eigenvector_centrality(nx_undir)
    cent_betw = networkx.betweenness_centrality(nx_undir)
    cent_close = networkx.closeness_centrality(nx_undir)

    if communities:
        ## get communities using greedy modularity maximization
        comm = networkx.algorithms.community.greedy_modularity_communities(nx_undir) # get communities
        comm = [sorted(comm[i]) for i in range(len(comm))] # make list of lists
    
    ## create pyvis network --------------------------------------------------  
    
    ## layout
    # network layout is passed to visjs as a string in json format
    # "transform" input to match with expected visjs input
    buttons = "true" if nav_buttons == True else "false"
    zoom = "true" if zoom_mousewheel == True else "false"
    physics = "true" if physics == True else "false"
    stab = "true" if load_stable == True else "false"
    grav = -4 * nx.number_of_nodes() if gravitation is None else gravitation
    #add_seed = """ "layout": {"randomSeed": %d}, """ % seed if seed is not None else ""
    # define layout
    nw_layout = """
        var options = {
            "nodes": {
                "borderWidth": 2,
                "borderWidthSelected": 4,
                "font": {
                    "size": %d,
                    "background": "rgba(255,255,255,0.5)"
                    },
                "shadow": {
                    "enabled": true
                    }
                },
            "edges": {
                "color": {
                    "inherit": true
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
                "navigationButtons": %s,
                "zoomView": %s
                },
            "physics": {
                "enabled": %s,
                "forceAtlas2Based": {
                    "gravitationalConstant": %d,
                    "springLength": %d,
                    "avoidOverlap": 1
                    },
                "minVelocity": 0.75,
                "solver": "forceAtlas2Based",
                "stabilization": {
                    "enabled": %s,
                    "fit": true
                    }
                }
            }
        """ % (font_size, buttons, zoom, physics, grav, spring_len, stab)
    
    ## initialize pyvis network and populate with nx data 
    nw = Network(directed=True, height=h, width=w, 
                 bgcolor=bg_color, font_color=font_color, heading=head)
    nw.from_nx(nx) # load nx data
    nw.set_options(nw_layout)
    
    # set css style for hover text fields with record info
    hover_style = "background-color:white; font-family:arial; margin:0"
    # colorblind-friendly colors to automatically color communities (9 max)
    colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', 
              '#984ea3', '#999999', '#e41a1c', '#dede00']
    colors_light = ["#c0d8ec", "#ffd8b2", "#c8e7c7", "#fcd9eb", "#edc9b5",
                    "#e1c8e4", "#e0e0e0", "#f7baba", "ffffa8"]

    ## iterate over records, find respective node, add metadata & set details
    for i, record in enumerate(oc_dat):
        
        ## get record doi & author(s)
        authoryear = get_authoryear(oc_dat[record]["author"], oc_dat[record]["year"])

        ## add metadata to respective node
        nw.get_node(record)["label"] = authoryear
        nw.get_node(record)["title"] = """
            <p style="{s}"> 
                <b>{a}</b>: <br> 
                <i>{t}</i> <br> 
                <a target="_blank" href=https://doi.org/{d}> {d} </a> <br>
                NW citations: {c_nw} <br>
                NW references: {r_nw} <br>
                OC citations: {c_all} <br>
                OC references: {r_all} <br>
                Degree: {deg} ({deg_p}) <br>
                Eigenvector: {eig} ({eig_p}) <br>
                Betweenness: {bet} ({bet_p}) <br>
                Closeness: {clo} ({clo_p})
            </p>
            """.format(s=hover_style, a=authoryear, t=oc_dat[record]["title"], d=record, 
                       c_nw=nx.in_degree[record], r_nw=nx.out_degree[record], 
                       c_all=oc_dat[record]["citation_count"], r_all=len(oc_dat[record]["reference"]),
                       deg=round(cent[record], 2), deg_p=perc(cent, record),
                       eig=round(cent_eig[record], 2), eig_p=perc(cent_eig, record),
                       bet=round(cent_betw[record], 2), bet_p=perc(cent_betw, record),
                       clo=round(cent_close[record], 2), clo_p=perc(cent_close, record),)
        
        ## set node size
        nw.get_node(record)["value"] = nx.in_degree[record]
        ## set node shape
        if shapes:
            nw.get_node(record)["shape"] = shapes[i]
        ## set node color
        if node_colors and border_colors:
            nw.get_node(node)["color"] = {
                "border": border_colors[i],
                "background": node_colors[i],
                "highlight": {
                    "border": border_colors[i],
                    "background": node_colors[i]
                }  
            } 
        
    ## color-code communities
    # iterate over nodes in communites and set colors
    if communities:
        for i, c in enumerate(comm):
            for node in c:
                nw.get_node(node)["color"] = {
                    "border": colors[i],
                    "background": colors_light[i],
                    "highlight": {
                        "border": colors[i],
                        "background": colors_light[i]
                    }
                }
    
    ## return network
    if communities:
        print(f"Return network; {nx.number_of_nodes()} nodes, {nx.number_of_edges()} edges, {len(comm)} communities.")
    else:
        print(f"Return network; {nx.number_of_nodes()} nodes, {nx.number_of_edges()} edges.")
    return(nw, nx)    

