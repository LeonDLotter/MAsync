#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:36:41 2022

@author: Leon D. Lotter
"""

from nimare.meta import ALE, MKDAChi2
from glob import glob
from datetime import datetime
from os.path import join, exists
from pprint import pprint
from nimare.extract import download_abstracts, fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare.dataset import Dataset

import logging
lgr = logging.getLogger(__name__)
lgr.setLevel(logging.INFO)


def download_neurosynth(save_path, version=7, abstract=True, vocab='terms', 
                        email=None, overwrite=True):
    """
    Downloads neurosynth database version {version} and abstracts to 
    {save_path}, converts to NiMARE database and returns database
    """
    

    
    ## Download
    lgr.info('Downloading Neurosynth database...')
    files = fetch_neurosynth(
        data_dir=save_path,
        version=f'{version}',
        overwrite=overwrite,
        source='abstract',
        vocab=vocab,
    )
    # files are saved to a new folder within "save_path" named "neurosynth".
    pprint(files)
    neurosynth_db = files[0]
    
    ## Convert to NiMARE dataset
    dset_file = join(save_path, 'neurosynth_dataset.pkl.gz')
    if not exists(dset_file) or overwrite == True:
        lgr.info('Converting to NiMARE dataset...')
        neurosynth_dset = convert_neurosynth_to_dataset(
            coordinates_file=neurosynth_db["coordinates"],
            metadata_file=neurosynth_db["metadata"],
            annotations_files=neurosynth_db["features"]
        )
        neurosynth_dset.save(dset_file)
    else:
        lgr.info(f'Loading existing NiMARE dataset from {dset_file}')
        neurosynth_dset = Dataset.load(dset_file)    
    print(neurosynth_dset)
    
    ## Download abstracts
    if abstract == True:
        dset_abstract_file = join(save_path, 'neurosynth_dataset_with_abstracts.pkl.gz')
        if not exists(dset_abstract_file) or overwrite == True:      
            lgr.info('Downloading abstracts...')
            if email is not None:
                neurosynth_dset = download_abstracts(neurosynth_dset, email)
                neurosynth_dset.save(join(save_path, 
                                          "neurosynth_dataset_with_abstracts.pkl.gz"))
            elif email is None:
                lgr.error('Provide email adress to download abstracts.')
        else:
            lgr.info('Loading NiMARE dataset with abstracts from '
                     f'{dset_abstract_file}')
            neurosynth_dset = Dataset.load(dset_abstract_file)    
            neurosynth_dset = Dataset.load(dset_file)    
    
    ## Return
    lgr.info('Finished. Returning NiMARE dataset')
    return(neurosynth_dset)


#=============================================================================


def create_neurosynth_topic_maps(ds, save_path, estimator='mkdachi2', topics=None, 
                                 topic_pref='LDA200_abstract_weight__', 
                                 topic_thresh=0.001, sample_size=10, 
                                 save_prefix='neurosynth_topic_',
                                 maps=['z'], overwrite=False):

    # get topics
    if topics is None:
        topics = ds.get_labels()

    lgr.info(f'Using dataset with {len(ds.ids)} studies and {len(topics)} topics.')
    if sample_size is not None:
        # overwrite sample size in ds
        ds.metadata = ds.metadata.assign(sample_sizes=sample_size)
        lgr.info(f'Setting sample size to constant n={sample_size}.')

    # get estimator
    if estimator == 'ale':
        meta = ALE()
    elif estimator == 'mkdachi2':
        meta = MKDAChi2(kernel__r=10)
    else: 
        lgr.error(f"'estimator' must be 'ale' or 'mkdachi2' not '{estimator}'!")

    # loop over topics
    for topic in topics:

        # topic name and map path
        #topic_name = topic.split("__", 1)[1]
        topic_maps = glob(join(save_path, f'{save_prefix}{topic}*.nii.gz'))

        # check if exist
        if overwrite == True or len(topic_maps) == 0:

            # get studies with {topic} > {topic_thresh}
            ids_topic = ds.get_studies_by_label(labels=topic_pref+topic, label_threshold=topic_thresh)
            # slice new dataset
            ds_topic = ds.slice(ids_topic)
            lgr.info(f'Topic {topic}: {len(ds_topic.ids)}/{len(ds.ids)} studies selected. '
                     f'Creating {estimator} map...')

            # create meta map
            now = datetime.now()
            if estimator == 'ale':
                meta_topic = meta.fit(ds_topic)
            else:
                ids_NOT_topic = list(set(ds.ids).difference(ids_topic))
                ds_NOT_topic = ds.slice(ids_NOT_topic)
                meta_topic = meta.fit(ds_topic, ds_NOT_topic)

            # save
            save_name = f'{save_prefix}{topic}_{len(ds_topic.ids)}'
            meta_topic.save_maps(output_dir=save_path, prefix=save_name, names=maps)
            lgr.info(f'Duration: {datetime.now() - now}. Saved {estimator} map(s) to '
                    f'{save_path}/{save_name}_*.nii.gz')

        else:
            lgr.info(f'Topic {topic}: file exists and overwrite is True -> skipping file.')


    