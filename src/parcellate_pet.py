
# %%
from nilearn.input_data import NiftiLabelsMasker
from scipy.stats import zscore 
from os.path import join
import pandas as pd

pet_dir = '/Users/leonlotter/MAsync/PET/maps'
atlas = '/Users/leonlotter/MAsync/project/data/atlases/Schaefer100-7_TianS1_2mm.nii.gz'
parcellater = NiftiLabelsMasker(atlas)

def get_data(maps, n_subj=None):
    dat_list = list()
    for map in maps:
        dat = parcellater.fit_transform(map)
        dat = zscore(dat, axis=1)
        dat_list.append(dat)
    if len(dat_list) > 1:
        dat_weighted = [dat * n_subj[i] for i, dat in enumerate(dat_list)]
        dat_weighted = sum(dat_weighted) / sum(n_subj)
        return(dat_weighted.tolist()[0])
    else:
        return(dat_list[0].tolist()[0])
        

# %%
pet_dict = {
    '5HT1a': get_data([join(pet_dir, '5HT1a-way100635-36-savli2012.nii.gz')]),
    '5HT1b': get_data([join(pet_dir, '5HT1b-p943-23-savli2012.nii.gz'), 
                       join(pet_dir, '5HT1b-p943-65-gallezot2010.nii.gz')], 
                      [23, 65]),
    '5HT2a': get_data([join(pet_dir, '5HT2a-cimbi36-29-beliveau2017.nii.gz')]),
    '5HT4': get_data([join(pet_dir, '5HT4-sb207145-59-beliveau2017.nii.gz')]),
    '5HT6': get_data([join(pet_dir, '5HT6-gsk215083-30-radhakrishnan2018.nii.gz')]),
    '5HTT': get_data([join(pet_dir, '5HTT-dasb-100-beliveau2017.nii.gz'), 
                      join(pet_dir, '5HTT-dasb-18-savli2012.nii.gz')], 
                      [100, 18]),
    'a4b2': get_data([join(pet_dir, 'A4B2-flubatine-30-hillmer2016.nii.gz')]),
    'CB1': get_data([join(pet_dir, 'CB1-omar-77-normandin2015.nii.gz')]),
    'D1': get_data([join(pet_dir, 'D1-sch23390-13-kaller2017.nii.gz')]),
    'D2': get_data([join(pet_dir, 'D2-flb457-37-smith2019.nii.gz'), 
                    join(pet_dir, 'D2-flb457-55-sandiego2015.nii.gz')], 
                    [37, 55]),
    'DAT': get_data([join(pet_dir, 'DAT-fpcit-174-dukart2018.nii.gz'),
                     join(pet_dir, 'DAT-fpcit-30-garciagomez2013.nii.gz')],
                     [174, 30]),
    'FDOPA': get_data([join(pet_dir, 'FDOPA-fluorodopa-12-garciagomez2018.nii.gz')]),
    'GABAa': get_data([join(pet_dir, 'GABAa-flumazenil-16-norgaard2020.nii.gz'), 
                       join(pet_dir, 'GABAa-flumazenil-6-dukart2018.nii.gz')], 
                      [16, 6]),
    'H3': get_data([join(pet_dir, 'H3-gsk189254-8-gallezot2017.nii.gz')]),
    'M1': get_data([join(pet_dir, 'M1-lsn3172176-24-naganawa2021.nii.gz')]),
    'mGluR5': get_data([join(pet_dir, 'mGluR5-abp688-22-rosaneto.nii.gz'), 
                        join(pet_dir, 'mGluR5-abp688-28-dubois2015.nii.gz'),
                        join(pet_dir, 'mGluR5-abp688-73-smart2019.nii.gz')], 
                       [22, 28, 73]),
    'MU': get_data([join(pet_dir, 'MU-carfentanil-204-kantonen2020.nii.gz'), 
                    join(pet_dir, 'MU-carfentanil-39-turtonen2021.nii.gz')], 
                   [204, 39]),
    'NET': get_data([join(pet_dir, 'NET-mrb-10-hesse2017.nii.gz'), 
                    join(pet_dir, 'NET-mrb-77-ding2010.nii.gz')], 
                   [10, 77]),
    'NMDA': get_data([join(pet_dir, 'NMDA-ge179-29-galovic2021.nii.gz')]),
    'VAChT': get_data([join(pet_dir, 'VAChT-feobv-18-aghourian2017.nii.gz'), 
                       join(pet_dir, 'VAChT-feobv-4-tuominen.nii.gz'),
                       join(pet_dir, 'VAChT-feobv-5-bedard2019.nii.gz')], 
                      [18, 4, 5])
    }

#%%
pet_data = pd.DataFrame.from_dict(pet_dict)
pet_data.to_csv('/Users/leonlotter/MAsync/project/data/datasets/pet_parcellated_data.csv', index=False)

# %%
