
# %% Import

from glob import glob
from os.path import join, split, splitext
from nilearn.image import load_img, resample_img
import numpy as np
from datetime import datetime

# %% Get data
hcp_dir = '/Volumes/Data_4TB/RSFC/HCP_data'
hcp_rest = glob(join(hcp_dir, '*/**/rfMRI_REST1_*_hp2000_clean.nii.gz'), 
                recursive = True)
hcp_struct = glob(join(hcp_dir, '*/**/T1w_restore_brain.nii.gz'), 
                recursive = True)

# %% Iterate and resample rest data
starttime = datetime.now()
for i, rest in enumerate(hcp_rest):
    # file parts
    pathfile = split(rest)
    dir_old = pathfile[0]
    file_old = pathfile[1]
    file_new = file_old.split('.nii.gz')[0] + '_3mm.nii'
    path_new = join(dir_old, file_new)

    # resample
    print(f'Rest file {i}: Resampling {rest} to {path_new}...')
    now = datetime.now()
    res = resample_img(rest, target_affine=np.diag([3,3,3]), interpolation='linear')
    res.to_filename(path_new)
    print(f'Duration: {datetime.now() - now}.')
print(f'Finished rest data, total duration: {datetime.now() - starttime}.')

# %% Iterate and unpack structural data
starttime = datetime.now()
for i, struct in enumerate(hcp_struct):
    # file parts
    pathfile = split(struct)
    dir_old = pathfile[0]
    file_old = pathfile[1]
    file_new = file_old.split('.gz')[0]
    path_new = join(dir_old, file_new)

    # unpack
    print(f'Structural file {i}: Unzipping {struct} to {path_new}...')
    now = datetime.now()
    struct_vol = load_img(struct)
    struct_vol.to_filename(path_new)
    print(f'Duration: {datetime.now() - now}.')
print(f'Finished structural data, total duration: {datetime.now() - starttime}.')

