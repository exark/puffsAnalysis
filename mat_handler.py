from scipy.io import loadmat
import numpy as np
import argparse
import cv2
import numpy as np
import glob
import math
from itertools import product
import re
import h5py

def import_scored_file_nonhdf5(scored_file):
    dat = loadmat(scored_file, struct_as_record=False, squeeze_me=True)
    r = re.compile('.*tracks.*', re.IGNORECASE)
    dat = dat[list(filter(r.match, dat.keys()))[0]]
    new_dat = np.concatenate([np.column_stack((np.repeat(ind, len(frame.f)), np.repeat(frame.isPuff, len(frame.f)),
                                            frame.f, frame.x, frame.y,
                                           np.repeat(1, len(frame.f)))) for ind, frame in enumerate(dat)])
    # uncomment if you only want to extract puffs:
    # new_dat = new_dat[np.where(new_dat[:,1] == 1)]
    return new_dat

def import_scored_file_hdf5(scored_file, with_buffer = False):
    dat = h5py.File(scored_file, 'r')
    r = re.compile('.*tracks.*', re.IGNORECASE)
    dat = dat[list(filter(r.match, dat.keys()))[0]]
    # store data with one row for each puff-frame
    # the values in each row are a puff id, isPuff, frame number, x location, y location, and an
    # indicator for whether the puff-frame was one originally identified (not one of the additional frames
    # we later add)
    new_dat = np.concatenate([np.column_stack((np.repeat(ind, dat[dat['y'][ind,0]].shape[0]),
                                               np.repeat(dat[dat['isPuff'][ind,0]][0,0],
                                                         dat[dat['y'][ind,0]].shape[0]),
                                               dat[dat['f'][ind,0]][:,0],
                                               dat[dat['x'][ind,0]][:,0],
                                               dat[dat['y'][ind,0]][:,0],
                                               dat[dat['c'][ind,0]][:,0],
                                               dat[dat['A'][ind,0]][:,0],
                                               np.repeat(1, dat[dat['y'][ind,0]].shape[0])))
                         for ind in range(dat['y'].shape[0])])
    # uncomment if you only want to extract puffs:
    # new_dat = new_dat[np.where(new_dat[:,1] == 1)]
    if with_buffer:
        buffers = np.concatenate([np.column_stack((
                np.repeat(ind,3),
                np.array([0, 1, 2]),
                (np.repeat(0,3) if isinstance(dat[dat['startBuffer'][ind,0]], h5py.Dataset) else dat[dat['startBuffer'][ind,0]]['A'][:,0]),
                (np.repeat(0,3) if isinstance(dat[dat['endBuffer'][ind,0]], h5py.Dataset) else dat[dat['endBuffer'][ind,0]]['A'][:,0])))
            for ind in range(dat['y'].shape[0])])
        return new_dat, buffers
    else:
        return new_dat

# takes a matlab struct containing the scored info and converts it
# to a numpy array with one row per puff-frame
def import_scored_file(scored_file):
    try:
        dat = import_scored_file_nonhdf5(scored_file)
    except:
        print("HDF5 file, retrying...")
        dat = import_scored_file_hdf5(scored_file)
    return dat

def load_mat_files(name, condition, fields, split = 'within_cell'):
    # split "within_cell" will split each cell in half, to training and heldout
    # split "in_half" will evenly split the total tracks in a condition across training and heldout
    # split "intelligently" will order cells by mean background and then split cells down the list between training and heldout
    shape = (0,len(fields))
    cargo = {'data' : pd.DataFrame(columns=fields),
             'training_data' : pd.DataFrame(columns=fields),
             'heldout_data' : pd.DataFrame(columns=fields)}

    data_list = [];
    training_list = [];
    heldout_list = [];

    bar = progressbar.ProgressBar()
    for file in bar(condition['files']):
        data = mat2py(condition['dir'] + file + '.mat', fields, condition['dir'], file)
        data = pd.DataFrame(np.nan_to_num(np.array(data.tolist())),
                            columns=fields,
                            index=pd.Index(np.r_[1:np.shape(data)[0]+1],name='track'))
        if split == 'within_cell':
            puffs = data[data.isPuff==1.]
            nonpuffs = data[data.isPuff==2.]
            r = np.random.RandomState(327)
            puff_mask = r.choice([True, False], puffs.shape[0], p=[0.5, 0.5])
            nonpuff_mask = r.choice([True, False], nonpuffs.shape[0], p=[0.5, 0.5])
            training_data = pd.concat([puffs[puff_mask], nonpuffs[nonpuff_mask]])
            heldout_data = pd.concat([puffs[[not x for x in puff_mask]], nonpuffs[[not x for x in nonpuff_mask]]])
            data_list.append(data)
            training_list.append(training_data)
            heldout_list.append(heldout_data)
        else:
            data_list.append(data)

    if split == 'within_cell':
        cargo['data'] = pd.concat(data_list, keys=condition['files'], names=['cell'])
        cargo['training_data'] = pd.concat(training_list, keys=condition['files'], names=['cell'])
        cargo['heldout_data'] = pd.concat(heldout_list, keys=condition['files'], names=['cell'])
    elif split == 'in_half':
        data = pd.concat(data_list, keys=condition['files'], names=['cell'])
        half_frame = int(np.shape(data)[0]/2)
        training_data = data.iloc[0:half_frame,:]
        heldout_data = data.iloc[half_frame+1:,:]
        cargo['data'] = data
        cargo['training_data'] = training_data
        cargo['heldout_data'] = heldout_data
    elif split == 'intelligently':
        data = pd.concat(data_list, keys=condition['files'], names=['cell'])
        global_background = {};
        for cell in np.unique(data.index.get_level_values('cell')):
            cell_data = data.xs(cell, level=0)
            global_background[cell] = np.mean(cell_data['global_background'])
        sorted_by_background = [data.xs(k, level=0) for k in sorted(global_background,
                                                                    key=global_background.get)]
        cargo['data'] = data
        cargo['training_data'] = pd.concat(sorted_by_background[::2], keys=condition['files'], names=['cell'])
        cargo['heldout_data'] = pd.concat(sorted_by_background[1::2], keys=condition['files'], names=['cell'])
    return cargo
