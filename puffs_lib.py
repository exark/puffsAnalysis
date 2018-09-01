from pathlib import Path
from scipy.io import loadmat
import numpy as np
import pandas as pd
import argparse
import cv2
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

def import_scored_file_hdf5(scored_file):
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


class PuffMovie:
    def __init__(self, track_file, frame_dir):
        self.track_file = track_file
        self.tracks = None

        frames = Path(frame_dir)
        filelist = list(frames.glob('*.tif*'))
        self.frames = sorted(filelist, key=lambda x: int(x.stem))

    def get_tracks(self):
        if self.tracks is None:
            tracks = import_scored_file(self.track_file)
            tracks = np.nan_to_num(tracks)
            track_columns = ['track','isPuff','f','x','y','c','A','original']
            tracks = pd.DataFrame(tracks, columns = track_columns)
            for column in ['track', 'isPuff', 'f', 'x', 'original']:
                tracks[column] = tracks[column].astype(int)
            self.tracks = pd.DataFrame(tracks, columns = track_columns)

        return self.tracks

    def get_event(self, track_num, window = 4, buffer = 3):
        tracks = self.get_tracks()
        track = tracks.loc[tracks['track'] == track_num, :]
        nframes = np.unique(track['f']).shape[0] + 2*buffer - 1
        cur_frame = iter(range(nframes))
        frames = np.zeros((2*window+1, 2*window+1, nframes))
        start_frame = np.min(track['f']) - buffer - 1
        end_frame = np.max(track['f']) + buffer - 1

        # read in front buffer frames
        row_start = int(np.round(track.loc[track.index[0], 'y']))
        col_start = int(np.round(track.loc[track.index[0], 'x']))
        for f in range(start_frame, np.min(track['f'])-1):
            img = cv2.imread(str(self.frames[f]), -1)
            img = np.nan_to_num(img)
            frame = img[(row_start-window-1):(row_start+window),
                    (col_start-window-1):(col_start+window)]
            frame = np.nan_to_num(frame)
            frames[:,:,next(cur_frame)] = frame

        # read in actual frames
        for f in track.itertuples():
            img = cv2.imread(str(self.frames[f.f-1]), cv2.IMREAD_UNCHANGED)
            img = np.nan_to_num(img)
            row_start = int(np.round(f.y))
            col_start = int(np.round(f.x))
            frame = img[(row_start-window-1):(row_start+window),
                    (col_start-window-1):(col_start+window)]
            frame = np.nan_to_num(frame)
            frames[:,:,next(cur_frame)] = frame

        # read in end buffer frames
        row_start = int(np.round(track.loc[track.index[-1], 'y']))
        col_start = int(np.round(track.loc[track.index[-1], 'x']))
        for f in range(np.max(track['f']), end_frame):
            img = cv2.imread(str(self.frames[f]), cv2.IMREAD_UNCHANGED)
            img = np.nan_to_num(img)
            frame = img[(row_start-window-1):(row_start+window),
                    (col_start-window-1):(col_start+window)]
            frame = np.nan_to_num(frame)
            frames[:,:,next(cur_frame)] = frame

        return puff_event(track, frames)

class puff_event():
    def __init__(self, track_info, frames):
        self.track_info = track_info
        self.frames = frames
