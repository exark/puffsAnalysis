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

def import_scored_file_hdf5(scored_file):
    dat = h5py.File(scored_file, 'r')
    r = re.compile('.*tracks.*', re.IGNORECASE)
    dat = dat[list(filter(r.match, dat.keys()))[0]]
    # store data with one row for each puff-frame
    # the values in each row are a puff id, isPuff, frame number, x location, y location, and an
    # indicator for whether the puff-frame was one originally identified (not one of the additional frames
    # we later add)
    new_dat = np.concatenate([np.column_stack((np.repeat(ind,
                                                     dat[dat['y'][ind,0]].shape[0]),
                                          np.repeat(dat[dat['isPuff'][ind,0]][0,0],
                                                   dat[dat['y'][ind,0]].shape[0]),
                                          dat[dat['f'][ind,0]][:,0],
                                          dat[dat['x'][ind,0]][:,0],
                                          dat[dat['y'][ind,0]][:,0],
                                          np.repeat(1,
                                                   dat[dat['y'][ind,0]].shape[0])))
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


# We add some extra frames to the beginning and the end of each event
# remove the call to this function from the extract_data function if you don't want extra frames
def add_extra_frames(scored_data):
    final_frame = max(scored_data[:,2])
    new_dat = scored_data.copy()
    for ind in np.unique(scored_data[:,0]):
        minframe = np.nanmin(scored_data[np.where(scored_data[:,0] == ind), 2])
        maxframe = np.nanmax(scored_data[np.where(scored_data[:,0] == ind), 2])
        puffval = scored_data[np.where(scored_data[:,0] == ind), 1][0,0]
        xvalmin = scored_data[np.where(scored_data[:,0] == ind), 3][0, 0]
        xvalmax = scored_data[np.where(scored_data[:,0] == ind), 3][0, -1]
        yvalmin = scored_data[np.where(scored_data[:,0] == ind), 4][0, 0]
        yvalmax = scored_data[np.where(scored_data[:,0] == ind), 4][0, -1]
        for j in range(5):
            new_dat = np.concatenate((new_dat, [[ind, puffval, minframe - j - 1, xvalmin, yvalmin, 0]]))
            new_dat = np.concatenate((new_dat, [[ind, puffval, maxframe + j + 1, xvalmax, yvalmax, 0]]))

    # this also serves to remove nan values
    new_dat = new_dat[np.where(new_dat[:,2] > 0)]
    new_dat = new_dat[np.where(new_dat[:,2] <= final_frame)]
    return new_dat

# link the rows in the numpy array of scored data to the frames in the video and extract intensity
# for each pixel in a 9x9 grid at each puff-frame.
# Returns a numpy array with 81 rows for each puff frame (one per pixel in the 9x9 grid) with each row
# containing values for puff id, isPuff, frame number, whether frame is original or one we added, relative
# x and y locations for the pixels, and the intensity in the pixel
def extract_events(videos, scored_data):
    filelist = glob.glob(videos)
    filelist.sort()

    frame_iter = enumerate(filelist)
    file_frame = -10
    sorted_dat = scored_data[np.argsort(scored_data[:,2])]
    res = [[0, 0, 0, 0, 0, 0, 0]]
    for eventid, isPuff, frameno, xloc, yloc, isOrig in sorted_dat:
        while (file_frame+1) < frameno:
            file_frame, fname = next(frame_iter)
            img = cv2.imread(fname, cv2.IMREAD_GRAYSCALE)

        img = np.nan_to_num(img)
        row_start = int(math.floor(yloc))
        col_start = int(math.floor(xloc))
        delta = 4
        block = img[(row_start-delta-1):(row_start+delta),
                (col_start-delta-1):(col_start+delta)]
        block = np.nan_to_num(block)

        side = 2*delta + 1
        row_rem = yloc - row_start
        col_rem = xloc - col_start

        tmp = [[eventid, isPuff, frameno, isOrig, c-delta-col_rem+0.5, r-delta-row_rem+0.5, block[r,c]] for
                r,c in product(range(side),repeat=2)]
        tmp = np.nan_to_num(tmp)
        res = np.concatenate((res, tmp))

    return res

# combines functions into one
def extract_data(scored_file, videos, outfile):
    dat = import_scored_file(scored_file)
    print('Imported data!')
    dat = add_extra_frames(dat)
    events = extract_events(videos, dat)
    print('Got events!')
    np.save(outfile, events)
    print('Saved file succcessfully!')

def main():
    scorefile_list = [#'/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRRachelCell8.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRRachelCell10.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRRachelCell14.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRRachelCell17.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRZYWCell2.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRZYWCell6.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/scoredTracksTfRZYWCell13.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/scoredTracksMORCell1.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/scoredTracksMORCell2.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/scoredTracksMORCell5.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/scoredTracksMORCell8.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/scoredTracksMORCell18.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/scoredTracksMORCell22.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/scoredTracksB2Cell11.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/scoredTracksB2Cell12.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/scoredTracksB2Cell14.mat',
                        '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/scoredTracksB2Cell3.mat']
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/scoredTracksB2Cell4.mat',
                        # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/scoredTracksB2Cell26.mat']
    video_list = [#'/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/No Treatment/Splits/Cell8_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/No Treatment/Splits/Cell10_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/No Treatment/Splits/Cell14_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/No Treatment/Splits/Cell17_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/20171116/Splits/Cell2_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/20171116/Splits/Cell6_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/TfR Puffs/20171116/Splits/Cell13_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/Splits/Cell1_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/Splits/Cell2_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/Splits/Cell5_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/Splits/Cell8_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/Splits/Cell18_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/2016 Washout Jenny/WTMOR/Splits/Cell22_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/WT/Crops/Splits/Cell11_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/WT/Crops/Splits/Cell12_0.1s/Ch1/*.tif*',
                    # '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/WT/Crops/Splits/Cell14_0.1s/Ch1/*.tif*',
                    '/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/WT/Crops/Splits/Cell3_0.1s/Ch1/*.tif*']
                    #'/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/WT/Crops/Splits/Cell4_0.1s/Ch1/*.tif*',
                    #'/Volumes/Coatamer/Users/weinberz/UM Drive/Data/Puffs Analysis/Rachel B2 Puffs/WT/Crops/Splits/Cell26_0.1s/Ch1/*.tif*']
    outfile_list = [#'TfRRachelCell8_events',
                    # 'TfRRachelCell10_events',
                    # 'TfRRachelCell14_events',
                    # 'TfRRachelCell17_events',
                    # 'TfRZYWCell2_events',
                    # 'TfRZYWCell6_events',
                    # 'TfRZYWCell13_events',
                    # 'MORCell1',
                    # 'MORCell2',
                    # 'MORCell5',
                    # 'MORCell8',
                    # 'MORCell18_events',
                    # 'MORCell22_events',
                    # 'B2Cell11_events',
                    # 'B2Cell12_events',
                    # 'B2Cell14_events',
                    'B2Cell3_events']
                    # 'B2Cell4_events',
                    # 'B2Cell26_events']


    for vid, name in enumerate(scorefile_list):
        print(name)
        for attempt in range(1):
            try:
                f = scorefile_list[vid]
                v = video_list[vid]
                o = outfile_list[vid]
                extract_data(f, v, o)
                break
            except Exception as e:
                print(e)
                print('failed')

if __name__ == '__main__': main()
