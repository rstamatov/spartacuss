""" Process the foci tracked manually with MtrackJ. Look around the focus and center it on the local maximum.
    Then measure the ROI around it and plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tifffile import imsave, imread
import pandas as pd
import os

def get_tracks_from_file(mtrackj_file):
    """ The input mtrackj_file is the result of MTrackJ.
        Output format is {t1: [x1, y1], t2: [x2, y2], ... tn: [xn, yn]} """
    
    tracks = []
    tracking_channel = None

    with open(mtrackj_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()

            if line[0] == 'Track':
                new_track = {}
                tracks.append(new_track)

            if line[0] == "Cluster":
                cluster = line[1]
            if line[0] == "Track":
                name = line[1]

            if line[0] == 'Point':

                full_name = cluster + "_" + name

                tracking_channel = int(float(line[6]))
                t = int(float(line[5]))
                x = float(line[2])
                y = float(line[3])
                z = float(line[4])
                
                tracks[-1][t] = [x, y, z, full_name]

    return tracks, tracking_channel

#########################################################################################################################

def mtrackj_to_trackmate(tracks_3d, output_filename):
    """ Generate a TrackMate XML file, so the tracks can be loaded in TrackMate.
        The idea is to use the 3d viewer and see the raw image with the tracks. """

    nTracks = len(tracks_3d)
    pixel_size_xy = 0.1205860
    pixel_size_z = 1 #0.0958333
    
    with open(output_filename, "w") as f:

        f.write("TRACK_ID,POSITION_X,POSITION_Y,POSITION_Z,FRAME\n")

        for i, track in enumerate(tracks_3d):

            for t in track.keys():
                point = track[t]
                x, y, z = point
                x_adj = x * pixel_size_xy
                y_adj = y * pixel_size_xy
                z_adj = z * pixel_size_z

                f.write(str(i) + "," + str(x_adj) + "," + str(y_adj) + "," + str(z_adj) + "," + str(t-1) + "\n" )


############################################################################################################################
            
def kymograph_aligned(stack, tracks_3d, numT, output_folder, channel, main_axis = "y", plus_minus_frames = None):

    """ Same as kymograph() but all starting from 0 """

    kymos = {}
    names = {}

    if main_axis == "y":
        dz, dy, dx = 0, 10, 2
    else:
        dz, dy, dx = 0, 2, 10
    
    # For each time point, load the Z stack
    for t in range(1, numT):

        raw_img = stack[t-1, :, :, :]
        raw_img = np.pad(raw_img, ((dz, dz), (dy, dy), (dx, dx)), 'empty')

        # some dimensions
        Nz, Ny, Nx = np.shape(raw_img)

        for i, track in enumerate(tracks_3d):

            if t == 1:
                kymos[i] = np.zeros((Nz, 0))

            kymo_slice = np.zeros((Nz, 5))

            if t in track.keys():       # get the current position 
                x, y, z, name = track[t]
                
            elif t < list(track.keys())[0]:   # if before the start, get the first position
                x, y, z, name = list(track.values())[0]
                
            else:                       # if after the end, get the last position
                x, y, z, name = list(track.values())[-1]

            x = dx + int(x)
            y = dy + int(y)
            z = dz + int(z)

            chunk = raw_img[:, y - dy : y + dy + 1, x - dx : x + dx + 1]

            # New functionality - untested
            if plus_minus_frames != None:                
                if t < list(track.keys())[0] - plus_minus_frames or t > list(track.keys())[-1] + plus_minus_frames:
                    chunk = np.zeros_like(chunk)            
            # end untested

            if main_axis == "y":
                try:
                    kymo_slice = np.max(chunk, axis = 1) # Max projection
                except Exception as e:
                    print(e)

            elif main_axis == "x":

                try:
                    kymo_slice = np.max(chunk, axis = 2) # Max projection
                except Exception as e:
                    print (e)
            		
            kymos[i] = np.hstack((kymos[i], kymo_slice))
            names[i] = name

    for i in range(len(tracks_3d)):
        imsave(output_folder + "/ch_" + str(channel) + "_track_" + names[i] + "_" + main_axis + ".tif", np.array(kymos[i], dtype = np.float64), bigtiff = True)


    
############################################################################################################################

def kymograph_aligned_z(stack, tracks_3d, numT, output_folder, channel, plus_minus_frames = None):

    """ Same as kymograph() but all starting from 0 """

    kymos = {}
    names = {}

    dy, dx, dz = 5, 2, 5
    
    # For each time point, load the Z stack
    for t in range(1, numT):
        
        raw_img = stack[t-1, :, :, :]
        raw_img = np.pad(raw_img, ((dz, dz), (dy, dy), (dx, dx)), 'empty')

        # some dimensions
        Nz, Ny, Nx = np.shape(raw_img)

        for i, track in enumerate(tracks_3d):

            if t == 1:
                kymos[i] = np.zeros((2 * dy + 1, 0))

            kymo_slice = np.zeros((2 * dy + 1, 2 * dx + 1))

            if t in track.keys():       # get the current position 
                x, y, z, name = track[t]
                
            elif t < list(track.keys())[0]:   # if before the start, get the first position
                x, y, z, name = list(track.values())[0]
                
            else:                       # if after the end, get the last position
                x, y, z, name = list(track.values())[-1]

            x = dx + int(x)
            y = dy + int(y)
            z = dz + int(z)


            chunk = raw_img[:, y - dy : y + dy + 1, x - dx : x + dx + 1]

            # New functionality - untested
            if plus_minus_frames != None:                
                if t < list(track.keys())[0] - plus_minus_frames or t > list(track.keys())[-1] + plus_minus_frames:
                    chunk = np.zeros_like(chunk)            
            # end untested


            try:
                kymo_slice = np.max(chunk, axis = 0) # Max projection
            except Exception as e:
                print(e)

            kymos[i] = np.hstack((kymos[i], kymo_slice))
            names[i] = name

    for i in range(len(tracks_3d)):
        imsave(output_folder + "/ch_" + str(channel) + "_track_" + names[i] + "_z" + ".tif", np.array(kymos[i], dtype = np.float64), bigtiff = True)

#########################################################################################################################

def add_z_coordinate(tracks_2d, stack, tracking_channel, tolerance, numT):
    """ Find the max pixel along Z for each x, y. Allow a tolerance around X, Y
        "tracks" is a list of [[x1, y1], [x2, y2], ...] for each track.
        Output: a list of [[x1, y1, z1], [x2, y2, z2], ...] for each track. """

    #print ("Adding the Z-coordinate of each track. This will take a minute. We've been working on VLPs for 2 years, so a few minutes won't make a difference. ")
    
    tracks_3d = tracks_2d.copy()
    signals = {}
    # For each time point, load the Z stack
    for t in range(numT-1):

        t_plus_one = t + 1 # MTrackJ starts counting from 1, not 0

        raw_img = stack[t, :, tracking_channel, :, :]

        # If the stack is 2D, add a singleton Z
        if len(np.shape(raw_img)) == 2:
            raw_img = raw_img[np.newaxis, :, :]

        # some dimensions
        Nz, Ny, Nx = np.shape(raw_img)
        Ny = np.arange(0, Ny)
        Nx = np.arange(0, Nx)

        for i, track in enumerate(tracks_2d):

            if i not in signals.keys():
                signals[i] = []
            if t_plus_one in track.keys():
                
                x, y, z, name = track[t_plus_one]

                x = int(x)
                y = int(y)

                # Find the Z of the maximum value in a cylinder around (X, y)      
                roi_pixels = (Ny[:,np.newaxis] - y)**2 + (Nx[np.newaxis,:] - x)**2 <= tolerance**2

                z_values = np.max(raw_img[:, roi_pixels], axis = 1)
                signals[i].append(np.sum(z_values))
                #z = np.argmax(raw_img[:, y, x])

                
                z = np.argmax(z_values)
 
                tracks_3d[i][t_plus_one][2] = z
                
    return tracks_3d, signals

    
##############################################################################################################################

def save_to_file(tracks_3d, original_file, new_filename):
    """ tracks_3d will contain the added Z coordinate. Rewrite the original_file with that Z,
        and save as new_filename. """

    new_lines = []
    with open(original_file, "r") as f_read:
        lines = f_read.readlines()
        track_counter = -1
        
        for line in lines:

            line_list = line.split()

            if line_list[0] == 'Track':
                track_counter += 1
                new_line = line

            elif line_list[0] == 'Point':                
                [name, number, x, y, z, t, c] = line_list
                coords = tracks_3d[track_counter][int(float(t))]

                if len(coords) == 4:
                    x, y, z, name = coords
                else:
                    x, y, name = coords
                    z = -1
                    
                new_line = "Point" + " " + number + " " + str(x) + " " + str(y) + " " + str(z) + " " + str(t) + " " + str(c) + "\n"
            else:
                new_line = line
                
            new_lines.append(new_line)

    with open(new_filename, "w") as f_write:
        for line in new_lines:
            f_write.write(line)

#########################################################################################################################

def measure_speed_mtrackj(tracks_3d, numT, output_folder, pixel_size, frame_time):

    from scipy.signal import savgol_filter

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    """ Add the speed of the particle to the tracks_3d object. Just 2D speeds for now """
    tracks_3d_copy = tracks_3d.copy()

    speeds = {}
    speeds_real_units = {}
    
    
    for i, track in enumerate(tracks_3d):
        last_position = None
        speeds[i] = []
        speeds_real_units[i] = []
        
        for t in range(numT):

            if t in track.keys():
                
                x, y, z, name = track[t]

                x = int(x)
                y = int(y)
                z = int(z)

                if last_position is None:
                    last_position = np.array([y, x])
                    delta = 0
                else:
                    delta = np.linalg.norm(np.array([y, x]) - last_position)
                    last_position = np.array([y, x])
                    speeds[i].append(delta)
                    speeds_real_units[i].append(convert_to_microns_per_second(delta, pixel_size, frame_time))
                    
                tracks_3d_copy[i][t] = x, y, z, name, delta

        plt.clf()
        x_real = [j * frame_time for j in range(len(speeds[i]))]
        plt.plot(x_real, speeds_real_units[i], 'p')

        if len(speeds_real_units[i]) > 30:
            plt.plot(x_real, savgol_filter(speeds_real_units[i], 15, 3))
        plt.ylim(0, 0.75)
        plt.xlabel("Time since start, s")
        plt.ylabel("speed, microns/s")
        plt.title(name)
    
        plt.savefig(output_folder + "/particle_" + str(i) + ".png")

    return tracks_3d_copy

#########################################################################################################################

def convert_to_microns_per_second(delta, pixel_size, frame_time):
    # delta is given as pixels per frame; pixel_size is microns per pixel; frame_time is seconds per frame
    return delta * pixel_size / frame_time
        
##############################################################################################################################

def load_stack(filename):
    stack = imread(filename)
    return stack


