""" Process the foci tracked manually with MtrackJ. Look around the focus and center it on the local maximum.
    Then measure the ROI around it and plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from tifffile import imsave, imread
import pandas as pd
import os
import matplotlib.gridspec as gridspec
import re
import sys
import traceback
from scipy.signal import savgol_filter

from mtrackj_postprocessing import *


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def merge_kymos(location):
    for particle in range(1, 200):
        combined_per_channel = []
        for channel in [1, 2, 3]:
            combined = []
            for axis in ["x", "y", "z"]:
                filename = location + "/kymos_aligned/ch_" + str(channel) + "_track_1_" + str(particle) + "_" + axis + ".tif"
                
                if os.path.exists(filename):
                    img = imread(filename)
                    img = np.flip(img, axis = 0)
                    combined.append(img)

            

            if len(combined) > 0:
                combined_per_channel.append(np.vstack(combined))

        if len(combined_per_channel) > 0:
            Nc = len(combined_per_channel)
            Ny, Nx = np.shape(combined_per_channel[0])
            merged = np.zeros((Nc, Ny, Nx))

            for c in range(Nc):
                merged[c, :, :] = combined_per_channel[c]
            imsave(location + "/kymos_aligned/Composite_track_" + str(particle) + "_xyz.tif", merged)
                    
    


def get_tracks_from_file_3d(mtrackj_file):
    """ The input mtrackj_file is the result of MTrackJ.
        Output format is {t1: [x1, y1], t2: [x2, y2], ... tn: [xn, yn]} """
    
    tracks = []

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
                
                t = int(float(line[5]))
                x = float(line[2])
                y = float(line[3])
                z = float(line[4])
                
                tracks[-1][t] = [x, y, z, full_name]

    return tracks

#########################################################################################################################

def plot_kymo_rgb(filename, ax, width, valid_frames, num_channels, color_order):
    
    image = imread(filename)
    plt.cla()


    hist = plt.hist(image[0, :, :].flatten(), 100)
    counts, values,_ = hist

    counts = counts[10:]
    values = values[10:]
    max_0 = values[np.argmax(counts)]

    image[0, :, :] -= max_0
    image[0, :, :] = image[0, :, :] / np.max(image[0, :, :])

    if num_channels > 1:
        plt.cla()
        hist = plt.hist(image[1, :, :].flatten(), 100)
        counts, values,_ = hist

        counts = counts[10:]
        values = values[10:]
        max_1 = values[np.argmax(counts)]

        plt.cla()

        image[1, :, :] -= max_1
        image[1, :, :] = image[1, :, :] / np.max(image[1, :, :])

    if num_channels == 3:
        hist = plt.hist(image[2, :, :].flatten(), 100)
        counts, values,_ = hist
        max_2 = values[np.argmax(counts)]

        counts = counts[10:]
        values = values[10:]
    

    # Create the RGB image for plotting
    Nz, Ny, Nx = np.shape(image)
    rgb_image = np.zeros((3, Ny, Nx))
    
    rgb_image[color_order[0], :, :] = image[0, :, :]

    if num_channels > 1:
        rgb_image[color_order[1], :, :] = image[1, :, :]

    if num_channels == 3:
        
        image[2, :, :] -= max_2
        image[2, :, :] = image[2, :, :] / np.max(image[2, :, :])
        rgb_image[color_order[2], :, :] = image[2, :, :]

   
    rgb_image = np.swapaxes(rgb_image, 0, 2)
    rgb_image = np.swapaxes(rgb_image, 0, 1)


    plt.cla()

    ax.tick_params(axis='both', which='both', colors = "white") 
    ax.set_ylabel("merged", labelpad = 0)

    rgb_image = np.flipud(rgb_image)
    
    ax.imshow(rgb_image[:, width * np.min(valid_frames) : width * np.max(valid_frames), :], aspect=1,  origin='lower')

#########################################################################################################################

def plot_kymos(filenames, ax, width, valid_frames, ylabel, color):

    images = []
    boundary_pixels = []
    for filename in filenames:
        
        image = imread(filename)
        plt.cla()
        hist = plt.hist(image.flatten(), 100)
        counts, values,_ = hist

        counts = counts[10:]
        values = values[10:]
        max_0 = values[np.argmax(counts)]

        
        image -= max_0      
        image = image / np.max(image)

        images.append(image)

        # add a white row for easier visual separation
        (Ny, Nx) = np.shape(image)
        images.append(np.ones((1, Nx)))
        boundary_pixels.append(Ny)

    image_stack = np.vstack(images)
    
    
    plt.cla()
    ax.set_ylabel(ylabel, labelpad = 1)


    Ny, Nx = np.shape(image_stack)
    image_stack_rgb = np.zeros((Ny, Nx, 3))
    image_stack_rgb[:, :, color] = image_stack

    for c in [0, 1, 2]:
        image_stack_rgb[boundary_pixels[0], :, c] = 1
        image_stack_rgb[boundary_pixels[0] + boundary_pixels[1] + 1, :, c] = 1
        image_stack_rgb[boundary_pixels[0] + boundary_pixels[1] + boundary_pixels[2] + 2, :, c] = 0


    #image_stack_rgb = np.flipud(image_stack_rgb)

    ax.tick_params(axis='both', which='both', colors = "white")
    
    ax.imshow(image_stack_rgb[:, width * np.min(valid_frames) : width * np.max(valid_frames), :], aspect=1,  origin='lower')
    ax.tick_params(labelright=True)
    

#########################################################################################################################

def check_clusters(filename):
    all_files = [x for x in os.listdir(filename) if "Composite_track" in x]
    all_numbers = []
    all_strings = []
    
    for file in all_files:
        all_numbers.append ("".join(str([int(s) for s in file.split("_") if s.isdigit()])))
    
    all_numbers = list(np.unique(all_numbers))

    all_numbers = [x[1:-1] for x in all_numbers]
    


    
        
    all_numbers.sort(key=natural_keys)

    all_numbers = [x.replace(", ", "_") for x in all_numbers]


    return all_numbers

#########################################################################################################################

def get_pixel_size_and_frame_time(metadata_file, hour):
    metadata_csv = pd.read_csv(metadata_file, sep=",",
                           header=0, encoding='utf-8', engine='python', encoding_errors = "ignore")

    key = ""
    for col in metadata_csv.columns:
        if hour in col:
            key = col

    time_interval = int(float(metadata_csv[key][4][:-2]))
    objective = metadata_csv[key][15]

    pixel_size  = 0.121
    if '60x' in objective:
        pixel_size = 0.201
    else:
        pixel_size = 0.121
    return int(time_interval), pixel_size
    
#########################################################################################################################
    
def measure_speed(tracks_3d,
                  cell_folder,
                  stack_filename,
                  frame_time,
                  cluster_numbers,
                  kymos_folder,
                  color_order,
                  which_plots,
                  log,
                  num_channels_given = None,
                  plot_height = 15,
                  plot_angles = False,
                  smooth_window = 21,
                  xy_pixel_size = 1,
                  z_pixel_size = 1):

    """ which_plots is a list of strings, specifying which kymos to include, e.g. ["ch_1", "ch_2", "combined"] """
    from scipy.signal import savgol_filter
        
    output_folder = cell_folder  + "/displacement"

    stack = imread(stack_filename)

    # If the stack is missing a dimension (e.g. only 1 channel or 1 Z slice), add a singleton axis
    if len(np.shape(stack)) == 3:
        stack = stack[:, np.newaxis, :, :]
        
    numT, num_channels, height, width = np.shape(stack)
    print ("Channels: " + str(num_channels))
    

    if num_channels_given is not None:
        num_channels = num_channels_given
        
    #numT = int(numT / 5)

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    """ Add the speed of the particle to the tracks_3d object. Just 2D speeds for now """
    tracks_3d_copy = tracks_3d.copy()

    width = 5
    cluster = 1
    num_in_previous_cluster = 0

    if not os.path.exists(cell_folder + "/combined_plots/"):
        os.mkdir(cell_folder + "/combined_plots/")


    speeds = {}
    speeds_z = {}
    speeds_3d = {}
    speeds_real_units = {}
    speeds_real_units_filter = {}

    num_tracks = len(tracks_3d)
    
        
    for i, track in enumerate(tracks_3d):

        try:
            

            cluster_string = cluster_numbers[i]
            last_position = None
            last_position_z = None
            speeds[i] = []
            speeds_z[i] = []
            speeds_3d[i] = []
            speeds_real_units[i] = []

            log.write(cluster_string + "\n")

            
            valid_frames = []
            positions = []
            angles = []
            
            for t in range(numT):

                if t in track.keys():
                    valid_frames.append(t)
                    
                    x, y, z, name = track[t]

                    

                    if last_position is None:
                        last_position = np.array([y, x])
                        delta = 0

                    if last_position_z is None:
                        last_position_z = np.array([z])
                        delta_z = 0
                        
                    else:
                        delta = np.linalg.norm(np.array([y, x]) - last_position)
                        delta_z = np.linalg.norm(np.array([z]) - last_position_z)

                        last_position = np.array([y, x])
                        last_position_z = np.array([z])
                        speeds[i].append(1000 * convert_to_microns_per_second(delta, xy_pixel_size, frame_time)) # nm/s
                        speeds_real_units[i].append(convert_to_microns_per_second(delta, xy_pixel_size, frame_time))

                        speeds_z[i].append(1000 * convert_to_microns_per_second(delta_z, z_pixel_size, frame_time)) # nm/s

                        speeds_3d[i].append((speeds[i][-1] ** 2 + speeds_z[i][-1] ** 2) ** 0.5)
                        
                 
                    tracks_3d_copy[i][t] = x, y, z, name, delta, delta_z
                    positions.append(np.array([z, y, x]))

            # Measure the angles

            if plot_angles:

                x = int(x)
                y = int(y)
                z = 0 #4 * int(z)
                positions = np.array(positions)

                if smooth_window > 4:
                    smooth_length = min(smooth_window, len(positions) - 1)
                    if smooth_length % 2 == 0:
                        smooth_length -= 1
                    
                    positions_x = savgol_filter(positions[:, 2], window_length = smooth_length, polyorder = 4)
                    positions_y = savgol_filter(positions[:, 1], window_length = smooth_length, polyorder = 4)
                    positions = np.vstack([positions_y, positions_x])
                    
                    positions = np.transpose(positions)
                
                
                for p in range(np.shape(positions)[0]-2):
                    v1 = positions[p+1, :] - positions[p, :]
                    v2 = positions[p+2, :] - positions[p+1, :]

                    unit_vector_1 = v1 / (0.001 + np.linalg.norm(v1))
                    unit_vector_2 = v2 / (0.001 + np.linalg.norm(v2))
                    dot_product = np. dot(unit_vector_1, unit_vector_2)
                    angles.append(180.0 * np.arccos(dot_product) / 3.14)

            


            ####################################################################################################################################################################
                    
            fig, ax = plt.subplots(2 + num_channels, 1, sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [1] * (2 + num_channels)}, figsize = (0.5 + len(valid_frames)/6, plot_height))
            if num_channels == 2:
                fig, ax = plt.subplots(2 + num_channels, 1, sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [1] * (2 + num_channels)}, figsize = (0.5 + len(valid_frames)/8, plot_height + 1))
            ####################################################################################################################################################################
           
            plot_kymo_rgb(kymos_folder + "/" + which_plots[0] + "_" + cluster_string + "_xyz.tif", ax[1], width, valid_frames, num_channels, color_order)

            first_color, second_color, third_color = color_order

            plot_kymos([kymos_folder + "/ch_1_track_1" + "_" + cluster_string + "_z.tif",
                        kymos_folder + "/ch_1_track_1" + "_" + cluster_string + "_y.tif",
                        kymos_folder + "/ch_1_track_1" + "_" + cluster_string + "_x.tif"], 
                        ax[2], width, valid_frames, "ch 1", color = first_color)

            if num_channels > 1:

                plot_kymos([kymos_folder + "/ch_2_track_1" + "_" + cluster_string + "_z.tif",
                            kymos_folder + "/ch_2_track_1" + "_" + cluster_string + "_y.tif",
                            kymos_folder + "/ch_2_track_1" + "_" + cluster_string + "_x.tif"], 
                            ax[3], width, valid_frames, "ch 2", color = second_color)

            if num_channels == 3:
                plot_kymos([kymos_folder + "/ch_3_track_1" + "_" + cluster_string + "_z.tif",
                            kymos_folder + "/ch_3_track_1" + "_" + cluster_string + "_y.tif",
                            kymos_folder + "/ch_3_track_1" + "_" + cluster_string + "_x.tif"], 
                            ax[4], width, valid_frames, "ch 3", color = third_color)

            ax[0].set_ylabel("speed, nm/s")
            ax[0].scatter(range(int(width/2)+1, 5 * len(speeds[i]), width), speeds[i], c = "black", s = 5)
            ax[0].plot(range(int(width/2)+1, 5 * len(speeds[i]), width), speeds[i], c = "black", label = "XY")

            ax[0].scatter(range(int(width/2)+1, 5 * len(speeds[i]), width), speeds_3d[i], c = "red", s = 5)
            ax[0].plot(range(int(width/2)+1, 5 * len(speeds[i]), width), speeds_3d[i], c = "red", label = "XYZ")

            ax[0].scatter(range(int(width/2)+1, 5 * len(speeds[i]), width), speeds_z[i], c = "cyan", s = 5)
            ax[0].plot(range(int(width/2)+1, 5 * len(speeds[i]), width), speeds_z[i], c = "cyan", label = "Z")

            if plot_angles:
                ax[0].scatter(range(int(width/2)+2, 5 * len(angles), width), angles, c = "green", s = 5)
                ax[0].plot(range(int(width/2)+2, 5 * len(angles), width), angles, c = "green", label = "angle")

            leg = ax[0].legend()

            for color,text in zip(["black", "red", "cyan"], leg.get_texts()):
                text.set_color(color)

            ax[-1].set_xlabel("image number (x " + str(frame_time) + "s)")

            fig.tight_layout()
            #plt.show()
            plt.savefig(cell_folder + "/combined_plots/particle_" + cluster_string + ".svg")
            plt.cla()

        except Exception as e:
            log.write(str(e) + "\n")
            log.write(str(sys.exc_info()[2]))
            print(traceback.format_exc())

    # save the speeds
    df = pd.DataFrame(data = speeds.values(), index=cluster_numbers)
    df = df.T
    df.to_excel(cell_folder + "/speeds.xlsx")

    # save the speeds in Z
    df = pd.DataFrame(data = speeds_z.values(), index=cluster_numbers)
    df = df.T
    df.to_excel(cell_folder + "/speeds_z.xlsx")

    # save the speeds in 3d
    df = pd.DataFrame(data = speeds_3d.values(), index=cluster_numbers)
    df = df.T
    df.to_excel(cell_folder + "/speeds_3d.xlsx")

        #############



    return tracks_3d_copy
    
#########################################################################################################################

def convert_to_microns_per_second(delta, pixel_size, frame_time):
    # delta is given as pixels per frame; pixel_size is microns per pixel; frame_time is seconds per frame
    pixel_size = float(pixel_size)
    return delta * pixel_size / frame_time

##############################################################################################################################

def infer_color_order(metadata_file, order):
    metadata_csv = pd.read_csv(metadata_file, sep=",",
                           header=0, encoding='utf-8', engine='python', encoding_errors = "ignore")

    key = ""
    for col in metadata_csv.columns:
        if hour in col:
            key = col

    channel_1 = metadata_csv[key][6][:3]
    channel_2 = metadata_csv[key][7][:3]
    channel_3 = metadata_csv[key][8][:3]

    channel_1 = int(channel_1)
    channel_2 = int(channel_2)

    if channel_3 != " ":
        channel_3 = int(channel_3)
    else:
        channel_3 = 1000

    color_code = {}
    for i in range(1001):

        if i < 450:
            color_code[i] = 2 #"blue"
        elif i < 540:
            color_code[i] = 1 #"green"
        elif i < 630:
            color_code[i] = 0 #"red"
        else:
            color_code[i] = 2 #"magenta" but shown as blue
        

    
    return color_code[channel_1], color_code[channel_2], color_code[channel_3]

##############################################################################################################################


cell_folder = None
raw_data = None
max_projection_filename = None
mtrackj_file = None


# Defaults
num_channels_given = None
time_interval = 1
smooth_window = 1
xy_pixel_size = 1
z_pixel_size = 1


plot_angles = False
for i in range(len(sys.argv)):
    if sys.argv[i] == "-out":
        cell_folder = sys.argv[i+1]

    if sys.argv[i] == "-in":
        raw_data = sys.argv[i+1]

    if sys.argv[i] == "-max":
        max_projection_filename = sys.argv[i+1]

    if sys.argv[i] == "-tracks" or sys.argv[i] == "-track":
        mtrackj_file = sys.argv[i+1]
        
    if sys.argv[i] == "-num_channels":
        num_channels_given = int(sys.argv[i+1])
        
    if sys.argv[i] == "-plot_angles":
        plot_angles = True

    if sys.argv[i] == "-smooth_window":
        smooth_window = int(sys.argv[i+1])

    if sys.argv[i] == "-t":
        time_interval = int(sys.argv[i+1])

    if sys.argv[i] == "-xy":
        xy_pixel_size = sys.argv[i+1]

    if sys.argv[i] == "-z":
        z_pixel_size = sys.argv[i+1]


stack = load_stack(raw_data)

if not os.path.exists(cell_folder):
    os.mkdir(cell_folder)

# If the stack is missing a dimension (e.g. only 1 channel or 1 Z slice), add a singleton axis
if len(np.shape(stack)) == 4 and num_channels_given > 1:
    stack = stack[:, np.newaxis, :, :, :]

elif len(np.shape(stack)) == 4 and num_channels_given == 1:
    stack = stack[:, :, np.newaxis, :, :]

numT, numZ, num_channels, Ny, Nx = np.shape(stack)

plus_minus_frames = 5

# Create results folders

if not os.path.exists(cell_folder + "/kymos_aligned"):
    os.mkdir(cell_folder + "/kymos_aligned")

mtrackj_file_3d = mtrackj_file[:-4] + "_3d.mdf"

if True:
    
    (tracks_2d, tracking_channel) = get_tracks_from_file(mtrackj_file)
    print ("Tracking channel: " + str(tracking_channel))
    tracks_3d, signals = add_z_coordinate(tracks_2d, stack, tracking_channel-1, 3, numT)
    save_to_file(tracks_3d,
                 mtrackj_file,
                 mtrackj_file_3d)
    
    # save the signal profile
    df = pd.DataFrame(data = signals.values(), index=signals.keys())
    df = df.T
    df.to_excel(cell_folder + "/signals.xlsx")


if True:
    tracks_3d, _ = get_tracks_from_file(mtrackj_file_3d)

    for axis in ["y", "x"]:
        for channel in range(num_channels):
            
            print ("Processing " + axis + " projection, channel " + str(channel + 1))
            kymograph_aligned(stack[:, :, channel, :, :], tracks_3d, numT, cell_folder + "/kymos_aligned", channel + 1, main_axis = axis, plus_minus_frames = plus_minus_frames)
            kymograph_aligned_z(stack[:, :, channel, :, :], tracks_3d, numT, cell_folder + "/kymos_aligned", channel + 1,  plus_minus_frames = plus_minus_frames)

if True:
    tracks_3d = measure_speed_mtrackj(tracks_3d, numT, cell_folder + "/displacement", 0.12, 30)

print ("Generated individual kymograms in kymos_aligned.")

#####################################################################################################################################################

plt.close('all')

merge_kymos(cell_folder)


color_order = [0, 1, 2] #infer_color_order("metadata_summaries/metadata_summary_" + date + ".csv", hour)

log_file = open("log.txt", "a")
log_file.write("Started processing " + cell_folder + "\n")

cluster_numbers = check_clusters(cell_folder + "/kymos_aligned")

tracks_3d = get_tracks_from_file_3d(mtrackj_file_3d)
tracks_3d = measure_speed(tracks_3d,
                          cell_folder,
                          max_projection_filename,
                          time_interval,
                          cluster_numbers,
                          cell_folder + "/kymos_aligned/",
                          color_order,
                          which_plots = ["Composite_track"],
                          log = log_file,
                          num_channels_given = num_channels_given,
                          plot_angles = plot_angles,
                          smooth_window = smooth_window,
                          xy_pixel_size = xy_pixel_size,
                          z_pixel_size = z_pixel_size)


plt.close('all')
print ("Processing complete.")

