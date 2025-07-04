# Requirements

1.	Generate a maximum intensity projection of the tif image stack and save it. Then track the particles using MTrackJ and save the tracks in a .mdf file.

2.	SPARTACUSS is python application. Please make sure you have python running on your system. This can be checked in a command terminal using the command

SPARTACUSS was tested on Windows 7, Windows 10, Windows 11, and Mac OS Sonoma.

python --version

Before using SPARTACUSS, please install the required python libraries. For example, using pip, type the following command in a command window:
pip install scipy matplotlib tifffile pandas openpyxl scikit-image
This installation step takes one minute on an average computer.

In case pip is not recognized, you need to install it as well, using these commands:

curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py

python get-pip.py

# Usage 

This section explains the usage in general. Refer to the section “How to reproduce” below for a specific example.
Please open a command window and start the script with the following arguments:
python spartacuss.py -in "path/to/input/stack.tif" -out "results_folder" -max "path/to/MAX intensity projection.tif" -tracks "path/to/mtrackj results.mdf"
Important: enclose filenames in quotation marks, so they are interpreted as a single string in case they contain spaces.
The above arguments are required. In addition, you can provide the following optional arguments:

-xy <pixel size in XY in microns per pixel>

-z <pixel size in Z in microns per pixel>

-t <time interval between frames, in seconds>

These are important for the correct scale of speed measurements. If not provided, they default to 1.

# How to reproduce

In the provided example, there are 5 files:
VLP example raw.tif – the raw image stack
VLP example MAX.tif – the maximum intensity projection along Z
spartacuss.py, mtrackj_postprocessing.py – the scripts
example_track.mdf

Please place these files in a folder and open a command window from inside this folder. Then use the following command:
python spartacuss.py -in "VLP example raw.tif" -max "VLP example MAX.tif" -out "results" -tracks "example_track.mdf" -xy 0.1 -z 0.2
A new folder named “results” will be generated and it will contain the output.
In case of only one channel, please provide the additional argument -num_channels 1
This is important because the program assumes two channels by default.

# Expected results

The output will be available in the folder designated by the -out argument. Inside, you will find the following results:

1.	combined_plots – this folder contains one .svg file per particle. On the top, the speed in 2D and 3D is displayed over time. On the bottom, the kymographs per channel, as well as the merged channels, are displayed.
2.	displacement – this folder contains one .png file per particle with a plot of the 3D speed over time, similar to the plot in combined_plots/ above, but also containing a smooth curve.
3.	Kymos_aligned – this folder contains intermediate files - the individual kymographs per channel and per X, Y, Z axis, which are used in combined_plots/ above. These are tif files and so they can be opened in another software, e.g. Fiji, to do further processing (for example, adjusting brightness and contrast, cropping).
4.	signals.xlsx – contains the intensity values of the particles over time, where each particle is a separate column.
5.	Speeds.xlsx, speeds_z.xlsx, speeds_3d.xlsx – the actual values of the particle speeds used in the plots, where each particle is a separate column.
   
The run time of the pipeline depends on the number of tracks and the size of the raw data. For the example provided, the processing takes around one minute on an average computer.
