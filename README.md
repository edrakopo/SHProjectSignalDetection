# SHProjectSignalDetection - Analysing DT5742 Digitiser pulses


All functions and programs are within the Code folder.
The main functions that allow this code to run are within functions.py (and functions.html), and that is where you should look if you want to find useful methods and scripts for processing signal data.

All other files are just scripts developed that use the functions from functions.py to complete a task, as has been described below. Some of the descriptions may be slightly out of date.

If you have any questions, please feel free to email me.

CODE FILE DESCRIPTIONS
======================

eventviewer.py
--------------
Allows one to choose any event from a designated file and plot it based on its timegate.

functions.py
------------
File that holds many functions used throughout the code in an organised manner. Plan is to move everything to this file in the near future as currently the code in each file is messy.

noiseFT.py
----------
Opens up ROOT file, roughly removes signal events using sigma from the mean.
Finds best straight line best fit for events.
Then removes baseline and trendline from events.
Applies Butterworth transform to these events to remove high frequency noise.
Calculates FT for baseline-trendline removed events, and adds all FTs across all events to get a general picture of the frequencies dominant in the data.
Lots of code is commented out of this, as it is the oldest and messiest.


noiseFTLCMS.py
--------------
Functions much the same as noiseFT.py, but instead of using baselinesubtraction and LINFIT, uses LCMS, so determine if the trendline removal has any impact on the Fourier representation of our data.

noisetrimmer.py
---------------
Opens up ROOT file. Removes baseline and trendline from events. Applies butterworth filter to events.
Applies rolling mean to Events.
Collects properties of events to determine signals.
(Properties: risetime, FWHM, length, depth)
Warning! Length only works if you change onesigvals (line102) to have some sort of sigma-from-mean sorting. As it relies on the sigma difference to calculate signal length. (Will be changed).
Prints distribution of all Properties.
Using cutoff parameters, determines possible signal values within the data


noisetrimmeropt.py
------------------
For larger data files.
Removes baseline and trendline from events using LCMS algorithm.
Removes unneeded property calculations.
Properties Calculated: FWHM and Depth.
Using cutoff parameters, determines possible signal values within the data.
Calculates the number of signal events found within range determined by Y

roottest.py
-----------
OBSOLETE CODE FROM START OF PROJECT.
Only kept due to usefulness of some of the operations within it..

signalident.py
---------------
Works in the same manner as noisetrimmer.py, but only uses signal values (using 1sigma from mean method) to allow for better visualisation of signal property distributions.


LCMStest.py
-----------
Testing the different LCMS trendline removal functions against the linfit function to see what is faster.
