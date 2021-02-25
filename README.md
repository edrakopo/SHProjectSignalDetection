# SHProjectSignalDetection


All functions and programs are within the Code folder.

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
Removes baseline from events.
Removes unneeded property calculations.
Properties Calculated: FWHM and Depth.
Using cutoff parameters, determines possible signal values within the data

roottest.py
-----------
OBSOLETE CODE FROM START OF PROJECT.
Only kept due to usefulness of some of the operations within it..

signalident.py
---------------
Works in the same manner as noisetrimmer.py, but only uses signal values (using 1sigma from mean method) to allow for better visualisation of signal property distributions.
