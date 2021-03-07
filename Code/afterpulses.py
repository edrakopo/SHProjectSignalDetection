import uproot
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc
import pandas as pd
import random
import pyfftw
import probfit
import statistics
from time import process_time

from scipy import signal
from scipy.fft import rfft, rfftfreq
from scipy import stats


###############################
#
# DESCRIPTION
#
# Signal spotting code, but every time it spots a signal, add all the ADC values after the signal to another 'event' that
# collects all post-signal event data. Then graph this along time and see if you can see any trends.
#
###############################


# Open the data, apply to variable
datafull = fnc.rootopen("E:\PMTsignals\Boulby_78_After.root")


# how many values in each event
samples = len(datafull[0])
# how long between samples within events
timegate = 2

# Time axis for
time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(samples):
    time.append(i*2)

# Event control, how many events do you want to process?
y = 50000 #len(datafull)

# take the y values from data, to stop array index mismatching
data = datafull[:y]


# Apply full data manipulation, all modifications applied
print("Sampling Data ")
modifieddata = fnc.scaledata(data, 200, 5, butter = True, rolling =  True, lcms =  True)

# Spot signal and collate data
cutdata, signalevents = fnc.signalspotter(modifieddata,timegate,-30,(7.5,20.5))

# Look at all signal events
#for i in range(len(signalevents)):
#    plt.plot(time,signalevents[i])
#    plt.xlabel("Sample Time (ns)")
#    plt.ylabel("ADC Value")
#    plt.title("Afterpulse signal event " + str(cutdata[i]))
#    plt.show()

# Afterpulse addition function

# Create the afterpulse array
# working under the assumption that samples after signal < total samples
afterpulsedata = [0] * samples

# Calculate the lowest ADC(y) values in each signal event


# loop across all events
for i in range(len(signalevents)):

    # variable that defines whether a sample should be added to the list
    passpeak = 0

    # pull out current event of interest
    event = signalevents[i]
    # set lowest value in event
    minval = np.amin(event)

    # loop across all samples
    for j in range(samples):
        # If you have reached the signal peak, copy data to afterpulsedata
        if (event[j] == minval):
            # to deal with array shape mismatches, make new array that is samples size but 0s everywhere except beyond j
            newarray = [0] * j
            newerarray = np.hstack((event[j:],newarray))
            # add element-wise the events onto afterpulsedata list
            afterpulsedata = np.add(afterpulsedata,newerarray)

plt.plot(time,np.abs(afterpulsedata))
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.title("Afterpulse Additive")
plt.show()
