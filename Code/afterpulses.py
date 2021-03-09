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
import xlsxwriter

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

file = "E:\PMTsignals\Boulby_78_After.root"
# Open the data, apply to variable
datafull = fnc.rootopen(file)


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
y = 5000 #segments
ys = 50000 # full list size

afterpulsetimes = []
# the code doesnt like being iterated over large lists, so take smaller lists and compile the results
for z in range(int(ys/y)):
    zzz = ys//y
    print("Loop " + str(z) + " of " + str(zzz))
    # take the y values from data, to stop array index mismatching
    data = datafull[y*z:y*(z+1)]
    print("Covering data from " + str(y*z) + " to " + str(y*(z+1)))

    # Apply full data manipulation, all modifications applied
    print("Sampling Data")
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



    print("Finding Afterpulses...")


    # modify signal events by shifting the each event to the left until initial signal is gone


    newerarray = []

    siglen = len(signalevents)
    # loop across all events
    sigarray = [0] * siglen
    for i in range(siglen):

        # peak finding reset variable
        peakreached = 0
        # pull out current event of interest
        event = signalevents[i]
        # set lowest value in event
        minval = np.amin(event)
        # loop across all samples
        for j in range(samples):
            # If you have reached the signal peak, look for signal cutoff
            if (event[j] == minval):
                peakreached = 1

            # if you've passed the peak, and have most likely passed , collect the new data
            if (peakreached == 1) and (event[j] > 0): # only enable when wanting more accuracy
                # to deal with array shape mismatches, make new array that is samples size but 0s everywhere except beyond j
                sigarray[i] = event[j:]
                sigarray[i].extend([0]*j)
                peakreached = 2

        # Progress
        if (i%10) == 0:
            print("Afterpulse Process: " + str(i) + "/" + str(siglen))


    #plt.plot(time,sigarray[0])
    #plt.xlabel("Sample Time (ns)")
    #plt.ylabel("ADC Value")
    #plt.title("Afterpulse event")
    #plt.show()

    # Collect afterpulses
    cutdata2, afterpulseevents = fnc.signalspotter(sigarray,timegate,-30,(7.5,20.5))

    #plt.plot(time,afterpulseevents[0])
    #plt.xlabel("Sample Time (ns)")
    #plt.ylabel("ADC Value")
    #plt.title("Afterpulse event")
    #plt.show()


    # Take the afterpulse data, and find when the peak is on the
    afterpulselen = len(afterpulseevents)
    # list of afterpulse times


    for k in range(afterpulselen):
            event = afterpulseevents[k]
            # find peak
            minval = np.amin(event)

            for l in range(samples):
                # if at the peak, take the value of time
                if (event[l] == minval):
                    afterpulsetimes.append(l*timegate)


plt.hist(afterpulsetimes, bins = 50)
plt.title("Afterpulse histogram across time for " + str(file))
plt.xlabel("Time (ns)")
plt.ylabel("Distribution")
plt.show()

np.savetxt('Afterpulse_values', afterpulsetimes)
