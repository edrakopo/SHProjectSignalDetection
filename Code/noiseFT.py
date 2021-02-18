import uproot
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc
import pandas as pd
import random
import pyfftw

from scipy import signal
from scipy.fft import rfft, rfftfreq
from scipy import stats

################ DESCRIPTION
#
# Opens up ROOT file, roughly removes signal events
# Finds best straight line best fit for events
# Then removes baseline and trendline from events
# Calculates FT for baseline-trendline removed events
################


# PMTsignals/Run203-PMT107.root
# PMTsignals/Run203-PMT78.root
# PMTsignals/Run103-noise-PMT78.root
# PMTsignals/Run103-noise-PMT107.root

#file = uproot.open("PMTsignals/Run203-PMT107.root")
#print(file.keys())

#print(file["Tree"].values())
#print(file["Tree"].keys())


# Open the data, apply to variable
file = "E:\PMTsignals\Run203-PMT78.root"

tree = uproot.open(file)["Tree"]
branches = tree.arrays()
print(branches['ADC'])

# how long between data taken
timegate = 2
# length of event
eventno = len(branches['ADC'][0])

time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(eventno):
    time.append(i*timegate)


print("Event number: ")
print(str(len(branches['ADC'])))
# Plot the first 50 events
for i in range(10000):

    datarepeated = branches['ADC'][i]
    plt.plot(time,datarepeated)
    plt.xlabel("Sample Time (ns)")
    plt.ylabel("ADC Value")
    plt.title(str(file) + " event " + str(i))
    plt.show()


# So we have in the file:

# -["Tree"]
#   -['ADC']
#     -[1]
#     -[2]
#     -[...]
#     -[10000]


# PLAN
# Take mean, sigma/stdev, min, max values of each event
# if min or max is 5 larger than stdev from mean, put in a new list as a "possible signal"

# baseline subtraction, done in the previous file (roottest.py), go through each file and find the 'linear trend'
# remove baseline (after making sure it isnt a signal) and trend data
# if it is a signal, use baseline from previous attempt to reduce it

# then use pyFFTW to do fourier transform of our data, possibly FT for multiple events all together (superimposed)


y = 100000
print("Calculating " + str(y) + " events")
# Collect important values about a certain number of events (y)
meanvals, stdvals, minvals, maxvals, sigvals, medvals = fnc.datacollate(branches['ADC'], y)

print("Data collation complete")
######### Show signal values and lists of all important values, will usually be commented out for faster runtimes
#for i in range(y):
#    if sigvals[i] == 1:
#        plt.plot(time,branches['ADC'][i])
#        plt.xlabel("Sample Time (ns)")
#        plt.ylabel("ADC Value")
#        plt.title(str(file) + " event " + str(i))
#        plt.show()

#print("Mean values: \n" + str(meanvals))
#print("Standard Deviations: \n" + str(stdvals))
#print("Minimum Values: \n" + str(minvals))
#print("Maximum Values: \n" + str(maxvals))
#print("Signal Values: \n" + str(sigvals))
#########

########### LINEAR FITTING ###########


# Line of best fit, to see modulation in baseline
pfit, stats, rms, c, m = fnc.linfit(time,branches,y)
# pfit is quite complicated, if given more time at the end REVISE THIS BIT OF CODE to not need c or m, but just pfit
print("Linear Fitting complete")
# Create distribution of 10th event - RAW DATA, signal removed
##datanosig = []
##for i in range(y):
    # if not a signal
##    if sigvals[i] == 0:
##        datanosig.append(branches['ADC'][i])
##eventdistr = fnc.adcdist(datanosig,y)

######## Create line to plot event RANDOM, commented out for run time
yline = []

# ensure random plot isn't a signal
while True:
    q = random.randint(0,y-1)
    if sigvals[q] == 0:
        break



for j in range(len(time)):
    yline.append(m[q]*time[j]+c[q])
print("Linear Fit Applied")
# Plot line of best fit over data
plt.plot(time,branches['ADC'][q],color='black',markersize=2)
plt.plot(time,yline,color='red',linewidth=3)
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.title("Trendline of " + str(file) + " event " + str(q) )
plt.show()


############# BASELINE/LINEAR TREND REMOVAL ###############
# PACK ALL OF THIS INTO FUNCTIONS.PY SOON!

# moving data to be adjusted, data[i] is the ith event

#scaleddata = [None] * y
lindata = [None] * y
newdata = []

# take y values from the array
sliceddata = branches['ADC'][:y]
# baseline subtraction
scaleddata = fnc.baselinesubtraction(sliceddata, c)

# Removal of baseline and lintrend
for i in range(y):
    # blocking removal of signal baseline, will be rewritten in the future

    # no signal, old baseline removal code
    #if sigvals[i] == 0:
        # remove baseline
        #print("event " + str(i) + ": \n" + str(branches['ADC'][i]))
        #print("subtract " + str(c[i]))
    #    scaleddata[i] = branches['ADC'][i] - c[i]
        #print("baseline removed event" + str(i) + ": \n" + str(scaleddata[i]))

        # Remove linear trend
        for j in range(len(time)):
            # Write new data to dummy list
            newdata.append(scaleddata[i][j] - m[i]*time[j])
        # Apply to real list
        lindata[i] = newdata
        # Refresh dummy list
        newdata = []
        #print("event " + str(i) + " complete. Continuing...")
    # signal. work on this later!
    #elif sigvals[i] ==1:
    #    print("event " + str(i) + " Signal! Do not process yet!")
print("Completed Trendline Removal")


# To ensure its not a signal variable
if sigvals[q] == 0:

    # Plot trendline/baseline removed data
    plt.plot(time,lindata[q],color='black',markersize=2)
    plt.xlabel("Sample Time (ns)")
    plt.ylabel("ADC Value")
    plt.title("Trendline and Baseline removed " + str(file) + " event " + str(q) )
    plt.show()

# Create distribution of 10th event - Trendline/baseline removed
# filter out signal values (Nonetype right now) from eventdistr
##lindatanosig = []
##for i in range(0,len(lindata)):
##    if lindata[i] != None:
##        lindatanosig.append(lindata[i])
##eventdistr = fnc.adcdist(lindatanosig,len(lindatanosig))




# check if linear trend for lindata[q] is now 0. test variables, no longer needed
#Polynomial = np.polynomial.Polynomial
#testpfit, teststats = (Polynomial.fit(time, lindata[q], 1, full=True, window=(0, 150), domain=(0, 150)))
#print("c: " + str(c[q]) + "\nm: " + str(m[q]))
#print(testpfit)

# finding DFT for noise events using pyFFTW
# pyFFTW creation
dftdata= [None] * y
freqdata = [None] * y

# Create FT data for random variable q, for testing.
dftdata[q] = pyfftw.interfaces.numpy_fft.rfft(lindata[q])
freqdata[q] = pyfftw.interfaces.numpy_fft.rfftfreq(len(time),1/500)

plt.plot(freqdata[q],np.abs(dftdata[q]))
plt.xlabel("Sample Frequency (MHz)")
plt.ylabel("Amplitude")
plt.title("Fourier transform of event " + str(q) + " - File " + str(file))
plt.show()

# Create FT of all data
dfthist = [0] * len(dftdata[q])

for i in range(y):
    # ensure no signal collected for now
    if sigvals[i] == 0:
        # create fourier transform here
        dftdata[i] = pyfftw.interfaces.numpy_fft.rfft(lindata[i])
        # add to summation array
        dfthist = np.add(dfthist,dftdata[i])

plt.plot(freqdata[q],np.abs((dfthist).imag))
plt.xlabel("Sample Frequency (MHz)")
plt.ylabel("Amplitude")
plt.title("Fourier Transform of all events (additive) of file " + str(file))
plt.show()
