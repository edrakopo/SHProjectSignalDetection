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
file = "E:\PMTsignals\Run103-noise-PMT78.root"

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
#for i in range(10000):

#    datarepeated = branches['ADC'][i]
#    plt.plot(time,datarepeated)
#    plt.xlabel("Sample Time (ns)")
#    plt.ylabel("ADC Value")
#    plt.title(str(file) + " event " + str(i))
#    plt.show()


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

# defines how many events you wish to sample across
y = 100000
print("Calculating " + str(y) + " events")


# Apply LCMS to data
lindata = fnc.LCMSlist(branches['ADC'][:y])

# finding DFT for noise events using pyFFTW
# pyFFTW creation
dftdata= [None] * y
freqdata = [None] * y

# Create FT data for random variable q, for testing.
dftdatapoint = pyfftw.interfaces.numpy_fft.rfft(lindata[0])
freqdata = pyfftw.interfaces.numpy_fft.rfftfreq(len(time),1/500)

# Create FT of all data
dfthist = [0] * len(dftdatapoint)

for i in range(y):
    # create fourier transform here
    dftdata[i] = pyfftw.interfaces.numpy_fft.rfft(lindata[i])
    # add to summation array
    dfthist = np.add(dfthist,dftdata[i])

plt.plot(freqdata,np.abs((dfthist).imag))
plt.xlabel("Sample Frequency (MHz)")
plt.ylabel("Amplitude")
plt.title("Fourier Transform of all events (additive) of file " + str(file))
plt.show()
