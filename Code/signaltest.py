import uproot
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc
import pandas as pd


from scipy import signal
from scipy.fft import rfft, rfftfreq
from scipy import stats

# PMTsignals/Run203-PMT107.root
# PMTsignals/Run203-PMT78.root
# PMTsignals/Run103-noise-PMT78.root
# PMTsignals/Run103-noise-PMT107.root



#file = uproot.open("PMTsignals/Run203-PMT107.root")
#print(file.keys())

#print(file["Tree"].values())
#print(file["Tree"].keys())


# Open the data, apply to variable
file = "PMTsignals/Run203-PMT78.root"

tree = uproot.open(file)["Tree"]
branches = tree.arrays()
print(branches['ADC'])


time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(150):
    time.append(i*2)



# Plot the first 50 events
#for i in range(50):
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

# baseline subtraction, do in the past step, go through each file and find the 'linear trend'
# remove baseline (after making sure it isnt a signal) and trend data
# if it is a signal, use baseline from previous attempt to reduce it

# then use pyFFTW to do fourier transform of our data, possibly FT for multiple events all together (superimposed)


# Mean,std, min, max, data collection

meanvals = []
stdvals = []
minvals = []
maxvals = []
sigvals = []
# Testing range of 50 values
for i in range(50):
    datarepeated = branches['ADC'][i]
    # Append data values to list
    meanvals.append(np.mean(datarepeated))
    stdvals.append(np.std(datarepeated))
    minvals.append(np.amin(datarepeated))
    maxvals.append(np.amax(datarepeated))

    # If mean - min val is larger than 5 sigma, mark as signal
    # If max val - mean is larger than 5 sigma, mark as signal
    if (meanvals[i] - minvals[i]) > (5*stdvals[i]) or (maxvals[i] - meanvals[i]) > (5*stdvals[i]):
        # Signal
        sigvals.append(1)

    else:
        # No signal
        sigvals.append(0)

    # Uncomment to show all the signal values
    #if sigvals[i] == 1:
    #    plt.plot(time,datarepeated)
    #    plt.xlabel("Sample Time (ns)")
    #    plt.ylabel("ADC Value")
    #    plt.title(str(file) + " event " + str(i))
    #    plt.show()






# Now collect all 'possible signals', will be useless when dealing with noise but useful right now
# Plan, keep array of same length at normal array with 0 or 1, determines if the array is

print("Mean values: \n" + str(meanvals))
print("Standard Deviations: \n" + str(stdvals))
print("Minimum Values: \n" + str(minvals))
print("Maximum Values: \n" + str(maxvals))
print("Signal Values: \n" + str(sigvals))
