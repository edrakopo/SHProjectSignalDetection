import uproot
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc
import pandas as pd
import random

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

# baseline subtraction, done in the previous file (roottest.py), go through each file and find the 'linear trend'
# remove baseline (after making sure it isnt a signal) and trend data
# if it is a signal, use baseline from previous attempt to reduce it

# then use pyFFTW to do fourier transform of our data, possibly FT for multiple events all together (superimposed)


y = 100
# Collect important values about a certain number of events (y)
meanvals, stdvals, minvals, maxvals, sigvals = fnc.datacollate(branches, y)


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





######## Create line to plot event RANDOM, commented out for run time
yline = []
q = random.randint(0,y)

for j in range(len(time)):
    yline.append(m[q]*time[j]+c[q])

# Plot line of best fit over data
plt.plot(time,branches['ADC'][q],color='black',markersize=2)
plt.plot(time,yline,color='red',linewidth=3)
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.title("Trendline of " + str(file) + " event " + str(q) )
plt.show()


############# BASELINE/LINEAR TREND REMOVAL ###############

# moving data to be adjusted, data[i] is the ith event
scaleddata = [None] * y
lindata = [None] * y
newdata = []

# Removal of baseline and lintrend
for i in range(y):
    # blocking removal of signal baseline, will be rewritten in the future

    # no signal
    if sigvals[i] == 0:
        # remove baseline
        print("event " + str(i) + ": \n" + str(branches['ADC'][i]))
        print("subtract " + str(c[i]))
        scaleddata[i] = branches['ADC'][i] - c[i]
        print("baseline removed event" + str(i) + ": \n" + str(scaleddata[i]))

        # Remove linear trend
        for j in range(len(time)):
            # Write new data to dummy list
            newdata.append(scaleddata[i][j] - m[i]*time[j])
        # Apply to real list
        lindata[i] = newdata
        # Refresh dummy list
        newdata = []
    # signal. work on this later!
    elif sigvals[i] ==1:
        print("event " + str(i) + " Signal! Do not process yet!")
print("Completed")


# To ensure its not a signal variable
if sigvals[q] ==0:

    # Plot trendline/baseline removed data
    plt.plot(time,lindata[q],color='black',markersize=2)
    plt.xlabel("Sample Time (ns)")
    plt.ylabel("ADC Value")
    plt.title("Trendline and Baseline removed " + str(file) + " event " + str(q) )
    plt.show()

# check if linear trend for lindata[q] is now 0. test variables, no longer needed
#Polynomial = np.polynomial.Polynomial
#testpfit, teststats = (Polynomial.fit(time, lindata[q], 1, full=True, window=(0, 150), domain=(0, 150)))
#print("c: " + str(c[q]) + "\nm: " + str(m[q]))
#print(testpfit)
