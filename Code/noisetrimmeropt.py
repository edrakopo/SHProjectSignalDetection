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

from iminuit import Minuit
from scipy import signal
from scipy.fft import rfft, rfftfreq
from scipy import stats

#collect data

################ DESCRIPTION
#
# Opens up ROOT file
# Removes baseline and trendline from events
# Applies butterworth filter to events
# Applies rolling mean to Events
# Collects properties of events to determine signals
# Shows all signals identified within the list (defined by y)
################


############### FILES
# PMTsignals/Run203-PMT107.root
# PMTsignals/Run203-PMT78.root
# PMTsignals/Run103-noise-PMT78.root
# PMTsignals/Run103-noise-PMT107.root


# Open the data, apply to variable
file = "E:\PMTsignals\Boulby_78.root"

tree = uproot.open(file)["Tree"]
branches = tree.arrays()
print(branches['ADC'])

# Move data to less nefarious sounding variable
datafull = branches['ADC']

# how many values in each event
samples = 150
# how long between samples within events
timegate = 2

# Time axis for
time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(samples):
    time.append(i*2)

# Event control, how many events do you want to process?
y = 1000000

# take the y values from data, to stop array index mismatching
data = datafull[:y]


# Find mean OBSOLETE, ONLY USE IF USING SIGMAEVENTS CODE
#meanvals = []
#for i in range(y):
#    datarepeated = data[i]
    # Append data values to list
#    meanvals.append(np.mean(datarepeated))

# timing each process
t0 = process_time()


# Apply LCMS baseline & trendline removal
scaleddata = fnc.LCMSlist(data)

# timing each process
t1 = process_time()
totalLCMS = t1-t0
print("[1/5] LCMS algorithm Applied. Runtime: " + str(totalLCMS) + "s")



# Butterworth Filter, removing everything above 2000MHz

# timing each process
t2 = process_time()

#filterdata = []
sos = signal.butter(5,200,'lp',fs=500,output='sos')
# Applying filter to data
for i in range(y):
    newdata = signal.sosfilt(sos,scaleddata[i])
    scaleddata[i] = newdata
    #filterdata.append(signal.sosfilt(sos,scaleddata[i]))

# timing each process
t3 = process_time()
totalbutterworth = t3-t2
print("[2/5] Butterworth Filter Applied. Runtime: " + str(totalbutterworth) + "s")

# Apply rolling mean to data, window of 5 samples currently

# timing each process
t4 = process_time()
#rollingdata = []
for i in range(len(scaleddata)):
    newerdata = fnc.rollmean(scaleddata[i],5)
    scaleddata[i] = newerdata
    # rollingdata.append(fnc.rollmean(scaleddata[i],5))

# timing each process
t5 = process_time()
totalrolling = t5-t4
print("[3/5] Rolling Mean Applied. Runtime: " + str(totalrolling) + "s")

# Collect new average values for signal properties calculations
#fmeanvals, fstdvals, fminvals, fmaxvals, fsigvals, fmedvals = fnc.datacollate(rollingdata, len(rollingdata))

# timing each process
t6 = process_time()

# To increase speed, moved this such that it is a much shorter and more concise piece of code
fminvals = []
fmedvals = []

# Collecting minimum and median for property evaluation
for o in range(y):
    datarepeater = scaleddata[o]
    # Append data values to list
    fminvals.append(np.amin(datarepeater))
    fmedvals.append(statistics.median(datarepeater))




# Finding properties of events

# sigma values
####sigmaval = 0
# lists
fwhmlst = []
sgdpthlst = []


for n in range(len(scaleddata)):


    ####onesigvals = fnc.sigmaevents(rollingdata[n],fmedvals[n],fstdvals[n],sigmaval) # 0 -> fmeanvals[n]

    # Takes difference from median to minimum value to find maximum of spike
    difference = fmedvals[n] + fminvals[n]

    # find half minimum (which is half maximum for us)
    halfmin = difference / 2
    # find nearest point within the array to half max(min). RETURNS THE ARRAY ELEMENT NUMBER! NOT THE ARRAY VALUE!
    nearest = (np.abs(scaleddata[n]-halfmin)).argmin()
    # Once spotted, start counting
    # count how many events occur between nearest and minimum
    spotted = 0
    count = 0
    for k in range(len(scaleddata[n])):
        # whichever component it comes into contact first, start the count
        if (scaleddata[n][k] == fminvals[n]) or (scaleddata[n][k] == scaleddata[n][nearest]):
            spotted += 1

        # if one component has been spotted, but not the other, keep counting
        if spotted == 1:
            count += 1

        # if both components have been spotted
        if spotted == 2:
            break

    # Apply length of sample size
    fwhmlength = 2*count*timegate
    #print("Event: " + str(n) + "\nFWHM Value: " + str(fwhmlength))
    # plotting purposes
    #fwhmlist = [FWHMVAL] * len(time)
    fwhmlist = [halfmin] * len(time)

    # to look at individual values for things ################## FWHM SIGNAL VIEWER ##################
    #check = 388
    #if n == check:
    #    print("Event " + str(388))
    #    print("halfmin, rollingdata["+str(check)+"][nearest] , fminvals["+str(check)+"], fwhmlength")
    #    print(halfmin, rollingdata[n][nearest], fminvals[check], fwhmlength)
    #    plt.plot(time,rollingdata[n], label="Normal")
    #    plt.plot(time,onesigvals, label=str(sigmaval) + " sigma")
    #    plt.plot(time,fwhmlist, label="FWHM")
    #    plt.legend()
    #    plt.title("Butterworth Rolling Filtered Signal with FWHM of " + str(file) + " event " + str(n))
    #    plt.xlabel("Sample Time (ns)")
    #    plt.ylabel("ADC Value")
    #    plt.show()

    # Apply to lists
    fwhmlst.append(fwhmlength)
    sgdpthlst.append(fminvals[n])

# timing each process
t7 = process_time()
totalproperties = t7-t6
print("[4/5] Properties Determined. Runtime: " + str(totalproperties) + "s")

# Construct data cutoff point here, remove all values with signal depth > -50 and see how many signals are left
cutdata = []
cutdepth = []
# signal depth cut off
sgdpthcutoff = -30
# fwhm cut offs
fwhmlowerbound = 7.5
fwhmupperbound = 20.5


for i1 in range(len(scaleddata)):
    if (sgdpthlst[i1] < sgdpthcutoff): #and (fwhmlst[i1] < fwhmupperbound) and (fwhmlst[i1] > fwhmlowerbound):
        #cutdepth.append(sgdpthlst[i1])
        cutdata.append(i1)

print("[5/5] Dark Counts Determined")

# Depth Cut off
#plt.hist(cutdepth,bins = 25)
#plt.title("Signal Depth Histogram for file " + str(file) + " with Cut-off at " + str(sgdpthcutoff))
#plt.ylabel("Count")
#plt.xlabel("Depth (ADC Value)")
#plt.show()

print("Number of signal events detected:")
print(len(cutdata))
