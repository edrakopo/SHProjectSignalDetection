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

# For moving to excel
import xlsxwriter

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
file = "E:\PMTsignals\Run103-noise-PMT78.root"

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
y = 100000

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
print("[1/4] LCMS algorithm Applied. Runtime: " + str(totalLCMS) + "s")



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
print("[2/4] Butterworth Filter Applied. Runtime: " + str(totalbutterworth) + "s")

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
print("[3/4] Rolling Mean Applied. Runtime: " + str(totalrolling) + "s")

# Collect new average values for signal properties calculations
#fmeanvals, fstdvals, fminvals, fmaxvals, fsigvals, fmedvals = fnc.datacollate(rollingdata, len(rollingdata))

# timing each process
t6 = process_time()


# Finding signal event no.s and the events themselves
cutdata, signalevents = fnc.signalspotter(scaleddata,timegate,-30,(7.5,20.5))

t7 = process_time()
totalproperties = t7-t6
print("[4/4] Dark Counts Determined. Runtime: " + str(totalproperties) + "s")

# Depth Cut off
#plt.hist(cutdepth,bins = 25)
#plt.title("Signal Depth Histogram for file " + str(file) + " with Cut-off at " + str(sgdpthcutoff))
#plt.ylabel("Count")
#plt.xlabel("Depth (ADC Value)")
#plt.show()

#print(cutdata[:200])

# Write cutdata to excel array, completely disposable. Just allows for easier transferring
# When calculating efficiencies
#workbook = xlsxwriter.Workbook('cutdata.xlsx')
#worksheet = workbook.add_worksheet()

#worksheet.write_column(0,0,cutdata[:200])
#workbook.close()
