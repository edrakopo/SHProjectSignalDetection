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
# Opens up ROOT file, roughly highlights signal events
# Removes baseline from signal events (complicated)
# Calculates FT for baseline removed events#
# Possibly superimpose them
################


############### FILES
# PMTsignals/Run203-PMT107.root
# PMTsignals/Run203-PMT78.root
# PMTsignals/Run103-noise-PMT78.root
# PMTsignals/Run103-noise-PMT107.root

# Open the data, apply to variable
file = "PMTsignals/Run203-PMT78.root"

tree = uproot.open(file)["Tree"]
branches = tree.arrays()
print(branches['ADC'])

# Move data to less nefarious sounding variable
data = branches['ADC']

# how many values in each event
samples = 150


# Time axis for
time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(samples):
    time.append(i*2)

# Event control, how many events do you want to process?
y = 10





# First, find signal events, use basic 5sigma method initially, can be changed at later date

# using basic mean/sig finder
meanvals, stdvals, minvals, maxvals, sigvals = fnc.datacollate(branches['ADC'], y)

print(sigvals)

# Two lists, one for holding signal values, one for holding event number
sigdata = []
sigevents = []


# Collect signal values
for k in range(y):
    if sigvals[k] == 1:
        sigdata.append(data[k])
        sigevents.append(k)
    else:
        print("No signal")

print("Finding plots")

# Baseline subtraction
#i=0
#scaleddata=[]
#for i in range(len(sigdata)):
    # remove baseline
    #print("event " + str(i) + ": \n" + str(branches['ADC'][i]))
    #print("subtract " + str(c[i]))
#    scaleddata.append(sigdata[i]-meanvals[i])
    #print("baseline removed event" + str(i) + ": \n" + str(scaleddata[i]))

# baseline subtraction
scaleddata = fnc.baselinesubtraction(sigdata,meanvals)

#plt.plot(time,scaleddata[0])
#plt.xlabel("Sample Time (ns)")
#plt.ylabel("ADC Value")
#plt.title(str(file) + " event " + str(sigevents[0]))
#plt.show()



# Fourier transform of signal

# Creating variables required for composite
#dftdata = [None] * len(sigdata)
#freqdata = pyfftw.interfaces.numpy_fft.rfftfreq(len(time),1/500)
dfthist = [0] * (samples//2+1)

dftdata,freqdata = fnc.fouriertransformsimple(time, scaleddata, 500, False)

# create FT of all data, and compilation
for q in range(len(sigdata)):
    #abs to remove negative components
    #dftdata[q] = np.abs(pyfftw.interfaces.numpy_fft.rfft(scaleddata[q]))
    dfthist = np.add(dfthist,dftdata[q])

# plot fourier values, commented out
#plt.plot(freqdata,dftdata[0])
#plt.xlabel("Frequency (MHz)")
#plt.ylabel("Amplitude")
#plt.title(str(file) + "FT of event " + str(sigevents[0]))
#plt.show()

# plot fourier values compilation, commented out
plt.plot(freqdata,dfthist)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Amplitude")
plt.title(str(file) + "FT compilation of event ")
plt.show()


## DISTRIBUTION OF 10TH EVENT WITH SIGNAL
# Create distribution of 10th event - Trendline/baseline removed
# filter out signal values (Nonetype right now) from eventdistr
#lindatasig = []
#for i in range(0,len(scaleddata)):
#    lindatasig.append(scaleddata[i])
#eventdistr = fnc.adcdist(lindatasig,len(lindatasig))






# looking at butterworth filter, removing everything above 200MHz
sos = signal.butter(5,200,'lp',fs=500,output='sos')

# FT of butterworth data
filterdata = []
#filtdftdata = [None] * len(sigdata)

for i in range(len(sigdata)):
    filterdata.append(signal.sosfilt(sos,scaleddata[i]))
    #filtdftdata[i] = np.abs(pyfftw.interfaces.numpy_fft.rfft(filterdata[i]))

# FT of butterworth filtered data
filtdftdata = fnc.fouriertransformsimple(time, filterdata, 500, True)

## DISTRIBUTION OF 10TH EVENT WITH SIGNAL - BUTTERWORTH
# Create distribution of 10th event - Trendline/baseline removed
# filter out signal values (Nonetype right now) from eventdistr
#butsig = []
#for i in range(0,len(filterdata)):
#    butsig.append(filterdata[i])
#eventdistr = fnc.adcdist(butsig,len(butsig))


# Butterworth Signal Plotting
plt.plot(time,filterdata[0])
plt.title("Butterworth Filtered Signal of " + str(file) + " event " + str(sigevents[0]))
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.show()

# Fourier Transform Signal Plotting
plt.plot(freqdata,filtdftdata[0])
plt.title("Butterworth Fourier Transform of " + str(file) + " event " + str(sigevents[0]))
plt.xlabel("Frequency (MHz)")
plt.ylabel("Amplitude")
plt.show()

# Study components of first signal event (sigevents[0]) will expand to all signal events soon
# Most likely will be split into functions in functions.py


# signal length CODE
    # find the point where the value goes 1sig/2sig/3sig from mean and plot over normal data, figure out what values work!
    # Then take the amount of time samples from the start to the end

# signal depth CODE


# integrated charge CODE
    # from time events, take how much the signal deviated from the mean additively

# rise time CODE
    # time to go from 0.1 to 0.9 of signal amplitude (SOMETHING ELSE NEEDED TO BE CALCULATED)
