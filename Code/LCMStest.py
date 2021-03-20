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

"""
Made to test LCMS function within functions.py
"""



# Open the data, apply to variable
file = "E:\PMTsignals\Boulby_78_Signal.root"

tree = uproot.open(file)["Tree"]
branches = tree.arrays()
print(branches['ADC'])

# Move data to less nefarious sounding variable
datafull = branches['ADC']

# how many values in each event, taken from first event
samples = len(datafull[0])
# how long between samples within events
timegate = 2

# Time axis for
time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(samples):
    time.append(i*2)

# take the y values from data, to stop array index mismatching

# select event 9
dataint = datafull[9]

#plt.plot(time,dataint)
#plt.xlabel("Sample Time (ns)")
#plt.ylabel("ADC Value")
#plt.title(str(file) + " event " + str(9))
#plt.show()

# LCMS usage
newdataint = fnc.LCMS(dataint)

#plt.plot(time,newdataint)
#plt.xlabel("Sample Time (ns)")
#plt.ylabel("ADC Value")
#plt.title(str(file) + " event " + str(9) + " LCMS scaled")
#plt.show()


# Testing linfit against LCMS

# Limiting size
y = 1000000
datanew = datafull[:y]




t0 = process_time()

# LINFIT CODE

# calculate mean
meani = []
for i in range(len(datanew)):
    meani.append(np.mean(datanew[i]))

# Baseline subtraction
scaleddata = fnc.baselinesubtraction(datanew,meani)

# Linear Fitting
pfit, stats, rms, c, m = fnc.linfit(time,scaleddata,y)

# Apply linear fit over data
lindatanew = []
yline = []
for i in range(y):
    for j in range(len(time)):
        # Write new data to dummy list
        lindatanew.append(scaleddata[i][j] - m[i]*time[j])

    yline.append(lindatanew)
    lindatanew=[]


t1 = process_time()

totalLINFIT = t1-t0






t2 = process_time()

# LCMS CODE
datanewLCMS = []
for i in range(len(datanew)):
    datanewLCMS.append(fnc.LCMS(datanew[i]))

t3 = process_time()

totalLCMS = t3 - t2





# LCMSfast code

t4 = process_time()
# Collect x values respective to number of y components
timer = []
for p in range(len(datanew[0])):
    timer.append(p)

# take division of time component
xhalf = np.divide(np.max(timer),2)
#print("xhalf: " + str(xhalf))

# Create new basis for time
xi = []
for j in range(len(timer)):
    xi.append(j - xhalf)
#print("xi: " + str(xi))

# Main code
datanewLCMSfast = []
for j in range(len(datanew)):
    datanewLCMSfast.append(fnc.LCMSfast(datanew[j],xhalf,xi))


t5 = process_time()

totalLCMSfast = t5 - t4

# LCMSlist code
t6 = process_time()

datanewLCMSlist = []
datanewLCMSlist = fnc.LCMSlist(datanew)

t7 = process_time()

totalLCMSlist = t7 - t6


print("Time taken for LINFIT: " + str(totalLINFIT) + "s")
print("Time taken for LCMS: " + str(totalLCMS) + "s")
print("Time taken for LCMSfast: " + str(totalLCMSfast) + "s")
print("Time taken for LCMSlist: " + str(totalLCMSlist) + "s")

# Plot LINFIT AND LCMSdata together
#plt.plot(time,yline[9],color='black',markersize=2, label = "Polynomial.fit (LS)", linewidth =2)
#plt.plot(time,datanewLCMS[9], color = 'red', markersize = 2, label = 'LCMS', linewidth = 3)
#plt.plot(time,datanewLCMSfast[9], color = 'yellow', markersize = 2, label = 'LCMSfast', linewidth = 0.5)
#plt.plot(time,datanewLCMSlist[9], color = 'blue', markersize = 2, label = 'LCMSlist', linewidth = 0.5)
#plt.legend(fontsize = 15)
#plt.xlabel("Sample Time (ns)", fontsize = 17)
#plt.ylabel("ADC Value", fontsize = 17)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.title("Least Squares Trendline Removal for event 9", fontsize = 22)
#plt.show()
