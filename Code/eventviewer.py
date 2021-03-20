import uproot
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc

# Viewing particular event of particular file.
# Open the data, apply to variable
file = "E:\PMTsignals\Run203-PMT78.root"
# File directory specific, just to remove the E:\PMTsignals\ rubbish
filename = file[14:]


tree = uproot.open(file)["Tree"]
branches = tree.arrays()
#print(branches['ADC'])

# how long between data taken
timegate = 2
# length of event
eventno = len(branches['ADC'][0])
time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(eventno):
    time.append(i*timegate)

# Input event
Nevent = int(input("What event do you wish to view? "))

plt.plot(time,branches['ADC'][Nevent])
plt.xlabel("Sample Time (ns)", fontsize = 17)
plt.ylabel("ADC Value", fontsize = 17)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title(str(filename) + " event " + str(Nevent), fontsize = 22)
plt.show()
