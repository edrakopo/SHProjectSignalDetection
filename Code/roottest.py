import uproot
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc
import pandas as pd


from scipy import signal
from scipy.fft import rfft, rfftfreq

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

# 10000 events
data = branches['ADC'][10]

time = []
# Creating list for sample times that are 2ns intervals, 150 samples
for i in range(150):
    time.append(i*2)




for i in range(50):
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


# Choose rolling average?
join = input("Use running mean data? [Y/N]")
if join.lower() == 'y':
    # Testing rolling averages with pandas
    df = pd.DataFrame({'A': data})
    # window of 3 for our rolling average
    df_rolling = df.rolling(5,min_periods=1).mean()
    print(df_rolling)
    data = df_rolling['A'].tolist()
elif join.lower() == 'n':
    print("Continuing without running mean data...")



# ACHTUNG, THIS CODE IS NOW OBSOLETE AND WILL MOST LIKELY NOT FUNCTION
# Line of best fit, to see modulation in baseline
yline, pfit, stats, rms = fnc.linfit(time,data,1)
c, m = pfit



# Plot line of best fit over data
plt.plot(time,data,color='black',markersize=2)
plt.plot(time,yline,color='red',linewidth=3)
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.title("Trendline of " + str(file))
plt.show()




# Now, remove baseline (and possibly also linear modulation, if it effects results significantly)
bslndata = data - c




# Plot removed baseline
plt.plot(time,bslndata,color='black',markersize=2)
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.title("Trendline of " + str(file) + " with baseline reduction")
plt.show()


# FFT experimentation
# will alter sample spacing if it causes issues

#sampling rate is 500MHz (2ns sampling time over 300ns)

yf = rfft(bslndata)
xf = rfftfreq(len(bslndata),1/500)

plt.plot(xf,yf)
plt.title("Rough Fourier Transform of " + str(file))
plt.xlabel("Frequency (MHz)")
plt.ylabel("Amplitude")
plt.show()


# looking at butterworth filter, removing everything above 220MHz

sos = signal.butter(10,200,'lp',fs=500,output='sos')
filtereddata = signal.sosfilt(sos,bslndata)


plt.plot(time,filtereddata)
plt.title("Butterworth Filtered Signal of " + str(file))
plt.xlabel("Sample Time (ns)")
plt.ylabel("ADC Value")
plt.show()

# Fourier Transform
yfbut = rfft(filtereddata)

plt.plot(xf,yfbut)
plt.title("Butterworth Fourier Transform of " + str(file))
plt.xlabel("Frequency (MHz)")
plt.ylabel("Amplitude")
plt.show()
