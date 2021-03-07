import matplotlib.pyplot as plt
import numpy as np
import math as m
import pyfftw
import pandas as pd
import statistics
import uproot

from scipy import signal
from scipy import signal
from scipy.fft import rfft, rfftfreq
from scipy import stats


# For writing pydoc
# py -m pydoc -w functions

# File Opening
def rootopen(file):
    """
    Opens a root file and pulls the appropriate data into an array, for easier use.
    data is now an array of all the events, in which sample data is held.
    Returns the event data.
    -> file
        -> Tree
            -> branches
                -> ['ADC']
                    -> [Events 0...n]

    Example usage;
    datafull = fnc.rootopen("E:\PMTsignals\Run203-PMT78.root")

    :param root file:
    :return data:
    """
    tree = uproot.open(file)["Tree"]
    branches = tree.arrays()


    # Move data to less nefarious sounding variable
    datafull = branches['ADC']

    return(datafull)


# Linear Fit
def linfit(x,y,z):
    """
    Method for finding the linear fit of a 2D plot for multiple non-signal events.
    Returns the numpy array of the y axis points for the best fit line.
    WARNING - Only applies to branches['ADC'][i], as hardcoded into 'data' variable.
              If you want this to change, you can alter it within the code on line 33

    OBSOLETE - This and baselinesubtraction have both been replaced by function LCMS

    :param x:           x values
    :param y:           y values
    :parad z:           number of events found
    :return pfit:       fitting parameters(offset value, gradient value)
    :return stats:      stats(residuals, rank, singular values of matrix, rcond)
    :return rms:        root mean square residual value
    :return c:          offset
    :return m:          mean
    """
    # time and anything related can be declared outwith the loop, it is unchanging
    time = x
    cmin, cmax = min(time), max(time)
    # Arrays generated, [None] * z method used for most just due to funkiness with appending a polynomial fit
    # Will return to this time permitting
    pfit = [None] * z
    stats = [None] * z
    rms = []
    c = [None] * z
    m = [None] * z

    for i in range(z):
        # Find equation of straight line through data to determine baseline & baseline modulation

        data = y[i]
        # https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyfit.html
        Polynomial = np.polynomial.Polynomial
        # Fit to range cmin, cmax
        pfit[i], stats[i] = (Polynomial.fit(time, data, 1, full=True, window=(cmin, cmax), domain=(cmin, cmax)))

        #print("Fitting Results: ", pfit[i], stats[i], sep='\n')

        # labelling variables from fitting for ease of use
        c[i], m[i] = pfit[i]
        # Not list based as they are not saved
        resid, rank, sing_val, rcond = stats[i]
        rms.append(np.sqrt(resid[0]/len(data)))

        #print("Fit: y = " + str(m[i]) + "x + " + str(c[i]))
        #print("RMS residual = " + str(rms[i]))



    return(pfit, stats, rms, c, m)

# Collation of different components of data
def datacollate(branches, y):
    """
    Method for collating important characeristics of each event within lists.
    Returns the means, standard deviations, minimum values, maximum values, and possible signal values of all events within (y) events registered

    :param branches:    Event data from ROOT file
    :param y:           number of events to be analysed
    :return meanvals:   Mean values of all events in a list
    :return stdvals:    Standard deviation values of all events in a list
    :return minvals:    Minimum values of all events in a list
    :return maxvals:    Maxmimum values of all events in a list
    :return sigvals:    List marking all events that had possible signals within them (0 or 1)
    :return medianvals: Median values of all events in a list
    """
    # Mean,std, min, max, data collection

    meanvals = []
    stdvals = []
    minvals = []
    maxvals = []
    sigvals = []
    medvals = []

    # Testing range of y values
    for i in range(y):

        datarepeated = branches[i]
        # Append data values to list
        meanvals.append(np.mean(datarepeated))
        stdvals.append(np.std(datarepeated))
        minvals.append(np.amin(datarepeated))
        maxvals.append(np.amax(datarepeated))
        medvals.append(statistics.median(datarepeated))


        # Create array holding values 0 or 1 for all possible signals
        # If mean - min val is larger than 5 sigma, mark as signal
        # If max val - mean is larger than 5 sigma, mark as signal
        if (meanvals[i] - minvals[i]) > (5*stdvals[i]) or (maxvals[i] - meanvals[i]) > (5*stdvals[i]):
            # Signal
            sigvals.append(1)

        else:
            # No signal
            sigvals.append(0)

    return(meanvals, stdvals, minvals, maxvals, sigvals, medvals)


# signal data for x sigms
def sigmaevents(data, meanval, stdval, x):
    """
    Takes data from one event and checks what values are beyong X sigma of the mean.
    Returns list of these Values, with all values not beyond sigma set at mean.

    :param data:        Data from one signal event
    :param meanvals:    Mean valuation of said event
    :param stdvals:     Standard deviation of event
    :param x:           The number of sigma needed to trigger listing
    :return sigdata:    List of values beyond x sigma
    """
    # Events beyond x sigma
    sigmasamples = []
    # Test variables, uncomment if needed
    #print("Data 0 "+ str(data[0]))
    #print("Meanval "+ str(meanval))
    #print("Stdval "+ str(stdval))
    #print("x " + str(x))

    for i in range(len(data)):
        if (abs(meanval - data[i])) > (x*stdval):
            sigmasamples.append(data[i])
        else:
            sigmasamples.append(meanval)

    return(sigmasamples)

# distribution of ADC Values
def adcdist(events,samplesize):
    """
    Takes all events, and records the 10th sample from each event.
    Then creating histogram of the samples, and returns the collection of 10th events.

    :param events:          Events from ROOT file
    :param samplesize:      No. of samples, to scale for different sample lengths
    :return PLOT:           Plot of histogram distribution
    :return eventdistr:     List of all 10th samples from all events
    """

    eventdistr = []
    for i in range(0,len(events)):
        # take the 10th sample of every event and add to list
        eventdistr.append(events[i][9])

    # Plotting histogram
    plt.hist(eventdistr,100)
    plt.title("Histogram of 10th value from all events")
    plt.show()

    return(eventdistr)


# baseline subtraction function
def baselinesubtraction(data,mean):
    """
    Takes data and mean value, and removes the baseline from all data passed in.
    Returns array of data with baseline removed

    OBSOLETE - This and LINFIT have both been replaced by function LCMS

    :param data:        Array of event data
    :param mean:        Array of all mean across all Events
    :return scaleddata: Array of event data baseline removed
    """
    scaleddata=[]
    for i in range(len(data)):
        # remove baseline
        #print("event " + str(i) + ": \n" + str(branches['ADC'][i]))
        #print("subtract " + str(c[i]))
        scaleddata.append(data[i]-mean[i])
        #print("baseline removed event" + str(i) + ": \n" + str(scaleddata[i]))
    return(scaleddata)


# Condensed FT creation
def fouriertransformsimple(time, scaleddata, samplefreq, freq):
    """
    Takes time, data, and the sample frequency. Asks if frequency is needed.
    Returns FT and FT-frequency of said data.

    :param time:        Array of all time Values
    :param scaleddata:  Array of events
    :param samplefreq:  Sampling frequency
    :param freq:        Boolean for whether or not freq needs to be calculated
                        True - Calculated already
                        False - Not yet calculated
    :return dftdata:    FT data of events
    :return freqdata:   FT-frequency data
    """
    # Creating variables required
    dftdata = [None] * len(scaleddata)
    # Create frequency data is needed (False)
    if freq == False:
        freqdata = pyfftw.interfaces.numpy_fft.rfftfreq(len(time),1/samplefreq)

    # Create FT data
    for q in range(len(scaleddata)):
        #abs to remove negative components
        dftdata[q] = np.abs(pyfftw.interfaces.numpy_fft.rfft(scaleddata[q]))
    # return
    if freq == False:
        return(dftdata,freqdata)
    else:
        return(dftdata)


def rollmean(data,window):
    """
    Takes data, and the window size.
    Returns the rolling mean data

    :param data:    Array of data Values
    :param window:  Size of window
    :return:        Array of rolling mean values
    """
    # Testing rolling averages with pandas
    df = pd.DataFrame({'A': data})
    # window of 3 for our rolling average
    df_rolling = df.rolling(window,min_periods=1).mean()
    rollingdata = df_rolling['A'].tolist()
    return(rollingdata)



def LCMS(data):
    """
    Use LCMS Algorithm to remove linear trendlines within data.
    Applies to one event, iterate within a list to apply to all events.

    :param data:        List of sample data (y component)
    :return:            Array of data with linear trendline removed
    """
    # Center time around 0, such that average is 0
    #time -= np.divide(np.max(time),2)
    # example of this, 300ns gate
    # -150 from all components, 0->-150, 300->150

    # Collect sum squared time
    #timesqsum = np.sum(np.square(time))

    ###### LCMS METHOD #####

    # Collect x values respective to number of y components
    time = []
    for p in range(len(data)):
        time.append(p)

    # First iteration ###################

    # LCMS time centering (which for some reason isn't centered, even though it could be using the above method)

    # All components outwith the w loop are constant throughout the function, due to their relation to x component

    # take division of time component
    xhalf = np.divide(np.max(time),2)
    #print("xhalf: " + str(xhalf))

    # Create new basis for time
    xi = []
    for j in range(len(time)):
        xi.append(j - xhalf)
    #print("xi: " + str(xi))

    # Collect xi^2
    xisq = np.sum(np.square(xi))
    #print("xisq: " + str(xisq))


    for w in range(1):

        # Baseline subtraction
        y = np.mean(data)
        a = data - y

        # Collect summed slope component
        slope = 0
        for i in range(len(a)):
            slope += (i - xhalf)*a[i]

        # Slope Calculation
        s = (1/xisq)*slope

        # apply new common-mode removed data to be reused in second loop
        data = []
        for k in range(len(a)):
            # Different CM subtractions based on whether first or second iteration
            if w == 0: # First iteration
                data.append(a[k]-s*(k-xhalf))
            if w == 1: # Second iteration
                data.append(a[k]-s*(k-xhalf)-y)

    # Return new data
    return(data)


def LCMSfast(data,xhalf,xi):
    """
    Use LCMS Algorithm to remove linear trendlines within data.
    Applies to one event, iterate within a list to apply to all events.

    Faster method, due to the x axis appropriate for this function
    being determined outwith the function. Please look at
    https://inspirehep.net/literature/928989
    more guidance on this issue, where xi is xi and xhalf is
    equivalent to 16 from section 2.2.

    :param data:        List of sample data (y component)
    :param xhalf:       Half the largest value in the x components
    :param xi:          New x axis based on xhalf, in which each point is one entry on the y axis
    :return datafn:     Array of data with linear trendline removed
    """
    # Center time around 0, such that average is 0
    #time -= np.divide(np.max(time),2)
    # example of this, 300ns gate
    # -150 from all components, 0->-150, 300->150

    # Collect sum squared time
    #timesqsum = np.sum(np.square(time))

    ###### LCMS METHOD #####

    # Collect xi^2
    xisq = np.sum(np.square(xi))
    #print("xisq: " + str(xisq))


    for w in range(1):

        # Baseline subtraction
        y = np.mean(data)
        a = data - y

        # Collect summed slope component
        slope = 0
        for i in range(len(a)):
            slope += (i - xhalf)*a[i]

        # Slope Calculation
        s = (1/xisq)*slope

        # apply new common-mode removed data to be reused in second loop
        data = []
        for k in range(len(a)):
            # Different CM subtractions based on whether first or second iteration
            if w == 0: # First iteration
                data.append(a[k]-s*(k-xhalf))
            if w == 1: # Second iteration
                data.append(a[k]-s*(k-xhalf)-y)

    # Return new data
    return(data)


def LCMSlist(datas):
    """
    Use LCMS Algorithm to remove linear trendlines within data.
    Applies to a list of events. Fastest method out of all LCMS methods

    :param data:        Array of events
    :return datafn:     Array of events with linear trendline removed
    """

    ###### LCMS METHOD #####
    # Abstraction from argument to allow for modification
    dataz = datas

    # All components outwith the w loop are constant throughout the function, due to their relation to x component

    # Collect x values respective to number of y components
    time = []
    for p in range(len(dataz[0])):
        time.append(p)

    # take division of time component
    xhalf = np.divide(np.max(time),2)
    #print("xhalf: " + str(xhalf))

    # Create new basis for time
    xi = []
    for j in range(len(time)):
        xi.append(j - xhalf)
    #print("xi: " + str(xi))

    # Collect xi^2
    xisq = np.sum(np.square(xi))
    #print("xisq: " + str(xisq))


    # loop across multiple events and apply to datafn
    datafn = []
    for o in range(len(dataz)):

        data = dataz[o]

        for w in range(1):

            # Baseline subtraction
            y = np.mean(data)
            a = data - y

            # Collect summed slope component
            slope = 0
            for i in range(len(a)):
                slope += (i - xhalf)*a[i]

            # Slope Calculation
            s = (1/xisq)*slope

            # apply new common-mode removed data to be reused in second loop
            data = []
            for k in range(len(a)):
                # Different CM subtractions based on whether first or second iteration
                if w == 0: # First iteration
                    data.append(a[k]-s*(k-xhalf))
                if w == 1: # Second iteration
                    data.append(a[k]-s*(k-xhalf)-y)

        if (o%(len(dataz)/50)==0):
            print(o)
        # append list
        datafn.append(data)


    # Return new data
    return(datafn)

def signalspotter(data, timegate, depthcutoff, fwhmcutoff):
    """
    Collects data given, and the cutoff values required.
    Calculates the FWHM and depth of all events.
    Determines if events are signals based on cutoff values.
    Returns the number of signals spotted and a list of signals events.

    :param data:                    Array of events
    :param timegate:                Length of time between each recorded sample (eg: 2ns -> 2)
    :param depthcutoff:             Value for signal depth cutoff
    :param fwhmcutoff:              Upper and Lower bounds for FWHM cutoff, written as (LB,UB)
    :return signalno, signallist:   List of the signal event numbers, as well as the signal events themselves
    """
    # Determine event list length
    eventno = len(data)

    # Declare minimum and median values for all events
    fminvals = []
    fmedvals = []

    # Using timegate, and event length. Create array of time variables



    # Collecting minimum and median for property evaluation
    for o in range(eventno):
        # select individual data points
        datarepeater = data[o]
        # Append data values to list
        fminvals.append(np.amin(datarepeater))
        fmedvals.append(statistics.median(datarepeater))


    # Finding properties of events

    # lists
    fwhmlst = []
    sgdpthlst = []


    for n in range(eventno):

        # FWHM FINDING

        # Takes difference from median to minimum value to find maximum of spike
        # Positive because fminvals is always negative, and median is basically 0 in every case
        # This could be altered to be better.
        difference = fmedvals[n] + fminvals[n]

        # find half minimum (which is half maximum for us)
        halfmin = difference / 2
        # find nearest point within the array to half max(min). RETURNS THE ARRAY ELEMENT NUMBER! NOT THE ARRAY VALUE!
        nearest = (np.abs(data[n]-halfmin)).argmin()
        # Once spotted, start counting
        # count how many events occur between nearest and minimum
        spotted = 0
        count = 0
        for k in range(eventno):
            # whichever component it comes into contact first, start the count
            if (data[n][k] == fminvals[n]) or (data[n][k] == data[n][nearest]):
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
        fwhmlist = [halfmin] * len(data[0])

        # Signal depth finding is very simple, as it requires just the lowest value within the event (not 100% true, but true for our linearised, smoothed data)

        # Apply property to lists
        fwhmlst.append(fwhmlength)
        # [n] applied here as fminvals is already an array that covers all events
        sgdpthlst.append(fminvals[n])

    # Signal spotting

    # Collating important variables before runtime
    cutdata = []
    datanumber = []
    sgdpthcutoff = depthcutoff
    fwhmlowerbound, fwhmupperbound = fwhmcutoff

    # If event passes the cutoff, select the event in its entirety as well as the number of the event
    for i1 in range(eventno):
        if (sgdpthlst[i1] < sgdpthcutoff) and (fwhmlst[i1] < fwhmupperbound) and (fwhmlst[i1] > fwhmlowerbound):
            cutdata.append(data[i1])
            datanumber.append(i1)

    print("Number of signal events detected")
    print(len(datanumber))

    # Return the number of each signal event, as well as all the data from each event
    return(datanumber, cutdata)


def scaledata(data, butterfreq, rollingwindow, lcms, butter, rolling):
    """
    Collects the standard signal data, and applies the different effects to the events
    to improve the ability to identify key components.

    Effects are:
    -LCMS pedestal and trendline removal
    -Butterworth low-pass frequency filter
    -Rolling mean across event data

    Includes modifiable butterworth frequency value, and rolling mean window size.

    :param data:            Array of events
    :param butterfreq:      Integer for the frequency limit of the butterworth filter (MHz)
    :param rollingwindow:   Integer for size of rolling mean window
    :param lcms:            Bool for whether LCMS is applied
    :param butter:          Bool for whether butterworth filter is applied
    :param rolling:         Bool for whether rolling mean is applied

    :return scaleddata:    Data with corresponding effects applied appropriately
    """

    # Length of data sample, used for loops
    y = len(data)

    # LCMS APPLICATION
    if (lcms == True):
        # Cost of forming a new full list is the memory issues that come with it.
        # The positive is that it runs faster than reapplying the function to all events over a list.
        scaleddata = LCMSlist(data)
        print("LCMS Algorithm applied...")
    else:
        print("Skipping LCMS Algorithm...")

    # BUTTERWORTH APPLICATION
    if (butter == True):
        # Creating the butterworth filter
        sos = signal.butter(5,butterfreq,'lp',fs=500,output='sos')
        # Applying filter to data
        for i in range(y):
            newdata = signal.sosfilt(sos,scaleddata[i])
            scaleddata[i] = newdata
        print("Butterworth Filter Applied at " + str(butterfreq) + " MHz...")
    else:
        print("Skipping Butterworth Filter...")

    # ROLLING MEAN APPLICATION
    if (rolling == True):
        for i in range(len(scaleddata)):
            newerdata = rollmean(scaleddata[i],rollingwindow)
            scaleddata[i] = newerdata
        print("Rolling Mean Applied...")
    else:
        print("Skipping Rolling Mean...")

    # If no modification to our data has been applied, write exception
    try:
        return(scaleddata)
    except:
        raise ValueError("No modifications were chosen. Please choose one modification to apply to your data.")
