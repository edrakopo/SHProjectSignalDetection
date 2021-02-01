import matplotlib.pyplot as plt
import numpy as np
import math as m
import pyfftw

# Linear Fit
def linfit(x,y,z):
    """
    Method for finding the linear fit of a 2D plot for multiple non-signal events.
    Returns the numpy array of the y axis points for the best fit line.
    WARNING - Only applies to branches['ADC'][i], as hardcoded into 'data' variable.
              If you want this to change, you can alter it within the code on line 33

    :param x: x values
    :param y: y values
    :parad z: number of events found
    :return: fitting parameters(offset value, gradient value), stats(residuals, rank, singular values of matrix, rcond), root mean square residual value
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

        data = y['ADC'][i]
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
    """
    # Mean,std, min, max, data collection

    meanvals = []
    stdvals = []
    minvals = []
    maxvals = []
    sigvals = []

    # Testing range of y values
    for i in range(y):
        datarepeated = branches[i]
        # Append data values to list
        meanvals.append(np.mean(datarepeated))
        stdvals.append(np.std(datarepeated))
        minvals.append(np.amin(datarepeated))
        maxvals.append(np.amax(datarepeated))

        # Create array holding values 0 or 1 for all possible signals
        # If mean - min val is larger than 5 sigma, mark as signal
        # If max val - mean is larger than 5 sigma, mark as signal
        if (meanvals[i] - minvals[i]) > (5*stdvals[i]) or (maxvals[i] - meanvals[i]) > (5*stdvals[i]):
            # Signal
            sigvals.append(1)

        else:
            # No signal
            sigvals.append(0)

    return(meanvals, stdvals, minvals, maxvals, sigvals)


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
