import matplotlib.pyplot as plt
import numpy as np
import math as m


# Linear Fit
def linfit(x,y,z):
    """
    Method for finding the linear fit of a 2D plot for multiple non-signal events.
    Returns the numpy array of the y axis points for the best fit line

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

        print("Fitting Results: ", pfit[i], stats[i], sep='\n')

        # labelling variables from fitting for ease of use
        c[i], m[i] = pfit[i]
        # Not list based as they are not saved
        resid, rank, sing_val, rcond = stats[i]
        rms.append(np.sqrt(resid[0]/len(data)))

        print("Fit: y = " + str(m[i]) + "x + " + str(c[i]))
        print("RMS residual = " + str(rms[i]))



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
        datarepeated = branches['ADC'][i]
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
