import matplotlib.pyplot as plt
import numpy as np
import math as m


# Linear Fit
def linfit(x,y):
    """
    Method for finding the linear fit of a 2D plot.
    Returns the numpy array of the y axis points for the best fit line

    :param time: x values
    :param data: y values
    :return: y values for straight line fit, fitting parameters(offset value, gradient value), stats(residuals, rank, singular values of matrix, rcond), root mean square residual value
    """
    # Find equation of straight line through data to determine baseline & baseline modulation
    time = x
    data = y
    # https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyfit.html
    Polynomial = np.polynomial.Polynomial

    # Fit to range cmin, cmax
    cmin, cmax = min(time), max(time)
    pfit, stats = Polynomial.fit(time, data, 1, full=True, window=(cmin, cmax), domain=(cmin, cmax))

    print("Fitting Results: ", pfit, stats, sep='\n')

    # labelling variables from fitting for ease of use
    c, m = pfit
    resid, rank, sing_val, rcond = stats
    rms = np.sqrt(resid[0]/len(data))

    print("Fit: y = " + str(m) + "x + " + str(c))
    print("RMS residual = " + str(rms))
    # Create line to plot
    yline = []
    print("Number of values: " + str(len(time)))
    for i in range(len(time)):
        yline.append(m*time[i]+c)

    return(yline, pfit, stats, rms)
