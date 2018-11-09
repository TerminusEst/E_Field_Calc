""" Read in the Lat vs. Maximum Eh and try to fit lines """

import os
import numpy as np
from scipy.optimize import curve_fit

################################################################################
# Functions

def rolling_average(X, Y, window):
    values = np.arange(-80, 80, 1)
    means = []
    for i in values:
        value_left, value_right = i-window, i+window
        indices = np.logical_and(X >= value_left, X <= value_right)
        means.append(np.mean(Y[indices]))
    return np.array(values), np.array(means)
    
################################################################################
# read in all of the data first
folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/lat_data/"
dates = sorted(os.listdir(folder))
dates.remove("2003_03_30.txt")

master_X, master_Y = [], []
for i in dates:
    X, Y = np.loadtxt(folder + i, unpack = True, usecols = (0, 1))
    master_X.append(X)
    master_Y.append(Y)
    
# maximum dst for the events
dst = array([-293, -354, -288, -288, -301, -387, -292, -383, -422, -374, -247])

################################################################################
# Select initial latitude range
clf()
for iii, bleh in enumerate(master_X[:8]):
    try:
        X, Y = master_X[iii], master_Y[iii]

        # select latitudes > 25
        for index, value in enumerate(X):
            if value >= 30:
                break

        X, Y = X[index:], Y[index:]
        X2, Y2 = rolling_average(X, Y, 2.5)

        index2 = list(Y2).index(np.nanmax(Y2))
        lat2 = X2[index2+4]

        # select latitudes > 25
        X3, Y3 = [], []
        for index, value in enumerate(X):
            if 45 <= value <= lat2:
                X3.append(X[index])
                Y3.append(Y[index])
        #X3 = X3 - np.mean(X3)

        ################################################################################
        # try to fit the function


        def fit_func(x, a, b, c, d): 
            return a * np.tanh(b*x - c) + d
          
        a = 0.8*np.mean(Y3)
        b = 1.0
        c = np.median(X3)
        d = np.mean(Y3)

        p0 = [a, b, c, d]
        params = curve_fit(fit_func, X3, Y3, p0, sigma = np.sqrt(Y3), maxfev = 100000)

        aa, bb, cc, dd = params[0]

        x = np.linspace(min(X3), max(X3), 100)
        y = fit_func(x, a, b, c, d)
        yy = fit_func(x, aa, bb, cc, dd)


        ax = subplot(8, 1, iii + 1)
        scatter(X3, Y3)
        #plot(x, y)
        plot(x, yy, 'r')
        if iii != len(master_X) - 1:
            plt.setp( ax.get_xticklabels(), visible=False)
            
        xlim([45, 80])
        grid(True)
    except:
        print(iii, dates[iii])
show()



    
"""
p0 = [0.4*np.max(Y3), 0.4*np.max(Y3), 0.3]

#params = curve_fit(fit_func, X3, Y3, p0, maxfev = 10000)
#aa, bb, cc = params[0]

#fitX = np.linspace(min(X3), max(X3), 1000)
#fitY = fit_func(fitX, aa, bb, cc)


x = np.linspace(min(X3), max(X3), 1000)

clf()

scatter(X3, Y3)
#plot(fitX, fitY)
#plot(x, fit_func(x, a, b, c))
                          
#plot(x, 0.5*np.max(Y3) + 0.5*np.max(Y3) * np.tanh(1.1*x - mean(X3)))

show()
"""




































