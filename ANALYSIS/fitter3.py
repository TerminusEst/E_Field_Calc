""" Read in the Lat vs. Maximum Eh and try to fit lines """

import seaborn as sns
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

################################################################################
# Functions

def rolling_average(X, Y, window):
    values = np.arange(min(X), max(X), 1)
    means = []
    for i in values:
        value_left, value_right = i-window, i+window
        indices = np.logical_and(X >= value_left, X <= value_right)
        means.append(np.mean(Y[indices]))
    return np.array(values), np.array(means)

def rolling_max(X, Y, window):
    values = np.arange(-80, 80, 1)
    means = []
    for i in values:
        value_left, value_right = i-window, i+window
        indices = np.logical_and(X >= value_left, X <= value_right)
        try:
            means.append(np.max(Y[indices]))
        except:
            means.append(np.nan)
    return np.array(values), np.array(means)

def rolling_average2(X, Y, window):
    means = []
    for i in X:
        value_left, value_right = i - window/2., i + window/2.
        indices = np.logical_and(X >= value_left, X <= value_right)
        means.append(np.mean(Y[indices]))
    return X, np.array(means)

def get_upper_limit(x, y):
    maxx = np.max(y)
    maxx_x = x[list(y).index(maxx)]
    
    ind = x <= maxx_x    
    return x[ind], y[ind]

def limit_lat(x2, y2, limit):
    ind = x2 >= limit
    return x2[ind], y2[ind]

def get_slope_section(X, Y, limit, window):
    x2, y2 = rolling_average2(X, Y, window/2.)
    x3, y3 = limit_lat(x2, y2, limit)
    x4, y4 = get_upper_limit(x3, y3)
    return x4, y4
    
def RMSD(m, o):
    return np.sqrt((1./len(m)) * np.sum((o - m)**2))
    
def limit_lat2(x, y, llimit, rlimit):
    ind1 = x >= llimit
    xx, yy = x[ind1], y[ind1]
    
    #rlimit = xx[yy.argmax()] + 5
    ind2 = xx <= rlimit
    
    return xx[ind2], yy[ind2]
    
################################################################################
# read in all of the data first
folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/lat_data/"
dates = sorted(os.listdir(folder))
#dates.remove("2003_03_30.txt")
#dates.remove("2000_07_15.txt")
#dates = dates[2:]
#dates = ["2003_03_30.txt"]

DST = []
master_X, master_Y = [], []
for i in dates:
    X, Y = np.loadtxt(folder + i, unpack = True, usecols = (0, 1), skiprows = 1)
    master_X.append(X)
    master_Y.append(Y)
    
    # maximum dst for the events
    f = open(folder + i)
    data = float(f.readlines()[0])
    f.close()
    DST.append(data)
################################################################################

fig = figure(1)
clf()

iii = 2
maxes, maxgrads = [], []

for iii in range(len(dates)):
    X, Y = master_X[iii], master_Y[iii]
    Y = np.log10(Y)

    max_x = X[Y.argmax()]
    #x2, y2 = get_slope_section(X, Y, 35, 10)
    x2, y2 = limit_lat2(X, Y, 20, 90)

    #a, b = rolling_average(x2, 10**y2, 5)
    #x2, y2 = a, np.log10(b)

    #------------------------------------------------------

    def func1(x, a, b, c, d):
        # tanh
        return a * np.tanh(b*x - c) + d

    a1 = 0.8*np.mean(x2)
    b1 = 1.0
    c1 = np.median(x2)
    d1 = np.mean(y2)

    p01 = [a1, b1, c1, d1]
    params1 = curve_fit(func1, x2, y2, p01, sigma = y2, absolute_sigma = False, maxfev = 100000)
    aa1, bb1, cc1, dd1 = params1[0]

    x3 = np.linspace(min(x2), max(x2), 10000)
    y3 = func1(x3, aa1, bb1, cc1, dd1)

    grad3 = np.gradient(y3)
    maxgrad = x3[grad3.argmax()]

    #------------------------------------------------------

    maxyvalue = (10**max(Y))/1000.
    #maxyvalue = np.mean(10**array(sorted(Y, reverse = True)[:3]))/1000.

    ax = subplot(len(dates), 2, (iii*2) + 1)

    semilogy(X, 10**Y, 'o')
    semilogy(x2, 10**y2, 'o', color = "orange")
    semilogy(x3, 10**y3, color = "green")
    xlim([0, 90])
    grid(True)
    axvline(maxgrad, linestyle = "dashed", color = "red")
    if iii != len(dates)-1:
        ax.set_xticklabels([])
    if iii == 0:
        title("Calculated Eh", fontsize = 24)
    text(0.2, 0.8, dates[iii][:-4], horizontalalignment='center', verticalalignment='center',
     transform = ax.transAxes, bbox=dict(facecolor='white', alpha=1))
    
    
    
    ax1 = subplot(len(dates), 2, (iii*2) + 2)
    plot(x3, np.gradient(y3), color = "green")
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")

    xlim([0, 90])
    grid(True)
    axvline(maxgrad, linestyle = "dashed", color = "red")
    
    if iii != len(dates)-1:
        ax1.set_xticklabels([])
    if iii == 0:
        title("Slope of fit Sigmoid", fontsize = 24)
    
    maxgrads.append(maxgrad)
    text(0.2, 0.8, "DST = " + str(DST[iii]) + " nT", horizontalalignment='center', verticalalignment='center',
     transform = ax1.transAxes, bbox=dict(facecolor='white', alpha=1))
    print(maxyvalue, dates[iii])
    
fig.text(0.5, 0.02, 'Magnetic Latitude', ha='center', fontsize = 24)

show()
    
slope, intercept, r_value, p_value, std_err = stats.linregress(DST, maxgrads)

fig2 = figure(2)
clf()
ax1 = subplot(111)
sns.regplot(DST, maxgrads)

for x, y, name in zip(DST, maxgrads, dates):
    text(x + 0.5, y + 0.25, name[:-4], horizontalalignment = "center", bbox=dict(facecolor='white', alpha=0.6))

ylabel("Mag Latitude", fontsize = 24)
xlabel("Minimum Dst (nT)", fontsize = 24)
ax1.legend()
grid(True)

text(0.1, 0.9, "Rsquared = %0.3f" % (r_value), fontsize = 22, horizontalalignment='center', verticalalignment='center',
     transform = ax1.transAxes, bbox=dict(facecolor='white', alpha=1))
     
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(DST, maxgrads)

    
    
    
    
    
    
