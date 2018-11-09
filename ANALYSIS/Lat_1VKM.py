
import os
import numpy as np
from scipy.optimize import curve_fit

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
    
def vkm_lat(x, y):
    for index, value in enumerate(y):
        if value >= 1000.:
            return x[index]
    
################################################################################
# read in all of the data first
folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/lat_data/"
dates = sorted(os.listdir(folder))
dates.remove("2003_03_30.txt")
dates.remove("2000_07_15.txt")

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

vkm_lats, maxe = [], []
for iii in range(len(dates)):
    X, Y = master_X[iii], master_Y[iii]
    
    x2, y2 = rolling_average(X, Y, 5)
    x3, y3 = limit_lat(x2, y2, 25)
    x4 = vkm_lat(x3, y3)
    vkm_lats.append(x4)
    maxe.append(max(y3))
    

def fit_line(x, y):
    p = np.polyfit(x, y, 1)
    
    return (p[0] * x) + p[1]

clf()

ax1 = subplot(311)

x, y = np.array(maxe)/1000., vkm_lats
#plot(x, fit_line(x, y))
sns.regplot(x, y)

xlabel("Max Calculated E-Field (V/km)", fontsize = 18)
ylabel("Latitude where 1VKM\n first seen", fontsize = 20)
grid(True)

ax2 = subplot(312)
#scatter(DST, vkm_lats)
sns.regplot(DST, vkm_lats)
xlabel("DST (nT)", fontsize = 18)
ylabel("Latitude where 1VKM\n first seen", fontsize = 20)
grid(True)

ax3 = subplot(313)
sns.regplot(np.array(maxe)/1000., DST)
xlabel("Max Calculated E-Field (V/km)", fontsize = 18)
ylabel("Min DST (nT)", fontsize = 20)
grid(True)
show()





























