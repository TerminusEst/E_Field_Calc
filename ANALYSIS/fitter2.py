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
    x2, y2 = rolling_average(X, Y, window/2.)
    x3, y3 = limit_lat(x2, y2, limit)
    x4, y4 = get_upper_limit(x3, y3)
    return x4, y4
    
def RMSD(m, o):
    return np.sqrt((1./len(m)) * np.sum((o - m)**2))
    
    
################################################################################
# read in all of the data first
folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/lat_data/"
dates = sorted(os.listdir(folder))
dates.remove("2003_03_30.txt")
dates.remove("2000_07_15.txt")

master_X, master_Y = [], []
for i in dates:
    X, Y = np.loadtxt(folder + i, unpack = True, usecols = (0, 1))
    master_X.append(X)
    master_Y.append(Y)
    
# maximum dst for the events
dst = array([-293, -354, -288, -288, -301, -387, -292, -383, -422, -374, -247])

################################################################################

iii = 7
X, Y = master_X[iii], master_Y[iii]
x2, y2 = get_slope_section(X, Y, 35, 10)

movx, movy = rolling_average(X, Y, 10)

#-------------------------------------------------------------------------------

def func1(x, a, b, c, d):
    # tanh
    return a * np.tanh(b*x - c) + d

def func2(x, a, b, c, d):
    #1/e
    return a * 1./(1+np.exp(-1*b*x+c)) + d

def func3(x, m, s, d, c):
    return d + (c/np.sqrt(2*np.pi*s))*np.exp((-(x-m)**2)/(2*s))
    
def func4(x, a, b, c, f):
    return a*(x-f)**2 + b*(x-f) + c


a1 = 0.8*np.mean(x2)
b1 = 1.0
c1 = np.median(x2)
d1 = np.mean(y2)

p01 = [a1, b1, c1, d1]
params1 = curve_fit(func1, x2, y2, p01, maxfev = 100000)
aa1, bb1, cc1, dd1 = params1[0]

#------------------------------------------------------
ind1 = X>40
ind2 = X < X[list(Y).index(max(Y))]
xqq = X[ind1*ind2]
yqq = Y[ind1*ind2]

a2 = 1.6 * np.mean(y2)
b2 = b1
c2 = c1
d2 = 0.25*d1

p02 = [a2, b2, c2, d2]
params2 = curve_fit(func2, x2, y2, p02, maxfev = 100000)
aa2, bb2, cc2, dd2 = params2[0]

#------------------------------------------------------

d = np.mean(y2)
m = 65
s = 1.
c = np.max(y2)/2.

p03 = [m, s, d, c]
#params3 = curve_fit(func3, x2, y2, p03, maxfev = 100000)
params3 = curve_fit(func3, X[X>40], Y[X>40], p03, maxfev = 100000)
m3, s3, d3, c3 = params3[0]

#------------------------------------------------------

def func4(x, a, b, c, d, f):

    return a*(x-f)**3 + b*(x-f)**2 + c*(x-f)+d
    
a = 0.
b = -1.
c = 0.
d = np.max(y2)
f = np.max(x2)

p04 = [a, b, c, d, f]
#params3 = curve_fit(func3, x2, y2, p03, maxfev = 100000)
params4 = curve_fit(func4, X[X>40], Y[X>40], p04, maxfev = 100000)
a4, b4, c4, d4, f4 = params4[0]

#------------------------------------------------------

x3 = np.linspace(40, 90, 1000)
Y1 = func1(x3, aa1, bb1, cc1, dd1)
Y2 = func2(x3, aa2, bb2, cc2, dd2)
Y3 = func3(x3, m3, s3, d3, c3)
Y4 = func4(x3, a4, b4, c4, d4, f4)
#print("TANH: ", RMSD(y2, Y1))
#print("1/e:  ", RMSD(y2, Y2))

grad1 = np.gradient(Y1)
max1 = x3[list(grad1).index(max(grad1))]
grad2 = np.gradient(Y2)
max2 = x3[list(grad2).index(max(grad2))]

minY1 = min(Y1)
for xxx, yyy in zip(x3, Y1):
    if yyy >= minY1 * e:
        break
"""
aaa = -3
bbb = 8.
ccc = np.max(Y) -600
fff = 68

clf()
plot(X, Y, 'o')
plot(x3, Y4)
plot(x3, Y2)
show()


"""
clf()
suptitle("Max TANH: %0.2f Max 1/e: %0.2f" % (max1, max2), fontsize = 24)
ax1 = subplot(211)
plot(X, Y, 'o', label = "original")
plot(x2, y2, 'o', label = "moving average")
plot(x3, Y1, 'r', label = "tanh")
#plot(x3, Y2, 'b', label = "1/e")
axvline(max1, color = 'r', linestyle = 'dashed')
#axvline(max2, color = 'b', linestyle = 'dashed')
axvline(xxx, color = "orange", linestyle = "dashed")

legend()
xlim([0, 90])
grid(True)

ax2 = subplot(212)
plot(x3, grad1, 'r', label = "tanh")
#plot(x3, grad2, 'b', label = "1/e")

axvline(max1, color = 'r', linestyle = 'dashed')
#axvline(max2, color = 'b', linestyle = 'dashed')
xlim([min(x2) - 5, max(x2) + 5])
grid(True)
xlim([0, 90])

show()

"""
x2, y2 = rolling_max(X, Y, 2.5)
xlog, ylog = rolling_max(X, np.log10(Y), 2.5)

clf()
ax1 = subplot(211)
plot(X, Y, 'o')
plot(x2, y2)

ax2 = subplot(212)
plot(X, np.log10(Y), 'o')
plot(xlog, ylog)

show()
"""






