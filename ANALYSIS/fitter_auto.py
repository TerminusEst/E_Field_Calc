import seaborn as sns
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
from gcvspline import GCVSmoothedNSpline
from draco import lat_stuff as ls

from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline
from gcvspline import GCVSmoothedNSpline, MSESmoothedNSpline, DOFSmoothedNSpline, SmoothedNSpline
################################################################################

    
################################################################################
# read in all of the data first
folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/lat_data/"
dates = sorted(os.listdir(folder))


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
    
DST, dates = (list(x) for x in zip(*sorted(zip(DST, dates), key=lambda pair: pair[0])))

xyz = 13
DST, dates = DST[:xyz], dates[:xyz]


################################################################################
fig1 = figure(1)
clf()
iii = 0

maxes, maxgrads = [], []

for iii in range(len(dates)):
    x, y = master_X[iii], master_Y[iii]
    x, y = ls.limit_lat2(x, y, 25, 90)
    ylog = np.log10(y)


    ################################################################################
    max_x = x[y.argmax()] + 50
    ind = x<max_x
    x = x[ind]
    y = ylog[ind]

    xs = np.linspace(x[0], x[-1], 1000)
    w = np.ones_like(x)

    GCV_auto = GCVSmoothedNSpline(x, y, w=w)
    p = 400#GCV_auto.p
    GCV_manual = SmoothedNSpline(x, y, w=w, p=p)
    ys = GCV_manual(xs)
    grady = np.gradient(ys)
    maxgrad = xs[grady.argmax()]

    ax = subplot(len(dates), 2, (iii*2) + 1)

    #plot(x, 10**y, 'o')
    #plot(xs, 10**ys, color = "red")
    semilogy(x, 10**y, 'o')
    semilogy(xs, 10**ys, color = "red")
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
    plot(xs, grady, color = "green")
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
    #print(maxyvalue, dates[iii])

################################################################################
#DST PLOT
slope, intercept, r_value, p_value, std_err = stats.linregress(DST, maxgrads)

show()
    
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
     
show()

