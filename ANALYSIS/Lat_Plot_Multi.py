import os

################################################################################

def rolling_average(X, Y, window):
    values = np.arange(-80, 80, 1)
    means = []
    for i in values:
        value_left, value_right = i-window, i+window
        indices = np.logical_and(X >= value_left, X <= value_right)
        means.append(np.mean(Y[indices]))
    return np.array(values), np.array(means)

def slope_changes(X, Y):    
    slope = np.gradient(Y)
    posneg = (slope > 0) * 1
    posneg2 = np.roll(posneg, 1)

    diff_posneg = abs(posneg - posneg2)
    diff_posneg[0] = diff_posneg[1]
    diff_posneg[-1] = diff_posneg[-2]
    
    indices = []
    for index, value in enumerate(X):
        if X[index] > 0:
            if diff_posneg[index] != 0:
                indices.append(index)
    
    return diff_posneg, indices
    
################################################################################
# read in all of the data first
folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/lat_data/"
dates = sorted(os.listdir(folder))
dates.remove("2003_03_30.txt")


################################################################################
# big old multiplot

master_X, master_Y = [], []
for i in dates:
    X, Y = np.loadtxt(folder + i, unpack = True, usecols = (0, 1))
    master_X.append(X)
    master_Y.append(Y)
    

fig = figure(1)
clf()    

for iii, vvv in enumerate(dates):

    X, Y = master_X[iii], master_Y[iii]
    XX, YY = rolling_average(X, Y, 7.5)


    ax = subplot(len(dates), 1, iii+1)

    plot(X, Y/1000., 'o', label = vvv[:-4])
    plot(XX, YY/1000., 'r')
    
    if iii != len(dates) - 1:
        plt.setp( ax.get_xticklabels(), visible=False)
    legend(loc = "center")

    xlim([-90, 90])
    grid(True)

suptitle("Moving Window (15 Degrees)", fontsize = 20)
fig.text(0.03, 0.5, 'Eh (V/km)', va='center', rotation='vertical', fontsize = 22)
xlabel("Magnetic Lat", fontsize = 22)
show()    

################################################################################
    
dst = [-293, -354, -288, -288, -301, -387, -292, -383, -422, -374, -247]

points = []

for iii in range(len(dates)):
    
    X, Y = master_X[iii], master_Y[iii]     # data
    XX, YY = rolling_average(X, Y, 7.5) # rolling mean

    # remove nans and interpolate
    indices = ~np.isnan(YY)

    XX2 = np.arange(-80, 80, 0.5)
    YY2 = np.interp(XX2, XX[indices], YY[indices])

    slope = np.gradient(YY2)

    maxslope = np.max(slope[XX2 > 0])
    maxslope_index = list(slope).index(maxslope)
    maxslope_lat = XX2[maxslope_index]

    points.append(maxslope_lat)
    
clf()

scatter(dst, points)

show()
    
"""

clf()

ax1 = subplot(211)
plot(X, Y, 'o')
plot(XX2, YY2)

ax2 = subplot(212)
plot(XX2, slope)
axhline(0, color = "black", linestyle = "dashed")

axvline(maxslope_lat, color = "red", linestyle = "dashed")

show()
"""
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    








