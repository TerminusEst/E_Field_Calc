from spacepy import pycdf
from mpl_toolkits.mplot3d import Axes3D
#import healpy.visufunc as hpvf
from mpl_toolkits.basemap import Basemap, addcyclic
import os

######################################################################################
"""
def xyz_to_lonlat(x, y, z):

    r = np.sqrt(x**2 + y**2 + y**2)

    # latitude:
    theta_temp = np.arccos(np.abs(z)/r)
    for i, v in enumerate(z):
        if v < 0:
            theta_temp[i]*=-1

    # longitude
    for xxx, yyy in zip(x, y):
        if (xxx >= 0):
            if (yyy >= 0):
                lon.append(np.rad2deg(np.arctan(yyy/xxx))) 
                lat.append(np.rad2deg(np.arccos(zzz/R)) - 90.) 
            else:
                lon.append(np.rad2deg(np.arctan(xxx/yyy)) + 270.)
                lat.append(np.rad2deg(np.arccos(zzz/R)) - 90.)
        else:
            if (yyy >= 0):
                lon.append(180. - np.rad2deg(np.arctan(yyy/xxx)))
                lat.append(np.rad2deg(np.arccos(zzz/R)) - 90.)
            else:
                lon.append(np.rad2deg(np.arctan(yyy/xxx)) + 180.)
                lat.append(np.rad2deg(np.arccos(zzz/R)) - 90.)

    return lon, lat, r

"""
def xyz_to_lonlat(x,y,z):

    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    elev = np.arctan2(z,np.sqrt(XsqPlusYsq))     # theta
    az = np.arctan2(y,x)                           # phi
    return np.rad2deg(az), np.rad2deg(elev)


######################################################################################
filename = "/home/blake/Drive/NASA_Shite/CDF_TEST/3d__var_1_e20120723-120000-000.out.cdf"

cdf = pycdf.CDF(filename)

x1 = cdf['x'][0]
y1 = cdf['y'][0]
z1 = cdf['z'][0]
bx1 = cdf['bx'][0]

R = np.sqrt(x1**2 + y1**2 + z1**2)

ind = np.logical_and(2.509 <= np.abs(R), np.abs(R) <= 2.7)
ind2 = z1 > 0
ind3 = ind * ind2

lon1, lat1 = xyz_to_lonlat(x1[ind3], y1[ind3], z1[ind3])
bx1 = bx1[ind3]

######################################################################################

filename2 = "/home/blake/Drive/NASA_Shite/CDF_TEST/3d_ASCII.txt"
data = np.loadtxt(filename2)

x2, y2, z2 = data[:,0], data[:,1], data[:,2]
bx2, by2, bz2 = data[:,3], data[:,4], data[:,5]

lon2, lat2 = xyz_to_lonlat(x2, y2, z2)

thingy = lat2 > 0 
lon2 = lon2[thingy]
lat2 = lat2[thingy]
bx2 = bx2[thingy]
######################################################################################
fig = figure(1)#, figsize = (20, 10))

minn = min(min(bx1), min(bx2))
maxx = max(max(bx1), max(bx2))

clf()

ax1 = subplot(221)
title("CDF File data", fontsize = 24)

m = Basemap(projection='ortho',lat_0=90.,lon_0=0,resolution='c')
m.drawmapboundary(fill_color='black')
m.drawmeridians(np.arange(0,360,30), color = 'w')
m.drawparallels(np.arange(-90,90,30), color = 'w')
#m.drawcoastlines(linewidth=0.25)
#m.fillcontinents(color='coral',lake_color='aqua')

xxx1, yyy1 = m(lon1, lat1)
#tripcolor(xxx1, yyy1, bx1, cmap = "jet", vmin = minn, vmax = maxx)
scatter(xxx1, yyy1, c = bx1, cmap = "jet", vmin = minn, vmax = maxx)
colorbar()

#-----------------------------------------------------------------------
ax2 = subplot(222)
title("CCMC Interface (Default)", fontsize = 24)

m = Basemap(projection='ortho',lat_0=90.,lon_0=0,resolution='c')
m.drawmapboundary(fill_color='black')
m.drawmeridians(np.arange(0,360,30), color = 'w')
m.drawparallels(np.arange(-90,90,30), color = 'w')

xxx2, yyy2 = m(lon2, lat2)
#tripcolor(xxx2, yyy2, bx2, cmap = "jet", vmin = minn, vmax = maxx)
scatter(xxx2, yyy2, c = bx2, cmap = "jet", vmin = minn, vmax = maxx)
colorbar()

#-----------------------------------------------------------------------
ax3 = subplot(223)

m = Basemap(projection='ortho',lat_0=90.,lon_0=0,resolution='c')
m.drawmapboundary(fill_color='black')
m.drawmeridians(np.arange(0,360,30), color = 'w')
m.drawparallels(np.arange(-90,90,30), color = 'w')

xxx1, yyy1 = m(lon1, lat1)
tripcolor(xxx1, yyy1, bx1, cmap = "jet", vmin = minn, vmax = maxx)
colorbar()

#-----------------------------------------------------------------------
ax4 = subplot(224)

m = Basemap(projection='ortho',lat_0=90.,lon_0=0,resolution='c')
m.drawmapboundary(fill_color='black')
m.drawmeridians(np.arange(0,360,30), color = 'w')
m.drawparallels(np.arange(-90,90,30), color = 'w')

xxx2, yyy2 = m(lon2, lat2)
tripcolor(xxx2, yyy2, bx2, cmap = "jet", vmin = minn, vmax = maxx)
colorbar()


show()







f = open("output.txt", 'w')
mystr = ""

for i in x1[ind3]:
    mystr += "%.4f" % (i) + ","

mystr += "\n"
for i in y1[ind3]:
    mystr += "%.4f" % (i) + ","

mystr += "\n"
for i in z1[ind3]:
    mystr += "%.4f" % (i) + ","


f.write(mystr)
f.close()



filename = "coord_output.txt"
np.savetxt(filename, list(zip(x1[ind3], y1[ind3], z1[ind3])))
















































