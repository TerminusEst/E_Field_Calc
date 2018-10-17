from spacepy import pycdf
from mpl_toolkits.mplot3d import Axes3D
#import healpy.visufunc as hpvf
from mpl_toolkits.basemap import Basemap, addcyclic
import os
import astropy
import astropy.coordinates as asc
import scipy.interpolate
######################################################################################

def xyz_to_lonlat(x,y,z):
    """Convert cartesian points to latitude longitude"""

    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    #elev = np.arccos(z/np.sqrt(XsqPlusYsq))     # theta
    elev = np.arccos(z/r)
    az = np.arctan2(y,x)                           # phi
    return np.rad2deg(az), 90 - np.rad2deg(elev)   # want 

def cart2sphvec(x, y, z, az, el, degrees = True):
    """Convert cartesian vectors to spherical vectors"""
    if degrees == True:
        el = np.deg2rad(el)
        az = np.deg2rad(az)

    Vr = (np.cos(el) * np.cos(az) * jx) + (np.cos(el) * np.sin(az) * jy) + (np.sin(el) * jz)
    Vaz = (-1 * np.sin(az) * x) + (np.cos(az) * y)
    Vel = (-1 * np.sin(el) * np.cos(az) * x) + (-1 * y * np.sin(el) * np.sin(az)) + (z * np.cos(el))

    return (Vaz, Vel, Vr)

def get_date(aaa):
    year = int(aaa[13:17])
    month = int(aaa[17:19])
    day = int(aaa[19:21])
    hour = int(aaa[22:24])
    minute = int(aaa[24:26])

    timedate = datetime.datetime(year, month, day, hour, minute)
    dayfrac = ((minute/60.) + hour)/24.
    loncorrect = ((1 - dayfrac)*360) - 180.

    return timedate, dayfrac, loncorrect


##########################################################################################
##########################################################################################
folder = "/home/blake/Drive/NASA_Shite/CDF_TEST/cdf_data/"
files_iono = sorted(os.listdir(folder))

out_folder = "/home/blake/Drive/NASA_Shite/CDF_TEST/images/"
file_counter = 1

cdf_file = "/home/blake/Drive/NASA_Shite/CDF_TEST/cdf_data/null.swmf.i_e20120723-153300-000.cdf"
mag_file = "/home/blake/Drive/NASA_Shite/CDF_TEST/mag_data/mag_grid_e20120723-155330..txt"

counter += 1.
cdf = pycdf.CDF(cdf_file)
x = cdf['x'][0]
y = cdf['y'][0]
z = cdf['z'][0]
jx = cdf['jx'][0]
jy = cdf['jy'][0]
jz = cdf['jz'][0]
jr = cdf['jr'][0]
jh = np.sqrt(jx**2 + jy**2 + jz**2)

lon, lat = xyz_to_lonlat(x, y, z)
Vaz, Vel, Vr = cart2sphvec(jx, jy, jz, lon, lat)
Vh = np.sqrt(Vaz**2 + Vel**2)

# get top/bottom hemispheres
lat_cutoff = 40.
ind1 = lat >= lat_cutoff
lon1, lat1 = lon[ind1], lat[ind1]
jr1 = Vh[ind1]

lon_points1 = [180., 280.]
lat_points1 = [80., 80.]

lons1 = np.linspace(lon_points1[0], lon_points1[1], 100)
lats1 = np.linspace(lat_points1[0], lat_points1[1], 100)

x, y = m1(lon1, lat1)
xx, yy = m1(lons1, lats1)

interp_data1 = scipy.interpolate.griddata(list(zip(x, y)), jr1, list(zip(xx, yy)), method = 'cubic')
interp_data2 = scipy.interpolate.griddata(list(zip(x, y)), jr1, list(zip(xx, yy)), method = 'linear')


def interp_spherical_data(x1, y1, x2, y2, lon, lat, J, m1, great_circle = False):
    """Interpolate J data in a line between two points (x1, y1) and (x2, y2) in lat-lon.
    Can be a straight line or a great circle, depending on need"""

    background_x, background_y = m1(lon, lat)

    if great_circle == False:
        interp_lons = np.linspace(x1, x2, 100)
        interp_lats = np.linspace(y1, y2, 100)
        interp_x, interp_y = m1(interp_lons, interp_lats)

    else:
        great_circle = m1.drawgreatcircle(x1, y1, x2, y2, del_s = 50., alpha = 1.0, color = "w", zorder = 102)
        interp_x = great_circle[0].get_data()[0]   
        interp_y = great_circle[0].get_data()[1]

    output_data = scipy.interpolate.griddata(list(zip(background_x, background_y)), jr1, list(zip(interp_x, interp_y)), method = 'cubic')

    output_lon, output_lat = m1(interp_x, interp_y, inverse = True)
    return output_data, output_lon, output_lat

#line = m1.drawgreatcircle(lon_points1[0],lat_points1[0],lon_points1[1],lat_points1[1],del_s=1000,color='w', zorder = 1e7)
##########################################################################################
x1, y1 = lon_points1[0], lat_points1[0]
x2, y2 = lon_points1[1], lat_points1[1]

fig1 = figure(1)

clf()
ax1 = subplot(121)
suptitle(str(timedate), fontsize = 24)

m1 = Basemap(projection='ortho',lat_0=90.,lon_0=0 + 180,resolution='c')
m1.drawmapboundary(fill_color='black')
m1.drawmeridians(np.arange(0,360,30), color = 'w', zorder = 101)
m1.drawparallels(np.arange(-90,90,30), color = 'w', zorder = 101)

xxx1, yyy1 = m1(lon1, lat1)
tripcolor(xxx1, yyy1, jr1, cmap = "jet", zorder = 100.)#, vmax = maxx, vmin = minn)
colorbar()


xs1, ys1 = m1(lons1, lats1)
plot(xs1, ys1, 'r', zorder = 102., lw = 3)
x, y, z = interp_spherical_data(x1, y1, x2, y2, lon1, lat1, jr1, m1, great_circle = True)


ax2 = subplot(122)
a, b, c = interp_spherical_data(x1, y1, x2, y2, lon1, lat1, jr1, m1, great_circle = False)

plot(c, a)
plot(z, x)
show()








