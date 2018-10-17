#ffmpeg -r 30 -f image2 -i %04d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 10 -pix_fmt yuv420p test3.mp4
#/usr/bin/ffmpeg -r 30 -f image2 -i %04d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 10 -pix_fmt yuv420p test3.mp4




from spacepy import pycdf
import matplotlib.dates as mdates

from mpl_toolkits.mplot3d import Axes3D
#import healpy.visufunc as hpvf
from mpl_toolkits.basemap import Basemap, addcyclic
import os
import astropy
import astropy.coordinates as asc

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

    Vr = (np.cos(el) * np.cos(az) * x) + (np.cos(el) * np.sin(az) * y) + (np.sin(el) * z)
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

def read_ELF_file(filename, lat_cutoff):
    """ Read in iono files downloaded via CCMC interface """
    data = np.loadtxt(filename)
    lat, mlt, Jh = data[:,0], data[:,1], data[:,4]
    lon = (np.around(mlt*15) + 180.0) % 360.

    ind1 = lat >= lat_cutoff
    lon1, lat1, Jh1 = lon[ind1], lat[ind1], Jh[ind1]

    ind2 = lat <= -1 * lat_cutoff
    lon2, lat2, Jh2 = lon[ind2], lat[ind2], Jh[ind2]
    return lon1, lat1, Jh1, lon2, lat2, Jh2

def get_time_ELF_file(iii):
    """Get datetime from iono file downloaded via CCMC interface"""
    ymd = iii.split("_")[2]
    year = int(ymd[:4])
    month = int(ymd[4:6])
    day = int(ymd[6:])

    hmn = iii.split("_")[3]
    hour = int(hmn[:2])
    minute = int(hmn[2:4])

    timedate = datetime.datetime(year, month, day, hour, minute)
    return timedate


def read_CDF_file(filename):
    cdf = pycdf.CDF(filename)

    x = cdf['x'][0]
    y = -1* cdf['y'][0]
    z = cdf['z'][0]
    jx = cdf['jx'][0]
    jy = cdf['jy'][0]
    jz = cdf['jz'][0]
    jr = cdf['jr'][0]
    jh = np.sqrt(jx**2 + jy**2 + jz**2)

    lon, lat = xyz_to_lonlat(x, y, z)
    Vaz, Vel, Vr = cart2sphvec(jx, jy, jz, lon, lat)
    Vh = np.sqrt(Vaz**2 + Vel**2)

    return lon, lat, Vh

##########################################################################################
##########################################################################################
input_folders = ["/media/blake/Elements/Events/Sean_Blake_071618_1/Ionosphere/",
                "/media/blake/Elements/Events/Sean_Blake_071718_1/Ionosphere/"]

image_folders = ["/media/blake/Elements/Events/Sean_Blake_071618_1/Ionosphere_Images/",
                "/media/blake/Elements/Events/Sean_Blake_071718_1/Ionosphere_Images/"]

for iono_folder, image_folder in zip(input_folders, image_folders):
    files_iono = sorted(os.listdir(iono_folder))
    #-----------------------------------------------------------------------------------------

    # read in all cdfs, get maximum
    maxes_jh, timedate_xaxis = [], []
    for i in files_iono:
        #lon1, lat1, Jh1, lon2, lat2, Jh2 = read_ELF_file(iono_folder + i, 40.)
        #Jh1, Jh2 = Jh1/1e6, Jh2/1e6
        #timedate = get_time_ELF_file(i)    

        lon, lat, Jh = read_CDF_file(iono_folder + i)
        maxes_jh.append(max(Jh)/1e6)
        timedate, dayfrac, loncorrect = get_date(i)
        timedate_xaxis.append(timedate)
        print(i)

    #-----------------------------------------------------------------------------------------

    counter = 0.
    for f in files_iono:
        print(counter)
        counter += 1.
        #lon1, lat1, Jh1, lon2, lat2, Jh2 = read_ELF_file(iono_folder + f, 40.)
        #Jh1, Jh2 = Jh1/1e6, Jh2/1e6
        #timedate = get_time_ELF_file(f)

        lon, lat, Jh = read_CDF_file(iono_folder + f)
        Jh = Jh/1e6
        timedate, dayfrac, loncorrect = get_date(f)

        ind1 = lat >= 40.
        lon1, lat1, Jh1 = lon[ind1], lat[ind1], Jh[ind1]

        ind2 = lat <= -1 * 40.
        lon2, lat2, Jh2 = lon[ind2], lat[ind2], Jh[ind2]

        ##########################################################################################
        ##########################################################################################

        fig1 = figure(1)
        clf()

        suptitle(str(timedate), fontsize = 24)

        ax1 = subplot2grid((5, 4), (0, 0), colspan = 2, rowspan = 2)

        m = Basemap(projection='ortho',lat_0=90.,lon_0=180.,resolution='c')
        m.drawmapboundary(fill_color='black')
        m.drawmeridians(np.arange(0,360,45), color = 'w', zorder = 101)
        m.drawparallels(np.arange(-90,90,45), color = 'w', zorder = 101)

        title("Noon", fontsize = 20)
        ylabel("Dawn", fontsize = 20)

        xxx1, yyy1 = m(lon1, lat1)
        tripcolor(xxx1, yyy1, Jh1, cmap = "jet", zorder = 100.)#, vmax = maxx, vmin = minn)
        #scatter(xxx1, yyy1, c = jr1, zorder = 1000.)

        #m.drawcoastlines(linewidth=0.25)
        #m.drawcountries(linewidth=0.25)
        #m.fillcontinents(color='coral',lake_color='aqua')
        #colorbar()
        #-------------------------------------------------------------------------------------

        ax2 = subplot2grid((5, 4), (0, 2), colspan = 2, rowspan = 2)

        m = Basemap(projection='ortho',lat_0=-90.,lon_0=0,resolution='c')
        m.drawmapboundary(fill_color='black')
        m.drawmeridians(np.arange(0,360,45), color = 'w', zorder = 101.)#, labels = [1, 0, 0, 1])
        m.drawparallels(np.arange(-90,90,45), color = 'w', zorder = 101.)#, labels = [0, 1, 1, 0])

        title("Noon", fontsize = 20)
        ylabel("Dawn", fontsize = 20)
        ax2.yaxis.set_label_position("right")

        xxx2, yyy2 = m(lon2, lat2)
        tripcolor(xxx2, yyy2, Jh2, cmap = "jet", zorder = 100.)#, vmax = maxx, vmin = minn)
        #scatter(xxx1, yyy1, c = jr2, zorder = 1000.)

        #m.drawcoastlines(linewidth=0.25)
        #m.drawcountries(linewidth=0.25)
        #m.fillcontinents(color='coral',lake_color='aqua')
        #plt.colorbar(orientation="horizontal", vmax = max(maxes_jh))#,fraction=0.07,anchor=(1.0,0.0))


        #-------------------------------------------------------------------------------------
        ax3 = subplot2grid((5, 4), (2, 0), colspan = 2, rowspan = 2)

        m = Basemap(projection='ortho',lat_0=90.,lon_0=180.,resolution='c')
        m.drawmapboundary(fill_color='black')
        m.drawmeridians(np.arange(0,360,45), color = 'w', zorder = 101)
        m.drawparallels(np.arange(-90,90,45), color = 'w', zorder = 101)

        ylabel("Dawn", fontsize = 20)
        mycmap = plt.cm.jet

        xxx1, yyy1 = m(lon1, lat1)
        cs1 = tripcolor(xxx1, yyy1, Jh1, cmap = mycmap, zorder = 100.)#, vmax = max(maxes_jh), vmin = minn)
        cs1.cmap.set_under('k')
        cs1.set_clim(0.1, max(maxes_jh))

        #m.drawcoastlines(linewidth=0.25)
        #m.drawcountries(linewidth=0.25)
        #m.fillcontinents(color='coral',lake_color='aqua')
        #colorbar()
        #-------------------------------------------------------------------------------------

        ax4 = subplot2grid((5, 4), (2, 2), colspan = 2, rowspan = 2)

        m = Basemap(projection='ortho',lat_0=-90.,lon_0=0,resolution='c')
        m.drawmapboundary(fill_color='black')
        m.drawmeridians(np.arange(0,360,45), color = 'w', zorder = 101.)#, labels = [1, 0, 0, 1])
        m.drawparallels(np.arange(-90,90,45), color = 'w', zorder = 101.)#, labels = [0, 1, 1, 0])

        #ylabel("Dawn", fontsize = 20)
        ax4.yaxis.set_label_position("right")

        xxx2, yyy2 = m(lon2, lat2)
        cs2 = tripcolor(xxx2, yyy2, Jh2, cmap = mycmap, zorder = 100.)#, vmax = max(maxes_jh)/1e6, vmin = minn)
        cs2.cmap.set_under('k')
        cs2.set_clim(0.1, max(maxes_jh))
        cb = plt.colorbar(fraction=0.046, pad=0.04)

        #-------------------------------------------------------------------------------------
        ax5 = subplot2grid((5, 4), (4, 0), colspan = 4, rowspan = 1)

        plot(timedate_xaxis, array(maxes_jh))
        axvline(timedate, lw = 1, linestyle = "dashed", color = 'r')
        grid(True)

        gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,2)))
        xlabel("29 October 2003", fontsize = 24)
        ylabel(r"Max $J_{H}$ x 1e6", fontsize = 24)
        savefig(image_folder + "%04d" % (counter) + ".png")
        #show()
        #break



"""
import os
path = '/home/blake/Downloads/Sean_Blake_092818_4_gif_images_293267079509/'
files = sorted(os.listdir(path))
i = 1

for file in files:
    endnumber = "%04d.gif" % (i)
    os.rename(os.path.join(path, file), os.path.join(path, endnumber))
    i = i+1
"""





