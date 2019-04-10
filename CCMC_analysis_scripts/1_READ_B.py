import os
from draco import EField as ef
import aacgmv2
import matplotlib.dates as mdates

#/usr/bin/ffmpeg -r 30 -f image2 -i %03d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 10 -pix_fmt yuv420p 2003_11_20_Scan.mp4

##########################################################################################

def time2float(x):
	"""converts datetime to float, so that interpolation/smoothing can be performed"""
	if (type(x) == numpy.ndarray) or (type(x) == list):
		emptyarray = []
		for i in x:
			z = (i - datetime.datetime(1970, 1, 1, 0)).total_seconds()
			emptyarray.append(z)
		emptyarray = numpy.array([emptyarray])
		return emptyarray[0]
	else:
		return (x - datetime.datetime(1970, 1, 1, 0)).total_seconds()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def float2time(x):
	"""converts array back to datetime so that it can be plotted with time on the axis"""
	if (type(x) == numpy.ndarray) or (type(x) == list):
		emptyarray = []
		for i in x:
			z = datetime.datetime.utcfromtimestamp(i)
			emptyarray.append(z)
		emptyarray = numpy.array([emptyarray])
		return emptyarray[0]
	else:
		return datetime.datetime.utcfromtimestamp(x)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def get_date(aaa, filetype):
    """get datetime from filename""" 
    if filetype == "iono":
        year = int(aaa[13:17])
        month = int(aaa[17:19])
        day = int(aaa[19:21])
        hour = int(aaa[22:24])
        minute = int(aaa[24:26])

    elif filetype == "mag":
        year = int(aaa[10:14])
        month = int(aaa[14:16])
        day = int(aaa[16:18])
        hour = int(aaa[19:21])
        minute = int(aaa[21:23])

    timedate = datetime.datetime(year, month, day, hour, minute)
    dayfrac = ((minute/60.) + hour)/24.
    loncorrect = ((1 - dayfrac)*360) - 180.

    #print(timedate, dayfrac, loncorrect)
    return timedate, dayfrac, loncorrect

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def calculate_magnetic_latitude(glon, glat, dtime):
    
    for height in arange(0, 100, 0.5):    
        calclat, calclon, z = aacgmv2.get_aacgm_coord(glat, glon, height, start)
        if np.isnan(calclat) != True:
            return calclat, calclon, height

    return nan, nan, nan
    
def get_td_1(a):
    b = a.split("_")
    year = int(b[2][:4])
    month = int(b[2][4:6])
    day = int(b[2][6:])
        
    hour = int(b[-1][:2])
    minute = int(b[-1][2:4])
    second = int(b[-1][4:6])
    
    td = datetime.datetime(year, month, day, hour, minute, second)
    return td

def get_td_2(a):
    b = a[10:]
    year = int(b[:4])
    month = int(b[4:6])
    day = int(b[6:8])
        
    hour = int(b[9:11])
    minute = int(b[11:13])
    second = int(b[13:15])
    td = datetime.datetime(year, month, day, hour, minute, second)
    return td

##########################################################################################
if 1:
    input_filename = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_1/DST_CCMC_output.txt"
    output_filename = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_1/Numpy_Data/DST"
    eff.read_write_dst(input_filename, output_filename)
##########################################################################################

folder_mag = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_1/Mag_Data/"

files_mag = sorted(os.listdir(folder_mag))
lenfiles = len(files_mag) * 1.0

geolon, geolat = np.loadtxt(folder_mag + files_mag[0], usecols = (0, 1), unpack = True, skiprows = 10)

#-------------------------------------------------------------------------------
# set up the master objects to be populated
# read in the files average the files for each minute
print("Reading in files")
timedate1 = []
masterBx = np.ones((int(len(files_mag)), len(geolon)))
masterBy, masterBz = np.copy(masterBx), np.copy(masterBx)

masterMx, masterMy, masterMz = np.copy(masterBx), np.copy(masterBx), np.copy(masterBx)
masterFx, masterFy, masterFz = np.copy(masterBx), np.copy(masterBx), np.copy(masterBx)
masterHx, masterHy, masterHz = np.copy(masterBx), np.copy(masterBx), np.copy(masterBx)
masterPx, masterPy, masterPz = np.copy(masterBx), np.copy(masterBx), np.copy(masterBx)

for i, v in enumerate(files_mag):
    dBn, dBe, dBd = np.loadtxt(folder_mag + v, skiprows = 4, usecols = (5, 6, 7), unpack = True)
    masterBx[i], masterBy[i], masterBz[i] = dBn, dBe, dBd
    timedate1.append(get_td_1(v))
    if i%10 == 0:
        print(i, v)

numpysavefolder = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_1/Numpy_Data/"
np.save(numpysavefolder + "BX", masterBx)
np.save(numpysavefolder + "BY", masterBy)

lonlat = np.concatenate((np.array([geolon]).T, np.array([geolat]).T), axis = 1)
np.save(numpysavefolder + "LONLAT", lonlat)
np.save(numpysavefolder + "TD", timedate1)

##########################################################################################























































