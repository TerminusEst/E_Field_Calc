import os
from draco import EField as ef
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

##########################################################################################
"""
resistivities, thicknesses = ef.Model_Profiles("Q")
folder_mag = "/media/blake/Elements/Events/Sean_Blake_092818_1/Surface_Mag/"
files_mag = sorted(os.listdir(folder_mag))

maglon, maglat = np.loadtxt(folder_mag + files_mag[0], skiprows = 4, unpack = True, usecols = (0, 1))
masterBx, masterBy = np.ones((int(len(files_mag)), len(maglon))), np.ones((int(len(files_mag)), len(maglon)))
masterEx, masterEy = np.ones((int(len(files_mag)), len(maglon))), np.ones((int(len(files_mag)), len(maglon)))

start = datetime.datetime(2012, 7, 23, 12, 0)
end = datetime.datetime(2012, 7, 24, 4, 47)

for i, v in enumerate(files_mag):
    print(v)
    data = np.loadtxt(folder_mag + v, skiprows=10)
    maglon, maglat, Bnn, Bee = data[:,0], data[:,1], data[:,3], data[:,4]

    masterBx[i] = Bnn
    masterBy[i] = Bee
"""

"""
for i in range(int(len(files_mag))):
    ct = start + datetime.timedelta(minutes = i)    # current time
    datestr = str(ct.year) + "%02d%02d-%02d%02d" % (ct.month, ct.day, ct.hour, ct.minute)
    print(datestr)

    Bn, Be, magcount = np.zeros(len(maglon)), np.zeros(len(maglon)), 0
    for f2 in files_mag:
        if datestr in f2:
            magcount += 1
            data = np.loadtxt(folder_mag + f2, skiprows = 10)
            maglon, maglat, Bnn, Bee = data[:,0], data[:,1], data[:,3], data[:,4]

            Bn += Bnn
            Be += Bee

    Bx = Bn/magcount
    By = Be/magcount

    masterBx[i] = Bx
    masterBy[i] = By

    break
"""
"""
##########################################################################################

for i in range(len(masterBy[0])):
    x, y = masterBx[:,i], masterBy[:,i]
    
    if i%100 == 0:
        print(i, "EFIELD")
    Ex, Ey = ef.E_Field_1D(x, y, resistivities, thicknesses, timestep = 60.)
    masterEx[:,i] = Ex
    masterEy[:,i] = Ey

masterBh = np.sqrt(masterBx**2 + masterBy**2)
masterdBh = np.gradient(masterBh, axis = 0)
masterEh = np.sqrt(masterEx**2 + masterEy**2)/1000.

##########################################################################################

lat_bins = sorted(list(set(maglat)))
lat_maxh = np.zeros(len(lat_bins))
lat_maxb = np.zeros(len(lat_bins))
lat_maxdb = np.zeros(len(lat_bins))

for i, v in enumerate(maglat):
    index = lat_bins.index(v)

    maxEh = max(masterEh.T[i][150:])
    maxBh = max(masterBh.T[i][150:])
    maxdBh = max(masterdBh.T[i][150:])

    if maxEh > lat_maxh[index]:
        lat_maxh[index] = maxEh
    if maxBh > lat_maxb[index]:
        lat_maxb[index] = maxBh
    if maxdBh > lat_maxdb[index]:
        lat_maxdb[index] = maxdBh

clf()

ax1 = subplot(311)
title("2012 Event, Inputs * 0.5", fontsize = 24)

plot(lat_bins, lat_maxh, label = "H")
legend()
ylabel("Eh (V/km)", fontsize = 20)
grid(True)

ax2 = subplot(312)
plot(lat_bins, lat_maxb, label = "Max Bh")
legend()
ylabel("Bh (nT)", fontsize = 20)
grid(True)

ax3 = subplot(313)
plot(lat_bins, lat_maxdb, label = "Max dBh")
legend()
ylabel("dBh (nT/min)", fontsize = 20)
xlabel("Latitude", fontsize = 20)

grid(True)
show()

##########################################################################################


"""

folder_mag2 = "/media/blake/Elements/Events/Sean_Blake_081318_1/Surface_Mag/"
files_mag2 = sorted(os.listdir(folder_mag2))

maglon2, maglat2 = np.loadtxt(folder_mag2 + files_mag2[0], skiprows = 4, unpack = True, usecols = (0, 1))
masterBx2, masterBy2 = np.ones((int(len(files_mag2)/2), len(maglon2))), np.ones((int(len(files_mag2)/2), len(maglon2)))
masterEx2, masterEy2 = np.ones((int(len(files_mag2)/2), len(maglon2))), np.ones((int(len(files_mag2)/2), len(maglon2)))

start = datetime.datetime(2012, 7, 23, 12, 0)
end = datetime.datetime(2012, 7, 24, 4, 47)

for i in range(int(len(files_mag2))):
    ct = start + datetime.timedelta(minutes = i)    # current time
    datestr = str(ct.year) + "%02d%02d-%02d%02d" % (ct.month, ct.day, ct.hour, ct.minute)
    print(datestr)

    Bn, Be, magcount = np.zeros(len(maglon2)), np.zeros(len(maglon2)), 0
    for f2 in files_mag2:
        if datestr in f2:
            magcount += 1
            data = np.loadtxt(folder_mag2 + f2, skiprows = 4)
            maglon, maglat, Bnn, Bee = data[:,0], data[:,1], data[:,2], data[:,3]

            Bn += Bnn
            Be += Bee

    Bx = Bn/magcount
    By = Be/magcount

    masterBx2[i] = Bx
    masterBy2[i] = By


masterBh2 = np.sqrt(masterBx2**2 + masterBy2**2)
masterdBh2 = np.gradient(masterBh2, axis = 0)
#masterEh = np.sqrt(masterEx**2 + masterEy**2)/1000.



for i in range(len(masterBy2[0])):
    x, y = masterBx2[:,i], masterBy2[:,i]
    
    if i%100 == 0:
        print(i, "EFIELD")
    Ex2, Ey2 = ef.E_Field_1D(x, y, resistivities, thicknesses, timestep = 60.)
    masterEx2[:,i] = Ex2
    masterEy2[:,i] = Ey2

masterEh2 = np.sqrt(masterEx2**2 + masterEy2**2)
maxes1 = np.max(masterEh, axis = 0)
maxes2 = np.max(masterEh2, axis = 0)

clf()

plot(maglat, maxes1)
plot(maglat2, maxes2/1000.)

show()



lat_bins2 = sorted(list(set(maglat2)))
lat_maxh2 = np.zeros(len(lat_bins2))
lat_maxb2 = np.zeros(len(lat_bins2))
lat_maxdb2 = np.zeros(len(lat_bins2))

lat_maxE, lat_maxB, lat_maxdB = [], [], []

maxesB2 = np.max(masterBh2, axis = 0)
maxesE2 = np.max(masterEh2, axis = 0)
maxesdB2 = np.max(masterdBh2, axis = 0)

for i, v in enumerate(lat_bins2):
    ind = maglat2==v
    lat_maxB.append(np.max(maxesB2[ind]))
    lat_maxE.append(np.max(maxesE2[ind]))
    lat_maxdB.append(np.max(maxesdB2[ind]))


clf()

ax1 = subplot(311)
title("2012 Event", fontsize = 24)

plot(lat_bins, lat_maxh, label = "2012 Event, Inputs * 0.5")
plot(lat_bins2, array(lat_maxE)/1000., label = "2012 Event")
legend()
ylabel("Eh (V/km)", fontsize = 20)
grid(True)

ax2 = subplot(312)
plot(lat_bins, lat_maxb, label = "2012 Event, Inputs * 0.5")
plot(lat_bins2, lat_maxB, label = "2012 Event")
legend()
ylabel("Bh (nT)", fontsize = 20)
grid(True)

ax3 = subplot(313)
plot(lat_bins, lat_maxdb, label = "2012 Event, Inputs * 0.5")
plot(lat_bins2, lat_maxdB, label = "2012 Event")
legend()
ylabel("dBh (nT/min)", fontsize = 20)
xlabel("Latitude", fontsize = 20)

grid(True)
show()

filename = "/home/blake/Desktop/Sean_Blake_081318_1_latoutput.txt"
np.savetxt(filename, list(zip(lat_bins2, array(lat_maxE)/1000., lat_maxB, lat_maxdB)))

filename = "/home/blake/Desktop/Sean_Blake_092818_1_latoutput.txt"
np.savetxt(filename, list(zip(lat_bins, lat_maxh, lat_maxb, lat_maxdb)))





















