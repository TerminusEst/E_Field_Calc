import os
import math
import cmath
from draco import draco as dr
from draco import plane as pl
import numpy.fft.fftpack as fft

from mpl_toolkits.basemap import Basemap
import matplotlib.patheffects as path_effects
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata

execfile("/home/blake/Drive/PaperMk3/Frequency_Analysis/MultiEfield_pre4.py")
################################################################################
################################################################################
# Read in Conductances

folder = "/home/blake/Drive/PaperMk3/MultiFreq/Conductances/Smooth/"
files = ["300m_SMOOTH_RHO.txt", "1km_SMOOTH_RHO.txt", "3km_SMOOTH_RHO.txt", "10km_SMOOTH_RHO.txt", "30km_SMOOTH_RHO.txt", "60km_SMOOTH_RHO.txt", "100km_SMOOTH_RHO.txt", "200km_SMOOTH_RHO.txt", "400km_SMOOTH_RHO.txt"]

resistances = np.ones((3600, len(files)))

for index, value in enumerate(files):
	res = np.loadtxt(folder + value)[:,2]
	
	resistances[:,index] = res	

################################################################################
# Read In B-fields

folder = "/home/blake/Drive/PaperMk3/FINAL/Single_Freq/DEC2015/1SECS/Output/"
filename = folder + "0001.txt"

lat, lon, bx, by = np.loadtxt(filename, usecols = (0, 1, 2, 3), skiprows = 0, unpack = True)
mega_bx, mega_by = [], []
mega_lats, mega_lons = [], []

for x, y, t, n in zip(bx, by, lat, lon):
	mega_bx.append([x])
	mega_by.append([y])
	mega_lats.append([t])
	mega_lons.append([n])


for i in range(2, 4321, 1):
	filename = "/home/blake/Drive/PaperMk3/FINAL/Single_Freq/DEC2015/1SECS/Output/" + "%04d" % (i) + ".txt"
	lat, lon, bx, by = np.loadtxt(filename, usecols = (0, 1, 2, 3), skiprows = 0, unpack = True)
	
	for index, value in enumerate(bx):
		mega_bx[index].append(bx[index])
		mega_by[index].append(by[index])
		mega_lats[index].append(lat[index])
		mega_lons[index].append(lon[index])

	print "READING: ", i

################################################################################
# Calculate E-field

maxxx = []
maxxy = []
#for jjj in [("5", 300., 70.), ("6", 250., 70.), ("7", 250, 100), ("8", 200, 100)]:


mega_ZZ1, mega_ZZ2 = [], []
mega_ex1, mega_ey1 = [], []
mega_ex2, mega_ey2 = [], []

for i in range(3600):

	bx = array(mega_bx[i]) * 1e-9
	by = array(mega_by[i]) * 1e-9

	signal = array(bx)
	n = signal.size
	timestep = 60.0

	freq = np.fft.fftfreq(n, d=timestep)
	freq[0] = 1./(n*timestep)

	res = resistances[i]
	new_resistivities, thicknesses = res_calculator4(res)

	appRes, phases, ZZ = Z_calculator(new_resistivities, thicknesses, freq)
	ex, ey = MultiECalc(ZZ, by, bx)

	mega_ex2.append(ex)
	mega_ey2.append(ey)

    print lat[i], lon[i]
	print "CALCULATING: ", i, len(new_resistivities), len(thicknesses)

################################################################################
# Write Out

for i in range(2880):
	numb = "%04d" % (i + 1,)
	filename = "/home/blake/Drive/PaperMk3/Frequency_Analysis/JOAN_MULTI/Complete/DEC2015_TAKE2/Efield" + numb + ".txt"

	exx = [x[i].real for x in mega_ex2]
	eyy = [x[i].real for x in mega_ey2]
	np.savetxt(filename, zip(lat, lon, exx, eyy))

	print "SAVING: ", filename


maxx = []
maxy = []


for i in mega_ex2:
	maxx.append(max(abs(i)))

for i in mega_ey2:
	maxy.append(max(abs(i)))

maxxx.append(max(maxx) * 1000.)
maxxy.append(max(maxy) * 1000.)










