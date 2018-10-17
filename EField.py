""" Functions for calculating 1d Z tensor and horizontal E-fields from Bx, By
Also for reading in INTERMAGNET magnetic files."""

import numpy as np
import math
import os
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

class Intermag:
    def __init__(self, filename):
        f = open(filename, 'r')
        data = f.readlines()
        f.close()

        for index, value in enumerate(data):
            if "IAGA CODE" in value:
                i = value.split(" ")
                self.name = [x for x in i if x != ""][2]

            if "Geodetic" in value:
                if "Latitude" in value:
                    j = value.split(" ")
                    self.lat = float([x for x in j if x != ''][2])

                if "Longitude" in value:
                    k = value.split(" ")
                    self.lon = float([x for x in k if x != ''][2])

            if "Reported" in value:
                l = value.split(" ")
                self.data_type = [x for x in l if x != ''][1]

            if "DATE" in value:
                skiprows = index + 1
                break

        timedate1, bx1, by1, bz1 = [], [], [], []
        for index, value in enumerate(data[skiprows:]):
            split_line = value.split(" ")
            split_line_no_space = [x for x in split_line if x != ""]

            # get timedate
            dates = split_line_no_space[0]
            year, month, day = int(dates[:4]), int(dates[5:7]), int(dates[8:])
            times = split_line_no_space[1]
            hour, minute, second = int(times[:2]), int(times[3:5]), int(times[6:8])
            timedate1.append(datetime.datetime(year, month, day, hour, minute, second))

            # get bx, by, bz
            bx1.append(float(split_line_no_space[3]))
            by1.append(float(split_line_no_space[4]))
            bz1.append(float(split_line_no_space[5]))

        # convert from HDZ to XYZ
        if self.data_type == "HDZF":
            H, D = np.array(bx1), array(by1)

            D = D * (pi/180.0) / 60.0 #convert D to degrees
            bx1 = cos(D)*H
            by1 = sin(D)*H

        self.timedate = array(timedate1)
        self.bx = array(bx1)
        self.by = array(by1)
        self.bz = array(bz1)
        self.bh = np.sqrt(self.bx**2 + self.by**2)

        self.timedate_float = time2float(self.timedate)

        # clean data of 99999's
        good_data1 = self.bz != 99999.
        good_data2 = self.bh - np.mean(self.bh) < 10000.
        good_data = good_data1|good_data2
 
        self.bad_data = np.sum(self.bz == 99999.)
        if self.bad_data == len(self.bz):
            print(self.name + " ALL DATA BAD")
            self.interp_bx, self.interp_by, self.interp_bz, self.interp_bh, self.dbh, self.eh = np.zeros(len(self.bz)), np.zeros(len(self.bz)), np.zeros(len(self.bz)), np.zeros(len(self.bz)), np.zeros(len(self.bz)), np.zeros(len(self.bz))
            return

        self.interp_bx = np.interp(self.timedate_float, self.timedate_float[good_data], self.bx[good_data])
        self.interp_by = np.interp(self.timedate_float, self.timedate_float[good_data], self.by[good_data])
        self.interp_bz = np.interp(self.timedate_float, self.timedate_float[good_data], self.bz[good_data])
        self.interp_bh = np.sqrt(self.interp_bx**2 + self.interp_by**2)
        self.dbh = np.gradient(self.interp_bh)
        self.eh = np.zeros(len(self.bz))


    def merge_days(self, temp, resistivities, thicknesses):   # merge two intermagnet files
        self.timedate = np.concatenate((self.timedate, temp.timedate))
        self.timedate_float = np.concatenate((self.timedate_float, temp.timedate_float))

        self.bx = np.concatenate((self.bx, temp.bx))
        self.by = np.concatenate((self.by, temp.by))
        self.bz = np.concatenate((self.bz, temp.bz))
        self.bh = np.concatenate((self.bh, temp.bh))
        self.interp_bx = np.concatenate((self.interp_bx, temp.interp_bx))
        self.interp_by = np.concatenate((self.interp_by, temp.interp_by))
        self.interp_bz = np.concatenate((self.interp_bz, temp.interp_bz))
        self.interp_bh = np.concatenate((self.interp_bh, temp.interp_bh))
        self.dbh = np.concatenate((self.dbh, temp.dbh))
        self.bad_data += temp.bad_data

        self.ex, self.ey = E_Field_1D(self.interp_bx, self.interp_by, resistivities, thicknesses, timestep = 60.)
        self.eh = np.sqrt(self.ex**2 + self.ey**2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def Model_Profiles(mystr):
    """Return Quebec or British Columbia model inputs
    """
    if mystr == "BC":
        resistivities = np.array([500., 150., 20., 300., 100., 10., 1.])
        thicknesses = 1000. * np.array([4., 6., 5., 65., 300., 200.])
        
    elif mystr == "Q":
        resistivities = np.array([20000., 200, 1000, 100, 3])
        thicknesses = 1000. * np.array([15, 10, 125, 200])
    else:
        print("Choose Either 'Q' for Quebec or 'BC' for British Columbia")
        return
    return resistivities, thicknesses

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def Z_Tensor_1D(resistivities, thicknesses, frequencies):
    """Calculate 1D Z Tensor for given thicknesses, frequencies
        
    Taken from:
    http://www.digitalearthlab.com/tutorial/tutorial-1d-mt-forward/  
    """
    mu = 4*np.pi*1E-7; #Magnetic Permeability (H/m)
    n = len(resistivities);
    master_Z, master_absZ, master_phase = [], [], []

    for frequency in frequencies:   
        w =  2*np.pi*frequency;       
        impedances = list(range(n));
        #compute basement impedance
        impedances[n-1] = np.sqrt(w*mu*resistivities[n-1]*1j);
       
        for j in range(n-2,-1,-1):
            resistivity = resistivities[j];
            thickness = thicknesses[j];
      
            # 3. Compute apparent resistivity from top layer impedance
            #Step 2. Iterate from bottom layer to top(not the basement) 
            # Step 2.1 Calculate the intrinsic impedance of current layer
            dj = np.sqrt((w * mu * (1.0/resistivity))*1j);
            wj = dj * resistivity;
            # Step 2.2 Calculate Exponential factor from intrinsic impedance
            ej = np.exp(-2*thickness*dj);                     
        
            # Step 2.3 Calculate reflection coeficient using current layer
            #          intrinsic impedance and the below layer impedance
            belowImpedance = impedances[j + 1];
            rj = (wj - belowImpedance)/(wj + belowImpedance);
            re = rj*ej; 
            Zj = wj * ((1 - re)/(1 + re));
            impedances[j] = Zj;    
    
        # Step 3. Compute apparent resistivity from top layer impedance
        Z = impedances[0];
        phase = math.atan2(Z.imag, Z.real)
        master_Z.append(Z)
        master_absZ.append(abs(Z))
        master_phase.append(phase)
        #master_res.append((absZ * absZ)/(mu * w))
    return np.array(master_Z)
  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def E_Field_1D(bx, by, resistivities, thicknesses, timestep = 60.):
    """Calculate horizontal E-field components given Bx, By, resistivities and thicknesses.
    Returns ex, ey in mV/km
    """
    mu0 = 4*np.pi * 1e-7
    freq = np.fft.fftfreq(bx.size, d = timestep)
    freq[0] = 1e-100

    Z = Z_Tensor_1D(resistivities, thicknesses, freq)
    bx_fft = np.fft.fft(bx)
    by_fft = np.fft.fft(by)

    exw = Z * by_fft/mu0; 
    eyw = -1 * Z * bx_fft/mu0

    ext = 1e-3 * np.fft.ifft(exw).real
    eyt = 1e-3 * np.fft.ifft(eyw).real

    return ext, eyt

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def PlaneE(bx, by, timedate_float, cond):
	""" Plane wave method of calculating E-field"""
	dbx, dby = [], []
	for t in range(len(bx)-1):
		dbx.append((bx[t + 1] - bx[t])/60.)
		dby.append((by[t + 1] - by[t])/60.)
	dbx.append(dbx[-1])
	dby.append(dby[-1])


	constant = (-1) * math.sqrt(1/(np.pi * 12.566e-7 * cond))
	L = 1438
	A = np.ones(L)
	B = 1/np.sqrt(np.arange(1438) + 2)#[::-1]

	for index, value in enumerate(A):
		d = index + 1
		if d%2 == 1:
			A[index] += 1

	A = A[::-1]
	ex_calc, ey_calc = [], []
	time_calc = []

	for i in range(1440, len(dby), 1):
		Dx = math.sqrt(60) * (((4/3.) * dby[i]) + dby[i-1])
		Cx = dby[i-1440:i-2][::-1]
		Ex = sum(A * B * Cx) * (2./3.) * math.sqrt(60)
		totalx = constant* (Ex + Dx)
		ex_calc.append(totalx)

		Dy = math.sqrt(60) * (((4/3.) * dbx[i]) + dbx[i-1])
		Cy = dbx[i-1440:i-2][::-1]
		Ey = sum(A * B * Cy) * (2./3.) * math.sqrt(60)
		totaly = constant* (Ey + Dy)
		ey_calc.append(totaly)


		time_calc.append(float2time(timedate_float[i - 1]) + datetime.timedelta(minutes = 2))

	return ey_calc, ex_calc, time_calc

##########################################################################################
# Example of how to read in 
Qres, Qthick = Model_Profiles("Q")
Unires = np.ones(len(Qres))

filename1 = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/2003_11_20/aae20031120dmin.min"
filename2 = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/2003_11_20/aae20031121dmin.min"

a = Intermag(filename1)
b = Intermag(filename2)

a.merge_days(b, Unires, Qthick)


ex1, ey1 = a.ex, a.ey
ey2, ex2, time_calc = PlaneE(a.interp_bx, a.interp_by, a.timedate_float, 1.)

clf()

ax1 = subplot(211)
plot(ex1)
plot(ey1)

ax2 = subplot(212)
plot(ex2)
plot(ey2)

show()
















