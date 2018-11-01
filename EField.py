import os
import numpy
import datetime
import aacgmv2

##########################################################################################

def Time2Float(x):
	"""Converts datetime to float, so that interpolation/smoothing can be performed"""
	if (type(x) == np.ndarray) or (type(x) == list):
		emptyarray = []
		for i in x:
			z = (i - datetime.datetime(1970, 1, 1, 0)).total_seconds()
			emptyarray.append(z)
		emptyarray = np.array([emptyarray])
		return emptyarray[0]
	else:
		return (x - datetime.datetime(1970, 1, 1, 0)).total_seconds()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def Float2Time(x):
	"""Converts array back to datetime so that it can be plotted with time on the axis"""
	if (type(x) == np.ndarray) or (type(x) == list):
		emptyarray = []
		for i in x:
			z = datetime.datetime.utcfromtimestamp(i)
			emptyarray.append(z)
		emptyarray = np.array([emptyarray])
		return emptyarray[0]
	else:
		return datetime.datetime.utcfromtimestamp(x)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

def Model_Profiles(mystr):
    """Return Quebec or British Columbia 1D resistivity models. 
    From Boteler & Pirjola (1998), 'The complex-image method for calculating the
    magnetic and electric fields produced at the surface of the Earth by the 
    auroral electrojet'
    
    Parameters
    -----------
    mystr = string. "Q" for Quebec, "BC" for British Columbia

    Returns
    -----------
    resistivities = array of resistivity values in Ohm.m
    thicknesses = array of thicknesses in m
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
    """Calculate 1D Z-Tensor for given ground resistivity profile.

    Parameters
    -----------
    resistivities = array or list of resistivity values in Ohm.m

    thicknesses = array or list of thicknesses in m.
        **len(resistivities) must be len(thicknesses) + 1**

    frequencies = array or list of frequencies to get response of
    
    Returns
    -----------
    Z = complex array of Z tensor values
    
    Taken from:
    http://www.digitalearthlab.com/tutorial/tutorial-1d-mt-forward/  
    """
    if len(resistivities) != len(thicknesses) + 1:
        print("Length of inputs incorrect!")
        return 
    
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
    
def E_Field_1D(bx, by, resistivities, thicknesses, timestep = 60., Z = None, calc_Z = True):
    """Calculate horizontal E-field components given Bx, By, resistivities and thicknesses.
    
    Parameters
    -----------
    bx, by = array of Bx, By timeseries in nT

    resistivities = array or list of resistivity values in Ohm.m

    thicknesses = array or list of thicknesses in m.
        **len(resistivities) must be len(thicknesses) + 1**

    timestep = time between samples (default is 60. for minute sampling)
    
    Z = complex Z-tensor array. If not supplied, Z will be calculated from input
        resistivities and thicknesses
    
    Returns
    -----------
    ext, eyt = arrays of electric field components in mV/km
    """
    mu0 = 4*np.pi * 1e-7
    freq = np.fft.fftfreq(bx.size, d = timestep)
    freq[0] = 1e-100

    if calc_Z == True:  # if you need to calculate Z
        Z = Z_Tensor_1D(resistivities, thicknesses, freq)
        
    bx_fft = np.fft.fft(bx)
    by_fft = np.fft.fft(by)

    exw = Z * by_fft/mu0; 
    eyw = -1 * Z * bx_fft/mu0

    ext = 1e-3 * np.fft.ifft(exw).real
    eyt = 1e-3 * np.fft.ifft(eyw).real

    return ext, eyt
    
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class Intermag:
    """Class for reading in INTERMAGNET magnetic data files easily.
    Can clean the data (crudely), merge multiple datafiles and calculate geoelectric fields
        
	Parameters
	-----------
	filename = location of INTERMAGNET file
	clean = set to True to clean data using .clean_data() (default = False)


    Attributes
    -----------    
    name = 3-letter string (e.g., VAL)
    glat, glon = geographic latitude and longitude
    mlat, mlon = aagcmv2 calculated magnetic latitude and longitude
    data_type = type of data in INTERMAGNET file (e.g., XYZF or HDZF)
    bx, by, bz = numpy arrays of raw data 
    bh = calculated horizontal component
    dbh = rate of change of bh
    timedate = array of datetime objects
    timedate_float = array of datetime floats (i.e., seconds from 1970)
        
    **IF YOU RUN .clean_data()**
    bx_clean, by_clean, bz_clean = cleaned bx, by, bz components
    bh_clean = cleaned horizontal component
    dbh_clean = rate of change of cleaned horizontal component
    good_data = boolean array indicating good or bad data
    
    **IF YOU RUN .calculate_efield(resistivities, thicknesses)**
    ex, ey, eh = calculated horizontal electric field components from cleaned bx, by

    Functions
    -----------    
    .clean_data() = cleans the data in place
    .calculate_efield() = calculates horizontal electric fields in mV/km
    
    Example use
    ----------- 
    >>> # read in first and second days of data, cleaning the data for each day
    >>> day1 = Intermag(filename1, clean = True)
    >>> day2 = Intermag(filename2, clean = True)
    >>> # merge the data
    >>> day1.merge_days(day2)
    >>> # clean the data with different std:
    >>> day1.clean_data(standard_deviations = 3, absolute_diff = 5000)
    >>> # calculate surface electric fields for 1D resistivity profile
    >>> day1.calculate_efield(resistivities, thicknesses)
    
	-----------------------------------------------------------------
    """
    def __init__(self, filename, clean = False):
        f = open(filename, 'r')
        data = f.readlines()
        f.close()

        for index, value in enumerate(data):
            if "IAGA CODE" in value:    # name of site
                i = value.split(" ")
                self.name = [x for x in i if x != ""][2]
            if "Geodetic" in value:     # geographic latitude
                if "Latitude" in value:
                    j = value.split(" ")
                    self.glat = float([x for x in j if x != ''][2])
                if "Longitude" in value: # geographic longitude
                    k = value.split(" ")
                    self.glon = float([x for x in k if x != ''][2])
            if ("Reported" in value) and ("#" not in value):
                l = value.split(" ")    # data type (XYZ or HDZ)
                self.data_type = [x for x in l if x != ''][1]
            if "DATE" in value: # start of data in file
                skiprows = index + 1
                break

        # read in the actual data
        timedate1, bx1, by1, bz1 = [], [], [], []
        for index, value in enumerate(data[skiprows:]):
            split_line = value.split(" ")
            split_line_no_space = [x for x in split_line if x != ""]

            # get list of datetimes
            dates = split_line_no_space[0]
            year, month, day = int(dates[:4]), int(dates[5:7]), int(dates[8:])
            times = split_line_no_space[1]
            hour, minute, second = int(times[:2]), int(times[3:5]), int(times[6:8])
            timedate1.append(datetime.datetime(year, month, day, hour, minute, second))

            # get bx, by, bz
            bx1.append(float(split_line_no_space[3]))
            by1.append(float(split_line_no_space[4]))
            bz1.append(float(split_line_no_space[5]))

        # convert from HDZ to XYZ if needed
        if self.data_type == "HDZF":
            H, D = np.array(bx1), np.array(by1)
            D = D * (pi/180.0) / 60.0 #convert D to degrees
            bx1 = cos(D)*H
            by1 = sin(D)*H

        self.timedate = np.array(timedate1)
        self.bx, self.by, self.bz = np.array(bx1), np.array(by1), np.array(bz1)
        self.bh = np.sqrt(self.bx**2 + self.by**2)
        self.dbh = np.gradient(self.bh)
        self.timedate_float = Time2Float(self.timedate)

        # get magnetic coordinates
        mag_pos = aacgmv2.wrapper.convert_latlon_arr(np.array([self.glat]), np.array([self.glon]), np.array([100]), self.timedate[0], code = "G2A")
        self.mlat = mag_pos[0][0]
        self.mlon = mag_pos[1][0]
        
        if clean == True:
            self.clean_data()
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
        
    def clean_data(self, standard_deviations = 8, absolute_diff = 10000):
        """crude cleaning function
        any points > standard_deviations away from mean
        or any points > absolute_diff from mean == removed
        a simple simple linear interpolation is then performed
        """
        
        good_data = np.ones(len(self.bx))   # assume all data good at start
        # now figure out which points to discard
        for component in [self.bx, self.by, self.bz]:
            good_std = np.abs(np.mean(component) - component) < standard_deviations * np.std(component)
            good_abs = np.abs(component - np.mean(component)) < 10000
            good_data = good_data * good_std * good_abs
        
        good_data = good_data.astype(bool)
        self.good_data = good_data
 
        # if all points are bad, tell the user, and give arrays of 0s.
        if sum(good_data) == 0:
            print(self.name, " is all bad data!")
            len_data = len(self.bx)
            self.bx_clean, self.by_clean, self.bz_clean = np.zeros(len_data), np.zeros(len_data), np.zeros(len_data)
            self.bh_clean, self.dbh_clean = np.zeros(len_data), np.zeros(len_data)
            
        # otherwise, interpolate the data
        self.bx_clean = np.interp(self.timedate_float, self.timedate_float[good_data], self.bx[good_data])
        self.by_clean = np.interp(self.timedate_float, self.timedate_float[good_data], self.by[good_data])
        self.bz_clean = np.interp(self.timedate_float, self.timedate_float[good_data], self.bz[good_data])
        self.bh_clean = np.sqrt(self.bx_clean**2 + self.by_clean**2)
        self.dbh_clean = np.gradient(self.bh_clean)

        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
     
    def merge_days(self, temp):
        """Merge the timeseries of two Intermag class instances
        If your two Intermag instances are day1, day2:
        >>> day1.merge_days(day2)
        """
        
        labels = ["bx", "by", "bz", "bh", "dbh", "timedate", "timedate_float", "good_data",
                "bx_clean", "by_clean", "bz_clean", "bh_clean", "dbh_clean"]
        
        for label in labels:
            concat_data = np.concatenate((getattr(self, label), getattr(temp, label)))
            setattr(self, label, concat_data)

        
    def calculate_efield(self, resistivities, thicknesses, timestep = 60., Z = None, calc_Z = True):
        """Calculate horizontal electric field using 1D resistivity profile in mV/km
        If your Intermag instance is day1:
        >>> Qres, Qthick = Model_Profiles("Q") # get Quebec resistivity profile
        >>> day1.calculate_efield(Qres, Qthick)
        """
        
        self.ex, self.ey = E_Field_1D(self.bx_clean, self.by_clean, resistivities, thicknesses, timestep,  Z)
        self.eh = np.sqrt(self.ex**2 + self.ey**2)

    
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# get quebec 
Qres, Qthick = Model_Profiles("Q")
freq = np.fft.fftfreq(4320, d = 60) 
Ztens = Z_Tensor_1D(Qres, Qthick, freq)

folder = "/home/blake/Drive/NASA_Shite/Storm_Analysis/INTERMAGNET/2003_10_29/"
sitename = "aae"
a = Intermag(folder + sitename + "20031029dmin.min", clean = True)
b = Intermag(folder + sitename + "20031030dmin.min", clean = True)
c = Intermag(folder + sitename + "20031031dmin.min", clean = True)

a.merge_days(b)
a.merge_days(c)

a.calculate_efield(Qres, Qthick)

    
