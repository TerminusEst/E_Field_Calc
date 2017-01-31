import math
import cmath
import numpy as np
import draco as dr
import datetime

# assume you have data minute_bx, minute_by and Ex and Ey
def res_calculator(old_resistivities, second, last):
	old_resistivities = np.log10(old_resistivities)

	depths_given = [300., 1000., 3000., 10000., 30000.]
	depths_wanted = [700, 2000, 7000, 20000]

	# Find Resistivities from averages
	new_resistivities = np.zeros(len(depths_wanted)+2)
	new_resistivities[0] = old_resistivities[0]

	for index, value in enumerate(depths_wanted):
		a = depths_given[index + 1]	
		b = old_resistivities[index +1]
		c = depths_given[index]	
		d = old_resistivities[index]
		e = depths_wanted[index]

		new_resistivities[index + 1] = ((a*b) - (c*d))/e


	res_new = 10**array(new_resistivities)

	depths_wanted.insert(0, 300.)
	depths_wanted.append(70000.)
	depths_wanted.append(100000.)

	# "floor" resistivity equal to the deepest measured
	res_new[-1] = 100.
	res_new = list(res_new)
	res_new.append(second)
	res_new.append(last)

	return res_new, depths_wanted

def res_calculator4(zzz):
	zzz = log10(zzz)
	depths_given = [300., 1000., 3000., 10000., 30000., 60000., 100000., 200000., 400000.]
	depths_wanted = [300, 700, 2000, 7000, 20000, 30000., 40000., 100000., 200000.]

	output = [zzz[0]]

	for index, value in enumerate(depths_wanted):
		if index == 0:
			continue
		a = zzz[index]
		b = depths_given[index]
		c = zzz[index - 1]
		d = depths_given[index - 1]
		e = depths_wanted[index]

		x = ((a*b) - (c*d))/e

		output.append(x)

	output.append(2)

	output = 10**array(output)
	return output, depths_wanted

def Z_calculator(resistivities, thicknesses, frequencies):
	mu = 4*math.pi*1E-7; #Magnetic Permeability (H/m)
	n = len(resistivities);

	appRes, phases, Zs = [], [], []

	for frequency in frequencies:   
		w =  2*math.pi*frequency;       
		impedances = list(range(n))
		#compute basement impedance
		impedances[n-1] = cmath.sqrt(w*mu*resistivities[n-1]*1j)

		for j in range(n-2,-1,-1):
		    resistivity = resistivities[j]
		    thickness = thicknesses[j]
	  
		    # 3. Compute apparent resistivity from top layer impedance
		    #Step 2. Iterate from bottom layer to top(not the basement) 
		    # Step 2.1 Calculate the intrinsic impedance of current layer
		    dj = cmath.sqrt((w * mu * (1.0/resistivity))*1j)
		    wj = dj * resistivity
		    # Step 2.2 Calculate Exponential factor from intrinsic impedance
		    ej = cmath.exp(-2*thickness*dj);                    
		
		    # Step 2.3 Calculate reflection coeficient using current layer
		    #          intrinsic impedance and the below layer impedance
		    belowImpedance = impedances[j + 1]
		    rj = (wj - belowImpedance)/(wj + belowImpedance)
		    re = rj*ej
		    Zj = wj * ((1 - re)/(1 + re))
		    impedances[j] = Zj

		# Step 3. Compute apparent resistivity from top layer impedance
		Z = impedances[0]

		realZ = real(Z)
		absZ = abs(Z)
		apparentResistivity = (absZ * absZ)/(mu * w)
		phase = math.atan2(Z.imag, Z.real)

		appRes.append(apparentResistivity)
		phases.append(phase)

		Zs.append(Z)

	return appRes, phases, Zs

def MultiECalc(ZZ, minute_bx, minute_by):

	mu = 4*math.pi*1E-7; # Magnetic Permeability (H/m) # 

	by_freq = numpy.fft.fftpack.fft(array(minute_by))
	fff = (1./mu)*array(ZZ)*by_freq
	ex = numpy.fft.fftpack.ifft(fff)

	bx_freq = numpy.fft.fftpack.fft(array(minute_bx))
	fff = (-1./mu)*array(ZZ)*bx_freq
	ey = numpy.fft.fftpack.ifft(fff)

	return ex[1440:], ey[1440:]


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


		time_calc.append(dr.float2time(timedate_float[i - 1]) + datetime.timedelta(minutes = 2))

	return ey_calc, ex_calc, time_calc

def dbdt(bx, by):
	""" Assumes minute values, returns db/dt nt min^-1"""
	dbx, dby = [], []
	for t in range(len(bx)-1):
		dbx.append((bx[t + 1] - bx[t]))
		dby.append((by[t + 1] - by[t]))
	dbx.append(dbx[-1])
	dby.append(dby[-1])

	return dbx, dby














