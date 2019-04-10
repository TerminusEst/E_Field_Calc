from draco import functions as eff
from draco import EField as ef
from gcvspline import SmoothedNSpline
import seaborn as sns
sns.set()
################################################################################
datafolder = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_2/Numpy_Data/"

data = np.load(datafolder + "IMAG_E.npy")
IMAG_lat, IMAG_e = data[:,0], data[:,1]
SWMF_lat, SWMF_e = data[:,2], data[:,3]

################################################################################
################################################################################

#def ftran_spline_single_day(lats, maxe, bootstrap_num = 500, lat_cutoff = 25):
"""RETURNS: 
0: timedate_float, 1: min_DST, 2: absolute threshold, 3 :absolute gradient
4: bootsrap threshold, 5: bootstrap std, 6: bootstrap lower conf, 7: bootstrap upper conf
8: yval thresh, 9: yval std, 10: yval lower conf, 11: yval upper conf
12: gradients mean, 13: gradients std, 14: gradients lower conf, 15: gradients upper conf
"""

################################################################################
################################################################################

xthresh1, ythresh1, gradients1, AAA = eff.fit_spline_lat_e(IMAG_lat, IMAG_e)
IMAG_x = np.mean(xthresh1)
IMAG_lci = np.quantile(xthresh1, 0.025)
IMAG_uci = np.quantile(xthresh1, 0.975)

xthresh2, ythresh2, gradients2, BBB = eff.fit_spline_lat_e(SWMF_lat, SWMF_e)
SWMF_x = np.mean(xthresh2)
SWMF_lci = np.quantile(xthresh2, 0.025)
SWMF_uci = np.quantile(xthresh2, 0.975)

SWMF_y = np.mean(ythresh2)
SWMF_ylci = np.quantile(ythresh2, 0.025)
SWMF_yuci = np.quantile(ythresh2, 0.975)

print(IMAG_x, IMAG_uci, IMAG_lci, "\n")

print(SWMF_x, SWMF_uci, SWMF_lci)
print(10**SWMF_y, 10**(SWMF_yuci), 10**(SWMF_ylci))

################################################################################

clf()

ax1 = subplot(211)
plot(IMAG_lat, IMAG_e, 'o', color = "red")
plot(SWMF_lat, SWMF_e, 'o', color = "blue")
plot(AAA[0], 10**AAA[1], color = "red")
plot(BBB[0], 10**BBB[1], color = "blue")
#ax2.semilogy()

axvline(IMAG_x, color = 'red', lw = 3)
#axvline(IMAG_lci, color = "red", lw = 3, linestyle = "dashed")
#axvline(IMAG_uci, color = "red", lw = 3, linestyle = "dashed")

axvline(SWMF_x, color = 'blue', lw = 3)
#axvline(SWMF_lci, color = "blue", lw = 3, linestyle = "dashed")
#axvline(SWMF_uci, color = "blue", lw = 3, linestyle = "dashed")
ax1.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax1.tick_params(axis = 'both', which = 'minor', labelsize = 16)

    
ax2 = subplot(212)
plot(IMAG_lat, IMAG_e, 'o', color = "red")
plot(SWMF_lat, SWMF_e, 'o', color = "blue")
plot(AAA[0], 10**AAA[1], color = "red")
plot(BBB[0], 10**BBB[1], color = "blue")
ax2.semilogy()

axvline(IMAG_x, color = 'red', lw = 3)
#axvline(IMAG_lci, color = "red", lw = 3, linestyle = "dashed")
#axvline(IMAG_uci, color = "red", lw = 3, linestyle = "dashed")

axvline(SWMF_x, color = 'blue', lw = 3)
#axvline(SWMF_lci, color = "blue", lw = 3, linestyle = "dashed")
#axvline(SWMF_uci, color = "blue", lw = 3, linestyle = "dashed")

ax2.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax2.tick_params(axis = 'both', which = 'minor', labelsize = 16)


show()
















