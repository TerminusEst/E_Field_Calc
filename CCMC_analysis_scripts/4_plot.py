from draco import EField as ef
from draco import functions as eff
from scipy.interpolate import griddata

################################################################################
datafolder = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_032919_2/Numpy_Data/"
data = np.load(datafolder + "IMAG_E.npy")

if True:
    masterEx = np.load(datafolder + "EX.npy")
    masterEy = np.load(datafolder + "EY.npy")
    masterEh = np.sqrt(masterEx**2 + masterEy**2)
    
    lonlat = np.load(datafolder + "LONLAT.npy")
    lon, lat = lonlat[:,0], lonlat[:,1]

    td = np.load(datafolder + "TD.npy")

    unique_lats = np.unique(lat)[3:-3] # shave off polar points
    max_Eh = []
    for i in unique_lats:
        ind = lat == i
        temp = masterEh.T[ind]
        max_Eh.append(np.nanmax(temp))
   

td_swmf = np.load("/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_032919_2/Numpy_Data/TD.npy")

IMAG_td = datetime.datetime(2001, 3, 30)
#blank, IMAG_sites = eff.read_in_single_day(IMAG_td)

td_esk = IMAG_sites['ESK'].timedate
eh_esk = IMAG_sites['ESK'].eh






clf()

ax1 = subplot(211)
title("2003-11-20\nReal Dst: -422 nT\nSWMF Dst: -356 nT", fontsize = 24)
plot(unique_lats, max_Eh, label = "Max SWMF Data")
plot(data[:,0], data[:,1], 'o', label = "IMAG Data")
plot(data[:,2], data[:,3], 'o', label = "Interpolated SWMF Data")
legend()
ylabel("Max Eh (mV/km)", fontsize = 24)
ax1.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax1.tick_params(axis = 'both', which = 'minor', labelsize = 16)

ax2 = subplot(212)

plot(unique_lats, max_Eh, label = "Max SWMF Data")
plot(data[:,0], data[:,1], 'o', label = "IMAG Data")
plot(data[:,2], data[:,3], 'o', label = "Interpolated SWMF Data")
legend()
ylabel("Max Eh (mV/km)", fontsize = 24)
ax2.semilogy()
xlabel("Magnetic Latitude", fontsize = 24)


ax2.tick_params(axis = 'both', which = 'major', labelsize = 16)
ax2.tick_params(axis = 'both', which = 'minor', labelsize = 16)

show()










































