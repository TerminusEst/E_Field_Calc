from draco import EField as ef
from draco import functions as eff
from scipy.interpolate import griddata

################################################################################
datafolder = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_2/Numpy_Data/"

# Read in the data and calculate E-Field
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
   
   
# Read in the INTERMAGNET sites and interpolate the grid data on it
# Interpolate the Time-Series onto a grid of  
# INTERMAG sites
if True:    
    IMAG_td = datetime.datetime(2001, 10, 28)
    blank, IMAG_sites = eff.read_in_single_day(IMAG_td)
    
    IMAG_lats = np.array([IMAG_sites[x].mlat for x in IMAG_sites])
    IMAG_lons = np.array([IMAG_sites[x].mlon for x in IMAG_sites])
    IMAG_lons += (IMAG_lons < 0).astype(int) * 360
    IMAG_lat_scal, IMAG_Eh_scal = eff.lat_scaling(IMAG_sites)

    interpEh = np.zeros((len(masterEh), len(IMAG_lons)))

    for i, v in enumerate(masterEh):
        ehh = np.nan_to_num(v)
        points = list(zip(lon, lat))
        otherpoints = list(zip(IMAG_lons, IMAG_lats))

        eh_interp = griddata(points, ehh, (otherpoints), method = "cubic")
        interpEh[i] = eh_interp

        print("INTERP: ", i)

    # Now write out the data
    IMAG_Eh_interp = np.max(interpEh, axis = 0)
    IMAG_lats_interp, IMAG_Eh_interp = (list(t) for t in zip(*sorted(zip(IMAG_lats, IMAG_Eh_interp))))
        
    a, b = np.array([IMAG_lat_scal]).T, np.array([IMAG_Eh_scal]).T
    c, d = np.array([IMAG_lats_interp]).T, np.array([IMAG_Eh_interp]).T
    out = np.concatenate((a, b, c, d), axis = 1)
    np.save(datafolder + "IMAG_E", out)


clf()

#scatter(unique_lats, np.log10(max_Eh))
#scatter(IMAG_lat_scal, np.log10(IMAG_Eh_scal))
#scatter(IMAG_lats, np.log10(IMAG_Eh_interp))

scatter(unique_lats, max_Eh)
scatter(IMAG_lat_scal, IMAG_Eh_scal)
scatter(IMAG_lats, IMAG_Eh_interp)

show()







