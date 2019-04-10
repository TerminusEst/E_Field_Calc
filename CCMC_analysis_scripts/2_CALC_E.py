import os
from draco import EField as ef
import aacgmv2
import matplotlib.dates as mdates

#/usr/bin/ffmpeg -r 30 -f image2 -i %03d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 10 -pix_fmt yuv420p 2003_11_20_Scan.mp4

##########################################################################################
# Read in all of the B first
numpysavefolder = "/home/blake/Drive/NASA_Data/CCMC_Runs_With_Changes/outputs/Sean_Blake_040519_2/Numpy_Data/"

# Read in numpy binaries
masterBx = np.load(numpysavefolder + "BX.npy")
masterBy = np.load(numpysavefolder + "BY.npy")
lonlat = np.load(numpysavefolder + "LONLAT.npy")
geolon = lonlat[:,0]
geolat = lonlat[:,1]
timedate1 = np.load(numpysavefolder + "TD.npy")

##########################################################################################
Qres, Qthick = ef.Model_Profiles("Q")
freq = np.fft.fftfreq(masterBx[:,0].size, d = 60)
freq[0] = 1e-100
ZZ = ef.Z_Tensor_1D(Qres, Qthick, freq)

masterEx, masterEy = np.zeros(masterBx.shape), np.zeros(masterBx.shape)

for i, v in enumerate(masterBx.T):
    bx = masterBx[:,i]
    by = masterBy[:,i]
    
    ex, ey = ef.E_Field_1D(bx, by, Qres, Qthick, timestep = 60, Z = ZZ)
    
    masterEx[:,i] = ex
    masterEy[:,i] = ey

    if (i%100) == 0:
        print(i)

##########################################################################################

np.save(numpysavefolder + "EX", masterEx)
np.save(numpysavefolder + "EY", masterEy)
































