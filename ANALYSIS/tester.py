###the top two lines are required on my linux machine
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from scipy.optimize import curve_fit #we could import more, but this is what we need




###defining your fitfunction
def func(x, a, b, c):
    return a - b* np.exp(c * x) 
    
    
###OP's data
baskets = X3
scaling_factor = np.array([95.5, 57.7, 1.4, 21.9, 88.8])


###let us guess some start values
initialGuess=[100, 100,-.01]
guessedFactors=[func(x,*initialGuess ) for x in baskets]
###making the actual fit
popt,pcov = curve_fit(func, baskets, scaling_factor,initialGuess)
#one may want to



print(popt)
print(pcov)




###preparing data for showing the fit
basketCont=np.linspace(min(baskets),max(baskets),50)
fittedData=[func(x, *popt) for x in basketCont]
###preparing the figure
fig1 = plt.figure(1)
ax=fig1.add_subplot(1,1,1)
###the three sets of data to plot
ax.plot(baskets,scaling_factor,linestyle='',marker='o', color='r',label="data")
ax.plot(baskets,guessedFactors,linestyle='',marker='^', color='b',label="initial guess")
ax.plot(basketCont,fittedData,linestyle='-', color='#900000',label="fit with ({0:0.2g},{1:0.2g},{2:0.2g})".format(*popt))
###beautification
ax.legend(loc=0, title="graphs", fontsize=12)
ax.set_ylabel("factor")
ax.set_xlabel("baskets")
ax.grid()
ax.set_title("$\mathrm{curve}_\mathrm{fit}$")
###putting the covariance matrix nicely
tab= [['{:.2g}'.format(j) for j in i] for i in pcov]
the_table = plt.table(cellText=tab,
                  colWidths = [0.2]*3,
                  loc='upper right', bbox=[0.483, 0.35, 0.5, 0.25] )
plt.text(250,65,'covariance:',size=12)
###putting the plot
plt.show()
###done
