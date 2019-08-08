#imports
import numpy as np
import astropy.units as u
%matplotlib inline
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Use the Hubble ETC to find the exposure times to get an SNR of 3 for a single filter for magnitudes 23-28, every half-magnitude.

magnitudes = np.arange(23,28.5,0.5) #ab
exposure_time = np.array([18.53,29.56,47.33,76.25,124.01,204.69,345.80,605.24,1115.44,2201.54, 4692.60]) #s

print(magnitudes)
print(exposure_time)

#take the log
log_exp = np.log(exposure_time)

#fit a curve to the data
def func(x, a, b, c):
    return a*x**2 + b*x + c

xdata = magnitudes
ydata = log_exp

plt.scatter(xdata, ydata, label='data')


popt, pcov = curve_fit(func, xdata, ydata)
popt

#plot the data + the curve fit and print out the coefficients

plt.plot(xdata, func(xdata, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))


plt.xlabel('magnitude')
plt.ylabel('log of exposure time')
plt.legend()
plt.show()

print("a:", a, "b:", b, "c:", c)

#Calculate the limiting magnitude of an image in that filter given an exposure time

exptime = 2148

a = popt[0]
b = popt[1]
c = popt[2] - np.log(exptime)

lim_mag = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)


print("The limiting magnitude is:",lim_mag)
