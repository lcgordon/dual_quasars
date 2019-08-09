#import your stuff
import pandas as pd
import numpy as np
import astroquery
import astropy
import os
import shutil

from astroquery.mast import Observations

%matplotlib inline
from pylab import *
import matplotlib
import matplotlib.pyplot as plt

from astropy.io import ascii
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

import os, fnmatch
import subprocess
from subprocess import call
import re
from numpy import arange

from astropy import units as u
from astropy import cosmology, coordinates
from astropy.coordinates import Distance
from astropy.cosmology import WMAP5
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import z_at_value

import scipy.integrate as integrate
from scipy import stats

#set up your universe
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
Omega_M = 0.3
Omega_L = 0.7
H0 = 70 * u.km / u.s / u.Mpc
c = 300000 * u.km / u.s

#calculate the solid angles of each of the three main detectors used
#solid angle for wfc3 in the UVIS channel
wfc3_uvis_arcsec = 160* 160 * u.arcsec**2
wfc3_uvis = wfc3_uvis_arcsec.to(u.sr)

#solid angle for wfc3 in the NIR channel
wfc3_nir = 123 * 137 * u.arcsec**2
wfc3_nir = wfc3_nir.to(u.sr)

#solid angle for ACS wfc
wfc_acs = 202 * 202 * u.arcsec**2
wfc_acs = wfc_acs.to(u.sr)

#import the parameters you calculated to be used in finding the limiting magnitude of the image
params = pd.read_csv("/Users/lcgordon/Desktop/params_by_filter.csv")
print(params)

#import the list of quasars with coordinates, redshift, and R magnitude
quasar_list = pd.read_csv('/Users/lcgordon/Desktop/HMQ_with_mag')
print(quasar_list)

n =0
sum_veff = 0
m = len(quasar_list)

for n in np.arange(n,m):
    #create a list in which all of the available volumes for the target will be stored
    available_volumes = []
    try:
    #pull quantities from the list
        RA = quasar_list.iloc[n][0]
        RA = str(RA)
        DEC = quasar_list.iloc[n][1]
        DEC = str(DEC)
        coords = RA + "_" + DEC
        redshift_target = quasar_list.iloc[n][3] #redshift of the target.
        r_mag_target = quasar_list.iloc[n][2] #magnitude of target
    
    #calculate the absolute magnitude of the target from the red magnitude
        dL = Distance(unit = u.Mpc, z=redshift_target, cosmology = cosmo) / u.Mpc
        magnitude_target = r_mag_target - 5*np.log(dL) - 25

#print(r_magnitude_target, magnitude_target, redshift_target)

    #search for and download all data products associated with those coordinates    
        obsTable = Observations.query_region(coords, radius = 0.02)
        dataProductsByObservation = Observations.get_product_list(obsTable)
        manifest = Observations.download_products(dataProductsByObservation, download_dir=coords,
                                                     obs_collection = "HST",
                                                     dataproduct_type = "image",
                                                     productType = "SCIENCE")
        
        
    except TimeoutError:
        n = n + 1
        print("a timeout error occurred at", coords)
        
       #isolate drizzle files 
    findQSOdrz = 'find . -type f \( -name "*drz.fits" \)'
    out = subprocess.Popen(findQSOdrz, shell=True, stdin=subprocess.PIPE, 
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = out.communicate()
    filelist = stdout.decode().split()
    
    file_number = len(filelist)
    print(filelist)
    
    #calculate the available volume for each drizzle image in this block
    i = 0
    while i < file_number:
        path_name = filelist[i]
        print(path_name)
  
        image = fits.open(path_name,ignore_missing_end=True)
        instrument = image[0].header['INSTRUME']
        exptime = image[0].header['EXPTIME'] #seconds
        detector = image[0].header['DETECTOR']
        try:
             filter = image[0].header["FILTER"]
        except KeyError:
            filter = image[0].header["FILTER1"]
    
        print(instrument, exptime, filter)
        i = i + 1
    
        zero_point = 0
    
    #choose the correct parameters. currently only has parameters for the most commonly used filters. also sets omega.
        if instrument == 'WFC3':
            if filter == "F350LP":
                zero_point = 26.84 #ab mag
                omega = wfc3_uvis
                param0 = params.iloc[6][2]
                param1 = params.iloc[6][3]
                param2 = params.iloc[6][4]
            elif filter == "F814W":
                zero_point = 24.99
                omega = wfc3_uvis
                param0 = params.iloc[6][2]
                param1 = params.iloc[6][3]
                param2 = params.iloc[6][4]
            elif filter == "F475W":
                zero_point = 25.59
                omega = wfc3_uvis
                param0 = params.iloc[7][2]
                param1 = params.iloc[7][3]
                param2 = params.iloc[7][4]
            elif filter == "F606W":
                zero_point = 25.98
                omega = wfc3_uvis
                param0 = params.iloc[8][2]
                param1 = params.iloc[8][3]
                param2 = params.iloc[8][4]
            elif filter == "F438W":
                zero_point =24.73
                omega = wfc3_uvis
                param0 = params.iloc[9][2]
                param1 = params.iloc[9][3]
                param2 = params.iloc[9][4]
            elif filter == "F160W":
                zero_point = 25.94
                omega = wfc3_nir
                param0 = params.iloc[3][2]
                param1 = params.iloc[3][3]
                param2 = params.iloc[3][4]
            elif filter == "F140W":
                zero_point = 26.45
                omega = wfc3_nir
                param0 = params.iloc[4][2]
                param1 = params.iloc[4][3]
                param2 = params.iloc[4][4]
            elif filter == "F110W":
                zero_point = 26.82
                omega = wfc3_nir
                param0 = params.iloc[5][2]
                param1 = params.iloc[5][3]
                param2 = params.iloc[5][4]
            else:
            #if the filter is not one of the ones listed I chose a default parameter set that is about the average of all the calculated ones
                #obviously we would ideally have the parameters for every single filter.
                zero_point = 25.98
                omega = wfc3_uvis
                param0 = params.iloc[8][2]
                param1 = params.iloc[8][3]
                param2 = params.iloc[8][4]
                
        elif instrument == "ACS":
            omega = wfc_acs
            if filter == "F775W":
                zero_point = 25.66
                param0 = params.iloc[1][2]
                param1 = params.iloc[1][3]
                param2 = params.iloc[1][4]
            elif filter == "F555W":
                zero_point = 25.711
                param0 = params.iloc[0][2]
                param1 = params.iloc[0][3]
                param2 = params.iloc[0][4]
                
            else: 
                zero_point = 25.711
                param0 = params.iloc[0][2]
                param1 = params.iloc[0][3]
                param2 = params.iloc[0][4]
        else: 
            print("cannot ID instrument", detector, instrument)
    
    #calculate the limiting magnitude
        a = param0
        b = param1
        c1 = param2 - np.log(exptime)

        lim_mag = (-b + np.sqrt(b**2 - 4*a*c1)) / (2*a)
        print("The limiting magnitude is:",lim_mag)
    
    #find the distance to the object

        dist_lim = 10**((lim_mag - magnitude_target + 5) / 5) * u.parsec
        dist_L = dist_lim.to(u.Mpc)
        #print(dist_L)

    #find z of the limit 
        z_lim = z_at_value(cosmo.luminosity_distance, dist_L, zmin = 0, zmax = 10000000)
            
        z = z_lim #redshift
        if z > 7.5: #7.5 is the maximum redshift of known quasars. sometimes the program spits out really absurd numbers for
        #redshift and so this corrects it back down to 7.5
            z = 7.5
        print(z)
        
        
            #calculate the pieces of dv/dz
        dist_L = Distance(unit = u.Mpc, z=z, cosmology = cosmo) #Mpc 
            
        Denom = np.sqrt(Omega_L + Omega_M * (1+z)**3)

        dist_A = dist_L/((1+z)**2) #Mpc

        dist_H = c/H0 #Mpc

    #finds dv/dz in gpc^3 per steradians
        dv_dz = (dist_H * (1 + z)**2 * dist_A**2 / Denom) * (4 * np.pi / u.sr) #Mpc**3 per sr
        dv_dz_gpc_3 = dv_dz.to(u.Gpc**3/u.sr) #converts to gpc^3 per sr
        
        dv_dz_unitless = dv_dz_gpc_3 * (u.sr/u.Gpc**3) #removes the units so that the integral will happen


    #print("Luminosity distance is: ", dist_L)
    #print("Angular diameter distance is:", dist_A)
    #print("Hubble distance is:", dist_H)
        print("Dv/dz is:", dv_dz_gpc_3)

#integrate to get available volume
        V_eff = integrate.quad(lambda z: dv_dz_unitless, 0, z) * omega * u.Gpc**3 /u.sr
        print("The available volume is", V_eff[0])
        available_volumes.append(V_eff[0]) #adds the volume to the list for this target
    print("The largest available volume for this target is:", max(available_volumes))
    sum_veff = sum_veff + max(available_volumes) #adds the maximum volume for this target to the sum
    print("the total available volume is", sum_veff)    
    n = n + 1 #moves to the next coordinate
    print("N:",n) #prints n to keep track if it shuts down.
    if os.path.isdir(coords) == True: #deletes data
        shutil.rmtree(coords)
