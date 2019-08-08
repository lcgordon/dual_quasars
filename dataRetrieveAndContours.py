#import all your packages

import pandas as pd
import numpy as np
import astroquery
import astropy
import os
import shutil

from astroquery.mast import Observations
from astroquery.mast import Catalogs

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

from time import sleep

#read in the list of quasars
quasar_list = pd.read_csv('/Users/lcgordon/Desktop/MAST-crossmatches-reduced-files/crossmatched_with_mag')
quasar_list = quasar_list.drop_duplicates()

#create an empty object for later use
distVec = []

m = len(quasar_list) 

for n in np.arange(0,m):
    try:
        #pull the coordinate pairs from the csv file
        RA_float = quasar_list.iloc[n,0]
        RA = str(RA_float)
        DEC_float = quasar_list.iloc[n,1]
        DEC = str(DEC_float)
        coords = RA + "_" + DEC
    
        #search for and download all data products associated with those coordinates
        print("coordinates:", coords, " n:", n)     
        table_of_observations = Observations.query_region(coords, radius = 0.02)
        data_product_list = Observations.get_product_list(table_of_observations)
        downloads = Observations.download_products(data_product_list, download_dir=coords,
                                                 obs_collection = "HST",
                                                 dataproduct_type = "image",
                                                 productType = "SCIENCE")
    
        #Filter out all files that are not drz.fits files (drizzle reduced)--->
        findQSOdrz = 'find . -type f \( -name "*drz.fits" \)'
        out = subprocess.Popen(findQSOdrz, shell=True, stdin=subprocess.PIPE, 
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = out.communicate()
        filelist = stdout.decode().split()
    
        length = len(filelist)
        print("drz files found:", length)
    
        #produce contour plots for all available drz images
        if length > 0:
            imageDirName = coords + "Images"       #making image directory and pulling RA/Dec from file name 
            os.mkdir(imageDirName)
            for i in range(length):
                file_location = filelist[i]
                m = re.search('HST/(.+?)/', file_location)
                if m:
                    file_name = m.group(1)
                    print(file_name)
                p = re.search('HSThits/(.+?)_', file_location)
                if p:
                    RA = p.group(1)
                q = re.search('_(.+?)/mastDownload', file_location)
                if q:
                    DEC = q.group(1)
                world_coords = RA + "_" + DEC
    
                path_name = filelist[i]
                print(path_name)
    
                hdu = fits.open(path_name)[1]       #converting world coords to pixel coords
                wcs = WCS(hdu.header)
                pixCoords=SkyCoord([RA], [DEC], frame='fk5', unit=u.deg).to_pixel(wcs=wcs, mode='all')
    
                if 0< pixCoords[0][0] < hdu.header['NAXIS1'] and 0 < pixCoords[1][0] < hdu.header['NAXIS2']:
                    fig = plt.figure()                        #plotting only quasars that lie in image
                    fig.add_subplot(111, projection=wcs)
                    plt.imshow(hdu.data, origin='lower', cmap='gray',clim=(0,.99))
                    plt.xlabel('RA')
                    plt.ylabel('Dec')
                    plt.grid(color='white', ls='solid')
                    plt.contour(hdu.data, levels=[.1,.5,1,4,7,10,20], cmap='cool', alpha=0.5)
    
                    pixCoords=np.concatenate(pixCoords)            #cropping image to only quasar
                    plt.xlim(pixCoords[0]-50, pixCoords[0]+50)
                    plt.ylim(pixCoords[1]-50, pixCoords[1]+50)
                    plt.savefig(fname = file_name ,bbox_inches='tight')
                    pic_name = file_name + ".png"
                    plt.clf()
                    close()
                    shutil.move(pic_name, imageDirName)            #moving the image into its specific directory
                
                    distVec.append(((pixCoords[0]-hdu.header['CRPIX1'])/hdu.header['NAXIS1']+(pixCoords[1]-hdu.header['CRPIX2'])/hdu.header['NAXIS2'])/2)
        
                else:
                    print("Quasar not in image")
        
            if os.path.isdir(coords) == True:
                shutil.rmtree(coords)               #deletes ALL data to conserve space
        
        else: 
            print("No data")
            if os.path.isdir(coords) == True:           #deleting the folder in case no drz files exist
                shutil.rmtree(coords)
        
        if os.path.isdir(imageDirName) == True and len(os.listdir(imageDirName)) == 0:
            os.rmdir(imageDirName)               #deleting empty image directories
    
        n = n + 1
        sleep(4)
    except (ValueError, TimeoutError): #in case an error occurs, wait a few seconds, delete the folder created.
        print("error occurred")
        n=n+1
        sleep(4)
        if os.path.isdir(coords) == True:
            shutil.rmtree(coords)
