#imports
import pandas as pd
import numpy as np
import astroquery
import astropy
from astroquery.mast import Observations


quasar_list = pd.read_csv('/Users/lcgordon/Downloads/5000-hmq') #import list
n = 0 #counts through coordinate list
m = 0 #counter for coordinate pairs WITH DATA

#find all of the 
while n < 5000:
    #produce coordinate string from list
    RA_float = quasar_list.iloc[n,0]
    RA = str(RA_float)
    DEC_float = quasar_list.iloc[n,1]
    DEC = str(DEC_float)
    coords = RA + "_" + DEC
    
    #query MAST, produce table of all data products associated with those coordinates
    table_of_observations = Observations.query_region(coords, radius = 0.02)
    product_list = Observations.get_product_list(product_list) #makes list of all products
    
    i = 0 #counter for obsTable rows
    
    while i < len(table_of_observations): 
        obs_collection = str(table_of_observations[i][1]) #converts the type to a string (essential)
        if  obs_collection == "HST": #check to see if the observatory collection is HST, caNNOT USE 'and' comparison
            m = m + 1 #add one to number of coordinate pairs that have data
            print("Data at", coords) #tells you which coordinate pairs, optional
            i = len(obsTable) + 1 #makes i out of the range of obstable, ending this while loop
            
            #Remove the comments on the following lines to download the data for the coordinates with HST data.
            downloads = Observations.download_products(product_list, download_dir=coords,
                                            obs_collection = "HST",
                                            dataproduct_type = "image",
                                            productType = "SCIENCE") 
            
        else: #if not HST, add one to i to continue looping through obsTable
            i = i + 1
    
    n = n+1 #move to next set of coordinate pairs
    print(n)
print("m =",m, "n=", n, "i =", i)
