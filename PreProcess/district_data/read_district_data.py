# dependencies
import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import seaborn as sns

sns.set()
sns.set_style("dark")
# sns.set_style("whitegrid")
sns.set_context("paper", font_scale=4, rc={"lines.linewidth": 4})
mpl.rcParams['lines.linewidth'] = 3

import pyvista as pv
import geopandas as gpd
from shapely.geometry import Point
from descartes import PolygonPatch

from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

import random
import datetime

def get_one_string(ls):
    msg = ''
    N = len(ls)
    for i in range(N):
        msg += str(ls[i])
        if i < N - 1:
            msg += ' '
            
    return msg


def read_dist_map(fdir, county_nms, county_dist, county_dist_names, county_polys, verbosity = 0):
    
    ## read shapefile
    tx_dist_shp = gpd.read_file(fdir + 'shp/District.shp')
    
    ## extract key data
    dist_geom = tx_dist_shp['geometry']
    dist_nbr = tx_dist_shp['TXDOT_DI_1']
    dist_nm = tx_dist_shp['TXDOT_DI_2']
    
    ## Create district list where 1st element is the 
    ## 1st district, 2nd element is the 2nd district, etc
    N_districts = 25
    new_dist_geom = []
    new_dist_nm = []
    for i in range(N_districts):
        for j in range(N_districts):
            if dist_nbr[j] == i + 1:
                new_dist_nm.append(dist_nm[j])
                new_dist_geom.append(dist_geom[j])
    
    ## collect counties for each district
    new_dist_county = []
    new_dist_county_names = []

    N_counties = len(county_nms)
    for i in range(N_districts):
        # search for counties that have i+1 as their district
        cn = []
        nm = []
        for j in range(N_counties):
            if county_dist[j] == i:
                cn.append(j)
                nm.append(county_nms[j])
        new_dist_county.append(cn)
        new_dist_county_names.append(nm)
    
    ## write to file
    write_to_file = True
    if write_to_file:
        fdname = open(fdir + 'district_names.txt', 'w')
        fdcounty = open(fdir + 'district_counties.txt', 'w')
        for i in range(N_districts):
            #fcpoly.write('{}\n'.format(new_county_polys[i]))
            fdname.write("{}\n".format(new_dist_nm[i]))

            Nc = len(new_dist_county[i])
            for j in range(Nc):
                fdcounty.write('{}'.format(new_dist_county[i][j]))
                if j < Nc - 1:
                    fdcounty.write(' ')
                else:
                    fdcounty.write('\n')

        #fcpoly.close()
        fdname.close()
        fdcounty.close()
        
    if verbosity > 2:
        for i in range(N_districts):
            msg = '{0:{spaces}}'.format(new_dist_nm[i], spaces=16)
            msg += ' Num counties: {0:{spaces}}'.format(len(new_dist_county[i]), spaces=6)
            print(msg)
        
    return new_dist_nm, new_dist_county, new_dist_county_names, new_dist_geom 


def read_dist_data(fdir, county_nms, county_dist, county_dist_names, county_polys, verbosity = 0):
    
    N_counties = len(county_nms)
    N_districts = 25
    
    # list of polygons defining given district (this we don't read)
    new_dist_geom = [[] for i in range(N_districts)]
    
    # name of each district
    new_dist_nm = ['' for i in range(N_districts)]
    
    # id of counties in each district
    new_dist_county = [0 for i in range(N_districts)]
    
    # name of counties in each district
    new_dist_county_names = ['' for i in range(N_districts)]
    
    ## open files
    fdname = open(fdir + 'district_names.txt', 'r')
    fdcounty = open(fdir + 'district_counties.txt', 'r')
    
    lname = fdname.readline().split()
    lcounty = fdcounty.readline().split()
    cnt = 1
    
    while cnt <= N_districts:
        
       
        new_dist_nm[cnt-1] = get_one_string(lname)
        counties = []
        counties_nm = []
        for s in lcounty:
            counties.append(int(s))
            counties_nm.append(county_nms[int(s)])

        new_dist_county[cnt-1] = counties
        new_dist_county_names[cnt-1] = counties_nm
        
        lname = fdname.readline().split()
        lcounty = fdcounty.readline().split()
        cnt += 1
        
    if verbosity > 2:
        for i in range(N_districts):
            msg = '{0:{spaces}}'.format(new_dist_nm[i], spaces=16)
            msg += ' Num counties: {0:{spaces}}'.format(len(new_dist_county_nbr[i]), spaces=6)
            print(msg)
        
    return new_dist_nm, new_dist_county, new_dist_county_names, new_dist_geom 