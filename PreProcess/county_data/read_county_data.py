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


def read_county_map(fdir = 'data/', verbosity = 0):
    
    ## read shapefile
    tx_county_shp = gpd.read_file(fdir + 'shp/County.shp')
    
    ## extract key data
    county_nbr = tx_county_shp['CNTY_NBR']
    county_nm = tx_county_shp['CNTY_NM']
    county_dist = tx_county_shp['DIST_NBR']
    county_dist_names = tx_county_shp['DIST_NM']
    county_geom = tx_county_shp['geometry']
    
    ## data has scattered county data which 
    ## we gather in vector with number of counties
    N_polys = len(county_nm)
    N_counties = 254
    
    # list of polygons defining given county
    new_county_polys = [[] for i in range(N_counties)]
    
    # name of each county
    new_county_nms = ['' for i in range(N_counties)]
    
    # id of district the county belongs to
    new_county_dist = [0 for i in range(N_counties)]
    
    # name of district the county belongs to
    new_county_dist_names = ['' for i in range(N_counties)]
    
    # populate data
    for i in range(len(county_nbr)):
        # find its county id
        cid = county_nbr[i]-1

        # polygon goes to cid^th place
        new_county_polys[cid].append(county_geom[i])

        # modify the name of cid county
        new_county_nms[cid] = county_nm[i]

        # get district this county belongs to
        new_county_dist[cid] = county_dist[i] - 1
        new_county_dist_names[cid] = county_dist_names[i]

    if verbosity > 1:
        print('size of new_county_polys: {}'.format(len(new_county_polys)))
        print('size of new_county_nms: {}'.format(len(new_county_nms)))
    
    ## write to file
    write_to_file = False
    if write_to_file:
        # save county data to file for later use
        #fcpoly = open('county_data/county_polygons.txt', 'w')
        fcname = open(fdir + 'county_names.txt', 'w')
        fcdist = open(fdir + 'county_districts.txt', 'w')
        fcdistnm = open(fdir + 'county_district_names.txt', 'w')
        for i in range(N_counties):
            #fcpoly.write('{}\n'.format(new_county_polys[i]))
            fcname.write("{}\n".format(new_county_nms[i]))
            fcdist.write('{}\n'.format(new_county_dist[i]))
            fcdistnm.write("{}\n".format(new_county_dist_names[i]))

        #fcpoly.close()
        fcname.close()
        fcdist.close()
        fcdistnm.close()
        
    if verbosity > 2:
        for i in range(N_counties):
            msg = '{0:{spaces}}'.format(new_county_nms[i], spaces=16)
            msg += ' {0:{spaces}}'.format(new_county_dist[i], spaces=6)
            msg += ' {0:{spaces}}'.format(new_county_dist_names[i], spaces=16)
            print(msg)
        
    return new_county_nms, new_county_dist, new_county_dist_names, new_county_polys 


def read_county_data(fdir = "data/", verbosity = 0):
    
    N_counties = 254
    
    # list of polygons defining given county (this we don't read)
    new_county_polys = [[] for i in range(N_counties)]
    
    # name of each county
    new_county_nms = ['' for i in range(N_counties)]
    
    # id of district the county belongs to
    new_county_dist = [0 for i in range(N_counties)]
    
    # name of district the county belongs to
    new_county_dist_names = ['' for i in range(N_counties)]
    
    ## open files
    fcname = open(fdir + 'county_names.txt', 'r')
    fcdist = open(fdir + 'county_districts.txt', 'r')
    fcdistnm = open(fdir + 'county_district_names.txt', 'r')
    
    lname = fcname.readline().split()
    ldist = fcdist.readline().split()
    ldistnm = fcdistnm.readline().split()
    cnt = 1
    
    while cnt <= N_counties:
        
       
        new_county_nms[cnt-1] = get_one_string(lname)
        new_county_dist[cnt-1] = int(ldist[0])
        new_county_dist_names[cnt-1] = get_one_string(ldistnm)
        
        lname = fcname.readline().split()
        ldist = fcdist.readline().split()
        ldistnm = fcdistnm.readline().split()
        cnt += 1
        
    if verbosity > 2:
        for i in range(N_counties):
            msg = '{0:{spaces}}'.format(new_county_nms[i], spaces=16)
            msg += ' {0:{spaces}}'.format(new_county_dist[i], spaces=6)
            msg += ' {0:{spaces}}'.format(new_county_dist_names[i], spaces=16)
            print(msg)
        
    return new_county_nms, new_county_dist, new_county_dist_names, new_county_polys 