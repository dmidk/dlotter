#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.read
Called from dlotter.__main__
"""
import sys
import argparse
import pandas as pd
import eccodes as ec
import pygrib
import numpy as np

class grib2Read:

    def __init__(self, args:argparse.Namespace, files_to_read:list) -> None:
        if args.verbose:
            print("Reading GRIB2", flush=True)

        data = self.read(args, files_to_read)
        
        return


    def read(self, args:argparse.Namespace, files_to_read:list) -> pd.DataFrame:
        grid_idx = 0 # Set to 0
        lats, lons = self.get_latlons(files_to_read[grid_idx])

        return


    def get_latlons(self, gribfile:str) -> tuple:
        """Get latitudes and longitudes from file. Uses pygrib as eccodes have no easy interface for that.

        Parameters
        ----------
        gribfile : str
            Path to gribfile

        Returns
        -------
        tuple
            tuple of latitudes, longitudes
        """

        gr = pygrib.open(gribfile)
        lats, lons = gr[1].latlons()
        gr.close()

        return lats, lons


#  print(gribfile)
#         f = open(gribfile, 'rb')
#         msg_count = ec.codes_count_in_file(f)
#         gid_list = [ec.codes_grib_new_from_file(f) for i in range(msg_count)]
#         f.close()

    # gid = gid_list[0]


    # def get_grid(self, gid: int) -> None:


    #     Ni = ec.codes_get(gid, 'Ni')
    #     Nj = ec.codes_get(gid, 'Nj')

    #     gridtype = ec.codes_get(gid, 'gridType')

    #     test = ec.codes_grib_get_data(gid)
    #     print(test)
    #     # if gridtype == 'lambert':
    #     #     projparams['proj']='lcc'
    #     #     projparams['lon_0']=self['LoVInDegrees']
    #     #     projparams['lat_0']=self['LaDInDegrees']
    #     #     projparams['lat_1']=self['Latin1InDegrees']
    #     #     projparams['lat_2']=self['Latin2InDegrees']

    #     lat_first = ec.codes_get(gid, 'latitudeOfFirstGridPointInDegrees')
    #     lon_first = ec.codes_get(gid, 'longitudeOfFirstGridPointInDegrees')
    #     # lat_last  = ec.codes_get(gid, 'latitudeOfLastGridPointInDegrees')
    #     # lon_last  = ec.codes_get(gid, 'longitudeOfLastGridPointInDegrees')
        
    #     return 