#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.comeps
Called from dlotter.__main__
"""
import sys
import argparse
import xarray as xr
import eccodes as ec
import pygrib
import numpy as np
import datetime as dt
from dmit import ostools
from .read import grib2Read

from scipy.spatial import cKDTree

class comeps:

    def __init__(self, args:argparse.Namespace) -> None:
        self.data = self.read_comeps(args)
        return


    def get_location_index(self, location, tree) -> int:
        location = location.split(',')

        ref_lat = float(location[0])
        ref_lon = float(location[1])

        dist, idx = tree.query((ref_lat,ref_lon))

        return idx


    def read_comeps(self, args:argparse.Namespace) -> xr.Dataset:

        base_file = "{}/mbr{:03d}/{:03d}".format(args.directory,0,0)
        if not ostools.does_file_exist(base_file): 
            print("Base file: {} was not found").format(base_file)
        
        lats, lons = grib2Read.get_latlons(self, base_file)
        
        no_members = args.members
        nt = args.files_per_member
        Nt_coords = np.zeros(nt, dtype=dt.datetime)
        
        locations = args.latlon.split(':')
        no_locations = len(locations)

        data_locations = list(zip(lats.flatten(), lons.flatten()))
        tree = cKDTree(data_locations)

        data_indexes = np.empty(no_locations, dtype=int)

        i=0
        for loc in locations:
            location_idx = self.get_location_index(loc, tree)
            data_indexes[i] = location_idx
            i+=1
        
        precip = np.full([nt, no_members, no_locations], np.nan)
        rain = np.full([nt, no_members, no_locations], np.nan)
        cloudcover = np.full([nt, no_members, no_locations], np.nan)
        visibility = np.full([nt, no_members, no_locations], np.nan)

        for l in range(no_locations):
            loc_idx = data_indexes[l]

            for k in range(nt):
                for m in range(no_members):

                    epsfile = "{}/mbr{:03d}/{:03d}".format(args.directory,m,k)
                    if not ostools.does_file_exist(epsfile): continue

                    gids = grib2Read.get_gids(self, epsfile)

                    time_gid = gids[0]

                    ec.codes_set(time_gid, 'stepUnits', 'm')

                    date = ec.codes_get(time_gid, 'dataDate')
                    time = ec.codes_get(time_gid, 'dataTime')
                    lead = ec.codes_get(time_gid, 'step')

                    analysis = dt.datetime.strptime(('%i-%.2i')%(date,time),'%Y%m%d-%H%M')
                    forecast = analysis + dt.timedelta(minutes=lead)

                    Nt_coords[k] = forecast

                    for i, gid in enumerate(gids):
                        shortName = ec.codes_get(gid, 'shortName')
                        level = ec.codes_get(gid, 'level')
                        typeOfLevel = ec.codes_get(gid, 'typeOfLevel')
                        levelType = ec.codes_get(gid, 'levelType')
                        iop = ec.codes_get(gid, 'indicatorOfParameter')

                        Ni = ec.codes_get(gid, 'Ni')
                        Nj = ec.codes_get(gid, 'Nj')

                        if iop==61 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                            values = ec.codes_get_values(gid)
                            if k == 0:
                                precip[k,m,l] = values[loc_idx]
                            else:
                                precip[k,m,l] = values[loc_idx] - precip[k-1,m,l]
                        
                        if iop==71 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                            values = ec.codes_get_values(gid)
                            cloudcover[k,m,l] = values[loc_idx]

                        if iop==20 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                            values = ec.codes_get_values(gid)
                            visibility[k,m,l] = values[loc_idx]
                        
                        if iop==181 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                            values = ec.codes_get_values(gid)
                            if k == 0:
                                rain[k,m,l] = values[loc_idx]
                            else:
                                rain[k,m,l] = values[loc_idx] - rain[k-1,m,l]
                            
                    ec.codes_release(gid)


        ds_grib = xr.Dataset(coords={"location": (["l"], np.arange(no_locations)), 
                                     "member": (["m"], np.arange(no_members)),
                                     "time": (["t"], Nt_coords)})

        ds_grib['precip'] = (['time', 'member', 'locations'], precip )
        ds_grib['tcc'] = (['time', 'member', 'locations'], cloudcover )
        ds_grib['rain'] = (['time', 'member', 'locations'], rain )
        ds_grib['visibility'] = (['time', 'member', 'locations'], visibility )


        return ds_grib