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

class comeps:

    def __init__(self, args:argparse.Namespace) -> None:
        self.data = self.read_comeps(args)
        return


    def read_comeps(self, args:argparse.Namespace) -> xr.Dataset:

        base_file = "{}/mbr{:03d}/{:03d}".format(args.directory,0,0)
        if not ostools.does_file_exist(base_file): 
            print("Base file: {} was not found").format(base_file)
        
        lats, lons = grib2Read.get_latlons(self, base_file)
        
        no_members = args.members
        nt = args.files_per_member
        Nt_coords = np.zeros(nt, dtype=dt.datetime)
        
        precip = np.full([nt,no_members,lats.shape[0],lons.shape[1]], np.nan)
        rain = np.full([nt,no_members,lats.shape[0],lons.shape[1]], np.nan)
        cloudcover = np.full([nt,no_members,lats.shape[0],lons.shape[1]], np.nan)
        visibility = np.full([nt,no_members,lats.shape[0],lons.shape[1]], np.nan)

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
                            precip[k,m,:,:] = values.reshape(Nj, Ni)
                        else:
                            precip[k,m,:,:] = values.reshape(Nj, Ni) - precip[k-1,m,:,:]
                    
                    if iop==71 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                        values = ec.codes_get_values(gid)
                        cloudcover[k,m,:,:] = values.reshape(Nj, Ni)

                    if iop==20 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                        values = ec.codes_get_values(gid)
                        visibility[k,m,:,:] = values.reshape(Nj, Ni)
                    
                    if iop==181 and level==0 and typeOfLevel=='heightAboveGround' and levelType=='sfc':
                        values = ec.codes_get_values(gid)
                        if k == 0:
                            rain[k,m,:,:] = values.reshape(Nj, Ni)
                        else:
                            rain[k,m,:,:] = values.reshape(Nj, Ni) - rain[k-1,m,:,:]
                        
                ec.codes_release(gid)

        ds_grib = xr.Dataset(coords={"lat": (["x","y"], lats), 
                                     "lon": (["x","y"], lons),
                                     "member": (["m"], np.arange(no_members)),
                                     "time": (["t"], Nt_coords)})

        ds_grib['precip'] = (['time', 'member', 'lat', 'lon'], precip )
        ds_grib['tcc'] = (['time', 'member', 'lat', 'lon'], cloudcover )
        ds_grib['rain'] = (['time', 'member', 'lat', 'lon'], rain )
        ds_grib['visibility'] = (['time', 'member', 'lat', 'lon'], visibility )


        return ds_grib