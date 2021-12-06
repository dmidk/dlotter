#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.read
Called from dlotter.__main__
"""
import sys
import argparse
import xarray as xr
import eccodes as ec
import pygrib
import numpy as np
import datetime as dt
from dmit import regrot

class grib2Read:

    def __init__(self, args:argparse.Namespace, files_to_read:list) -> None:
        if args.verbose:
            print("Reading GRIB2", flush=True)

        self.parameters = args.parameters

        self.set_bools()

        self.data = self.read(args, files_to_read)

        return


    def set_bools(self) -> None:
        """Set bools for use in the module
        """

        self.search_t2m = False
        self.found_t2m = False

        if 't2m' in self.parameters: 
            self.search_t2m = True
            
        self.search_uv=False    
        self.found_u = False
        self.found_v = False
        self.found_uv=False
        if 'w10m' in self.parameters: 
            self.search_uv = True

        self.search_precip = False
        self.found_precip = False
        if 'precip' in self.parameters:
            self.search_precip = True
            

        return


    def read(self, args:argparse.Namespace, files_to_read:list) -> xr.Dataset:
        """Fetch data from the gribfile and return an xarray.Dataset()

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments to dlotter.__main__
        files_to_read : list
            List of str with paths to gribfiles

        Returns
        -------
        xr.Dataset
            Dataset with found parameters
        """

        grid_idx = 0 # Set to 0
        lats, lons = self.get_latlons(files_to_read[grid_idx])

        Nt = len(files_to_read)
        Nt_coords = np.zeros(Nt, dtype=dt.datetime)
        
        if self.search_t2m: t2m = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_uv: 
            u10 = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
            v10 = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_precip: precip = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        for k,f in enumerate(files_to_read):
            gids = self.get_gids(f)

            time_gid = gids[0]

            ec.codes_set(time_gid, 'stepUnits', 'm')

            date = ec.codes_get(time_gid, 'dataDate')
            time = ec.codes_get(time_gid, 'dataTime')
            lead = ec.codes_get(time_gid, 'step')

            analysis = dt.datetime.strptime("{:d}-{:04d}".format(date,time), '%Y%m%d-%H%M')
            forecast = analysis + dt.timedelta(minutes=lead)
            Nt_coords[k] = forecast

            for i, gid in enumerate(gids):
                shortName = ec.codes_get(gid, 'shortName')
                level = ec.codes_get(gid, 'level')
                typeOfLevel = ec.codes_get(gid, 'typeOfLevel')
                levelType = ec.codes_get(gid, 'levelType')

                Ni = ec.codes_get(gid, 'Ni')
                Nj = ec.codes_get(gid, 'Nj')

                if levelType=='103': levelType='sfc' # For grib2, leveltype 103 is surface

                if self.search_t2m and (shortName=='t' or shortName=='2t') and level==2 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_t2m = True
                    t2m[k,:,:] = values.reshape(Nj, Ni)

                if self.search_uv and (shortName=='u' or shortName=='10u') and level==10 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_u = True
                    u10[k,:,:] = values.reshape(Nj, Ni)

                if self.search_uv and (shortName=='v' or shortName=='10v') and level==10 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_v = True
                    v10[k,:,:] = values.reshape(Nj, Ni)

                if self.search_precip and (shortName=='tp' or shortName=='tprate') and level==0 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_precip = True
                    precip[k,:,:] = values.reshape(Nj, Ni)

                ec.codes_release(gid)

        ds_grib = xr.Dataset(coords={"lat": (["x","y"], lats), 
                                     "lon": (["x","y"], lons), 
                                     "time": (["t"], Nt_coords)})

        if self.found_t2m: ds_grib['t2m'] = (['time', 'lat', 'lon'], t2m - 273.15 )
        if self.found_u: ds_grib['u10m'] = (['time', 'lat', 'lon'], u10 )
        if self.found_u: ds_grib['v10m'] = (['time', 'lat', 'lon'], v10 )
        if self.found_precip: ds_grib['precip'] = (['time', 'lat', 'lon'], precip )

        if len(list(ds_grib.data_vars)) == 0:
            raise SystemExit('No variables found. This can be due to missing tables in ECCODES_DEFINITION_PATH or that the requested keys are not yet implemented')
 
        ds_grib = self.sort_by_time(ds_grib)

        return ds_grib

    
    def sort_by_time(self, dataarray:xr.Dataset) -> xr.Dataset:

        nt = dataarray.dims['time']
        parameters = list(dataarray.data_vars)

        time = dataarray['time'].values
        idx = np.argsort(time)

        da = dataarray.sortby('time')

        for p in parameters:
            for k in range(nt):
                da[p][k,:,:] = dataarray[p][idx[k],:,:]
        
        return da


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
        g = gr[1]

        if g['gridType'] == 'rotated_ll':
            print('Found rotated grid, discarding pygrib extraction')
            latdim = g.Nj
            londim = g.Ni

            latFirst = g.latitudeOfFirstGridPointInDegrees
            latLast  = g.latitudeOfLastGridPointInDegrees
            lonFirst = (g.longitudeOfFirstGridPointInDegrees % 180.)-180.
            lonLast  = g.longitudeOfLastGridPointInDegrees
            dy = g.jDirectionIncrementInDegrees
            dx = g.iDirectionIncrementInDegrees
            latPole = g.latitudeOfSouthernPoleInDegrees
            lonPole = g.longitudeOfSouthernPoleInDegrees

            lons, lats = np.meshgrid(np.linspace(
                lonFirst, lonLast, londim), np.linspace(latFirst, latLast, latdim))

            lons, lats = regrot.rot_to_reg(lonPole, latPole, lons, lats)
        else:
            lats, lons = gr[1].latlons()

        gr.close()

        return lats, lons

    
    def get_gids(self, gribfile:str) -> list:
        """Get GribIDs (gid) for all the messages in one gribfile

        Parameters
        ----------
        gribfile : str
            path to gribfile

        Returns
        -------
        list
            list of grib-ids
        """
        
        f = open(gribfile, 'rb')
        msg_count = ec.codes_count_in_file(f)
        gids = [ec.codes_grib_new_from_file(f) for i in range(msg_count)]
        f.close()

        return gids




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
