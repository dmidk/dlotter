#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.read
Called from dlotter.__main__
"""
import argparse
import xarray as xr
import eccodes as ec
import netCDF4 as nc
import pygrib
import numpy as np
import datetime as dt
from dmit import regrot
from typing import Union

class grib2Read:
    """Class for reading grib files
    """

    def __init__(self, args:argparse.Namespace, files_to_read:list) -> None:
        """Constructor for grib2Read class

        Parameters
        ----------
        args : argparse.Namespace
            Namespace containing all the arguments

        files_to_read : list
            List of files to read

        """
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

        self.search_td2m = False
        self.found_td2m = False
        if 'td2m' in self.parameters:
            self.search_td2m = True

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

        self.search_slp = False
        self.found_slp  = False
        if 'slp' in self.parameters:
            self.search_slp = True

        self.search_tcc = False
        self.found_tcc  = False
        if 'tcc' in self.parameters:
            self.search_tcc = True

        self.search_lmhc = False
        self.found_lcc = False
        self.found_mcc = False
        self.found_hcc = False
        if 'lmhc' in self.parameters:
            self.search_lmhc = True

        self.search_snow = False
        self.found_snow = False
        if 'snow' in self.parameters:
            self.search_snow = True

        self.search_ws = False
        self.found_ws = False
        if 'ws' in self.parameters:
            self.search_ws = True

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
        if self.search_td2m: td2m = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_uv:
            u10 = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
            v10 = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_precip: precip = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
        if self.search_slp: slp = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
        if self.search_tcc: tcc = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_lmhc:
            lcc = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
            mcc = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
            hcc = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_snow: snow = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        if self.search_ws: ws = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)


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

            for _, gid in enumerate(gids):
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

                if self.search_td2m and (shortName=='td' or shortName=='2td') and level==2 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_td2m = True
                    td2m[k,:,:] = values.reshape(Nj, Ni)

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

                if self.search_precip and (shortName=='tp' or shortName=='tprate') \
                                        and level==0 and typeOfLevel=='heightAboveGround' \
                                        and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_precip = True
                    precip[k,:,:] = values.reshape(Nj, Ni)

                if self.search_slp and (shortName=='pres') and level==0 and \
                                        typeOfLevel=='heightAboveSea' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_slp = True
                    slp[k,:,:] = values.reshape(Nj, Ni)

                if self.search_lmhc and (shortName=='lcc') and level==0 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_lcc = True
                    lcc[k,:,:] = values.reshape(Nj, Ni)

                if self.search_lmhc and (shortName=='mcc') and level==0 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_mcc = True
                    mcc[k,:,:] = values.reshape(Nj, Ni)

                if self.search_lmhc and (shortName=='hcc') and level==0 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_hcc = True
                    hcc[k,:,:] = values.reshape(Nj, Ni)

                if self.search_tcc and (shortName=='tcc') and level==0 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_tcc = True
                    tcc[k,:,:] = values.reshape(Nj, Ni)

                if self.search_snow and (shortName=='tpsolid') and level==0 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_snow = True
                    snow[k,:,:] = values.reshape(Nj, Ni)

                if self.search_ws and (shortName=='ws') and level==10 and \
                                        typeOfLevel=='heightAboveGround' and levelType=='sfc':
                    values = ec.codes_get_values(gid)
                    self.found_ws = True
                    ws[k,:,:] = values.reshape(Nj, Ni)

                ec.codes_release(gid)

        ds_grib = xr.Dataset(coords={"lat": (["x","y"], lats),
                                     "lon": (["x","y"], lons),
                                     "time": (["t"], Nt_coords)})

        if self.found_t2m: ds_grib['t2m'] = (['time', 'lat', 'lon'], t2m - 273.15 )
        if self.found_td2m: ds_grib['td2m'] = (['time', 'lat', 'lon'], td2m - 273.15 )
        if self.found_u: ds_grib['u10m'] = (['time', 'lat', 'lon'], u10 )
        if self.found_u: ds_grib['v10m'] = (['time', 'lat', 'lon'], v10 )
        if self.found_precip: ds_grib['precip'] = (['time', 'lat', 'lon'], precip )
        if self.found_slp: ds_grib['slp'] = (['time', 'lat', 'lon'], slp * 0.01)
        if self.found_tcc: ds_grib['tcc'] = (['time', 'lat', 'lon'], tcc )
        if self.found_lcc: ds_grib['lcc'] = (['time', 'lat', 'lon'], lcc )
        if self.found_mcc: ds_grib['mcc'] = (['time', 'lat', 'lon'], mcc )
        if self.found_hcc: ds_grib['hcc'] = (['time', 'lat', 'lon'], hcc )
        if self.found_snow: ds_grib['snow'] = (['time', 'lat', 'lon'], snow )
        if self.found_ws: ds_grib['ws'] = (['time', 'lat', 'lon'], ws )

        if len(list(ds_grib.data_vars)) == 0:
            raise SystemExit('No variables found. This can be due to missing tables \
                              in ECCODES_DEFINITION_PATH or that the requested keys \
                              are not yet implemented')

        ds_grib = self.sort_by_time(ds_grib)

        return ds_grib


    def sort_by_time(self, dataarray:xr.Dataset) -> xr.Dataset:
        """Sort dataarray by time.

        Parameters
        ----------
        dataarray : xr.Dataset
            Dataarray to sort.

        Returns
        -------
        xr.Dataset
            Sortet dataarray.
        """

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
        """Get latitudes and longitudes from file. Uses pygrib as
        eccodes have no easy interface for that.

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
            # dy = g.jDirectionIncrementInDegrees
            # dx = g.iDirectionIncrementInDegrees
            latPole = g.latitudeOfSouthernPoleInDegrees
            lonPole = g.longitudeOfSouthernPoleInDegrees

            lons, lats = np.meshgrid(np.linspace(
                lonFirst, lonLast, londim), np.linspace(latFirst, latLast, latdim))

            lons, lats = regrot.rot_to_reg(lonPole, latPole, lons, lats)
        else:
            lats, lons = gr[1].latlons()

        gr.close()

        return lats, lons


    def get_gids(self, gribfile:str, TextIOWrapper:bool=False) -> list:
        """Get GribIDs (gid) for all the messages in one gribfile

        Parameters
        ----------
        gribfile : str
            path to gribfile
        TextIOWrapper : bool
            if file is open send the object instead and set TextIOWrapper to True

        Returns
        -------
        list
            list of grib-ids
        """

        if not TextIOWrapper:
            f = open(gribfile, 'rb')
            msg_count = ec.codes_count_in_file(f)
            gids = [ec.codes_grib_new_from_file(f) for i in range(msg_count)]
            f.close()
        else:
            #gribfile has already been opened into a TextIOWrapper in this case
            msg_count = ec.codes_count_in_file(gribfile)
            gids = np.zeros(msg_count, dtype=int)
            for i in range(msg_count):
                gids[i] = ec.codes_grib_new_from_file(gribfile)

        return gids



class netcdf2read:
    """Class for reading netcdf files.
    """

    def __init__(self, args:argparse.Namespace, files_to_read:list) -> None:
        """Constructor for netcdf2read class.

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        files_to_read : list
            List of files to read
        """
        if args.verbose:
            print("Reading NetCDF", flush=True)

        self.parameters = args.parameters

        self.set_bools()

        self.data = self.read(args, files_to_read)

        return


    def set_bools(self) -> None:
        """Set booleans for which variables to read.
        """

        self.search_t2m = False
        self.found_t2m = False
        if 't2m' in self.parameters:
            self.search_t2m = True

        self.search_precip = False
        self.found_precip = False
        if 'precip' in self.parameters:
            self.search_precip = True

        return


    def sort_by_time(self, dataarray:xr.Dataset) -> xr.Dataset:
        """Sort xarray dataset by time.

        Parameters
        ----------
        dataarray : xr.Dataset
            Dataset to sort

        Returns
        -------
        xr.Dataset
            Sorted dataset
        """

        nt = dataarray.dims['time']
        parameters = list(dataarray.data_vars)

        time = dataarray['time'].values
        idx = np.argsort(time)

        da = dataarray.sortby('time')

        for p in parameters:
            for k in range(nt):
                da[p][k,:,:] = dataarray[p][idx[k],:,:]

        return da


    def read(self, args:argparse.Namespace, files_to_read:list) -> xr.Dataset:
        """Fetch data from the netcdf files and return an xarray.Dataset()

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments to dlotter.__main__
        files_to_read : list
            List of str with paths to netcdf files

        Returns
        -------
        xr.Dataset
            Dataset with found parameters
        """

        lats, lons = self.get_latlons(files_to_read[0])

        Nt = len(files_to_read)
        Nt_coords = np.zeros(Nt, dtype=dt.datetime)

        if self.search_t2m: t2m = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)
        if self.search_precip: precip = np.full([Nt,lats.shape[0],lons.shape[1]], np.nan)

        for k,f in enumerate(files_to_read):
            print(k,f)

            f = nc.Dataset(files_to_read[k])

            forecast = f.getncattr('ValidDate')
            forecast = dt.datetime.strptime(forecast, '%Y-%b-%d %H:%M:%S')

            Nt_coords[k] = forecast

            if self.search_t2m:
                t2m_key = self.find_relevant_key(f, 't2m')
                if t2m_key is not None:
                    t2m[k,:,:] = f[t2m_key][:,:]
                    self.found_t2m = True

            if self.search_precip:
                precip_key = self.find_relevant_key(f, 'precip')
                if precip_key is not None:
                    precip[k,:,:] = f[precip_key][:,:]
                    self.found_precip = True


            f.close()


        ds_grib = xr.Dataset(coords={"lat": (["x","y"], lats),
                                     "lon": (["x","y"], lons),
                                     "time": (["t"], Nt_coords)})

        if self.found_t2m: ds_grib['t2m'] = (['time', 'lat', 'lon'], t2m - 273.15 )
        if self.found_precip: ds_grib['precip'] = (['time', 'lat', 'lon'], precip)

        if len(list(ds_grib.data_vars)) == 0:
            raise SystemExit('No variables found. This can be due to missing tables in \
                              ECCODES_DEFINITION_PATH or that the requested keys are \
                              not yet implemented')

        ds_grib = self.sort_by_time(ds_grib)

        return ds_grib


    def get_latlons(self, netcdf_file:str) -> tuple:
        """Get latitudes and longitudes from file. Uses netcdf4-python as
        eccodes have no easy interface for that.

        Parameters
        ----------
        netcdf_file : str
            Path to netcdf file

        Returns
        -------
        tuple
            tuple of latitudes, longitudes
        """

        f = nc.Dataset(netcdf_file)

        latkey = self.find_relevant_key(f, 'latitude')
        lonkey = self.find_relevant_key(f, 'longitude')

        if latkey is None or lonkey is None:
            raise SystemExit('No lat/lon found in netcdf file')

        lats = f[latkey][:,:]
        lons = f[lonkey][:,:]

        f.close()

        return lats, lons


    def find_relevant_key(self, dataset:nc._netCDF4.Dataset, key:str) -> Union[str, None]:
        """Check if a key is in the netcdf file.

        Parameters
        ----------
        nc._netCDF4.Dataset
            Netcdf file
        key : str
            Key to search for

        Returns
        -------
        str
            None if no key is found, otherwise returns the best guess key
        """

        available_keys = dataset.variables.keys()

        return_key = None

        if key == 'latitude':
            if 'lat' in available_keys: return_key = 'lat'
            elif 'latitude' in available_keys: return_key = 'latitude'
            elif 'LAT' in available_keys: return_key = 'LAT'

        elif key == 'longitude':
            if 'lon' in available_keys: return_key = 'lon'
            elif 'longitude' in available_keys: return_key = 'longitude'
            elif 'LON' in available_keys: return_key = 'LON'

        elif key == 't2m':
            if 't2m' in available_keys: return_key = 't2m'
            elif 't2maboveground' in available_keys: return_key = 't2maboveground'
            elif 'T2M' in available_keys: return_key = 'T2M'

        elif key == 'precip':
            if 'precip' in available_keys: return_key = 'precip'
            elif 'precipitation' in available_keys: return_key = 'precipitation'
            elif 'PRECIP' in available_keys: return_key = 'PRECIP'

        return return_key
