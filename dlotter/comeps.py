#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.comeps
Called from dlotter.__main__
"""
import argparse
import xarray as xr
import eccodes as ec
import numpy as np
import datetime as dt
from dmit import ostools
import dmit.calc.sun as sun
from .read import grib2Read
import gc

from scipy.spatial import cKDTree

class comeps:
    """Class for reading COMEPS data
    """

    def __init__(self, args:argparse.Namespace) -> None:
        """Constructor for comeps class

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from argparse
        """
        self.data = self.read_comeps(args)
        return


    def get_location_index(self, location:str, tree:cKDTree) -> int:
        """_summary_

        Parameters
        ----------
        location : str
            String of location with latitude and longitude
        tree : cKDTree
            scipy.spatial.cKDTree object

        Returns
        -------
        int
            Index of location in tree
        """
        location = location.split(',')

        ref_lat = float(location[0])
        ref_lon = float(location[1])

        _, idx = tree.query((ref_lat,ref_lon))

        return idx


    def read_comeps(self, args:argparse.Namespace) -> xr.Dataset:
        """Read COMEPS data

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from argparse

        Returns
        -------
        xr.Dataset
            Dataset with COMEPS data
        """

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
        data_lats = np.empty(no_locations, dtype=np.single)
        data_lons = np.empty(no_locations, dtype=np.single)

        i=0
        for loc in locations:
            location_idx = self.get_location_index(loc, tree)
            data_indexes[i] = location_idx
            data_lats[i] = lats.flatten()[location_idx]
            data_lons[i] = lons.flatten()[location_idx]

            i+=1

        del tree

        precip = np.zeros([nt, no_members, no_locations], dtype=np.single) + np.nan
        precip_solid = np.zeros([nt, no_members, no_locations], dtype=np.single) + np.nan
        cloudcover = np.zeros([nt, no_members, no_locations], dtype=np.single) + np.nan
        visibility = np.zeros([nt, no_members, no_locations], dtype=np.single) + np.nan
        # night 0 if day, 1 if night (note not dependent on number of members)
        night = np.zeros([nt, no_locations], dtype=np.single) + np.nan

        for l in range(no_locations):
            loc_idx = data_indexes[l]

            for k in range(nt):

                for m in range(no_members):
                    if args.verbose: print("\nTime: {}/{}, Member: {}/{}".format(k,nt,m,no_members))
                    epsfile = "{}/mbr{:03d}/{:03d}".format(args.directory,m,k)
                    if not ostools.does_file_exist(epsfile): continue

                    f = open(epsfile, 'r')
                    if args.verbose: print("Reading: {}".format(epsfile))

                    gids = grib2Read.get_gids(self, f, TextIOWrapper=True)

                    time_gid = gids[0]

                    ec.codes_set(time_gid, 'stepUnits', 'm')

                    date = ec.codes_get(time_gid, 'dataDate')
                    time = ec.codes_get(time_gid, 'dataTime')
                    lead = ec.codes_get(time_gid, 'step')

                    analysis = dt.datetime.strptime(('%i-%.2i')%(date,time),'%Y%m%d-%H%M')
                    forecast = analysis + dt.timedelta(minutes=lead)

                    if k==0 and m==0:
                        # We'll save some time and only calculate one sunset and sunrise
                        sunrise, sunset = self.get_sunrise_sunset(data_lats[l], data_lons[l],
                                                                  forecast)
                    if m==0:
                        night[k,l] = self.is_it_night(sunrise, sunset, forecast)

                    Nt_coords[k] = forecast

                    for i, gid in enumerate(gids):
                        #shortName = ec.codes_get(gid, 'shortName')
                        level = ec.codes_get(gid, 'level')
                        typeOfLevel = ec.codes_get(gid, 'typeOfLevel')
                        levelType = ec.codes_get(gid, 'levelType')
                        iop = ec.codes_get(gid, 'indicatorOfParameter')

                        #Ni = ec.codes_get(gid, 'Ni')
                        #Nj = ec.codes_get(gid, 'Nj')

                        if iop==61 and level==0 and typeOfLevel=='heightAboveGround' \
                            and levelType=='sfc':

                            values = ec.codes_get_values(gid)
                            if k == 0:
                                precip[k,m,l] = values[loc_idx]
                            else:
                                precip[k,m,l] = values[loc_idx] - precip[k-1,m,l]

                        if iop==71 and level==0 and typeOfLevel=='heightAboveGround' \
                            and levelType=='sfc':

                            values = ec.codes_get_values(gid)
                            cloudcover[k,m,l] = values[loc_idx]

                        if iop==20 and level==0 and typeOfLevel=='heightAboveGround' \
                            and levelType=='sfc':

                            values = ec.codes_get_values(gid)
                            visibility[k,m,l] = values[loc_idx]

                        if iop==185 and level==0 and typeOfLevel=='heightAboveGround' \
                            and levelType=='sfc':

                            values = ec.codes_get_values(gid)
                            if k == 0:
                                precip_solid[k,m,l] = values[loc_idx]
                            else:
                                precip_solid[k,m,l] = values[loc_idx] - precip_solid[k-1,m,l]

                        ec.codes_release(gid)

                    f.close()

                # Member loop closed
                gc.collect()
            # Time loop closed
        # Location loop closed


        ds_grib = xr.Dataset(coords={"location": (["l"], np.arange(no_locations)),
                                     "member": (["m"], np.arange(no_members)),
                                     "time": (["t"], Nt_coords)})

        ds_grib['precip'] = (['time', 'member', 'locations'], precip)
        ds_grib['tcc'] = (['time', 'member', 'locations'], cloudcover)
        ds_grib['precip_solid'] = (['time', 'member', 'locations'], precip_solid)
        ds_grib['visibility'] = (['time', 'member', 'locations'], visibility)
        ds_grib['night'] = (['time','locations'], night)


        return ds_grib


    def get_sunrise_sunset(self, latitude:float, longitude:float, time:dt.datetime) -> tuple:
        """Get time of sunrise and sunset for a given location and time.

        Parameters
        ----------
        latitude : float
            Latitude of location.
        longitude : float
            Longitude of location.
        time : dt.datetime
            Time (day) of interest.

        Returns
        -------
        tuple
            Sunrise and sunset time.
        """
        sunrise = sun.calc_sun_time(longitude, latitude, time, setrise='sunrise')
        sunset = sun.calc_sun_time(longitude, latitude, time, setrise='sunset')

        return sunrise, sunset


    def is_it_night(self, sunrise:dt.datetime, sunset:dt.datetime, time:dt.datetime) -> int:
        """Get information if it is night or day.

        Parameters
        ----------
        sunrise : dt.datetime
            Time of sunrise
        sunset : dt.datetime
            Time of sunset
        time : dt.datetime
            current time

        Returns
        -------
        int
            0 if day, 1 if night
        """

        if time.hour >= sunrise.hour and time.hour < sunset.hour:
            binary_day = 0 # Day
        else:
            binary_day = 1 # Night

        return binary_day