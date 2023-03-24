#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.plotdiff
Called from dlotter.__main__
"""
import sys
import argparse
import xarray as xr
import numpy as np
import datetime as dt

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

import cartopy.crs as ccrs
import cartopy.feature as cfeature


class plotdiff:
    """Class for creating geographic plots of differences
    """

    def __init__(self, args:argparse.Namespace, data_list:xr.Dataset) -> None:
        """Constructor for plotdiff class

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        data : list of xr.Dataset
            Data to plot
        """

        self.projections(args)

        parameters = args.parameters
        avail_parameters = list(data_list[0].data_vars)

        self.nt = data_list[0].dims['time']

        print('- Plotting differences', flush=True)

        if 't2m' in parameters:
            if 't2m' not in avail_parameters:
                print('t2m was not found in available parameters: {},\
                       cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plotdiff_t2m(args, data_list)

        if 'td2m' in parameters:
            if 'td2m' not in avail_parameters:
                print('td2m was not found in available parameters: {},\
                       cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plotdiff_td2m(args, data_list)

        if 'w10m' in parameters:
            if 'u10m' not in avail_parameters or 'v10m' not in avail_parameters:
                print('u10m/v10m was not found in available parameters: {},\
                       cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plotdiff_w10m(args, data_list)

        if 'precip' in parameters:
            if 'precip' not in avail_parameters:
                print('precip was not found in available parameters: {},\
                       cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plotdiff_precip(args, data_list)

        if 'z' in parameters:
            if 'z' not in avail_parameters:
                print('z was not found in available parameters: {},\
                       cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plotdiff_z(args, data_list)

        return


    def projections(self,args:argparse.Namespace) -> None:
        """Set projection for plotting

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        """
        if args.area == 'dk':
            self.projection = ccrs.AlbersEqualArea(central_longitude=11.0, central_latitude=0.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [6, 16, 53, 59]

        if args.area == 'gl':
            self.projection = ccrs.TransverseMercator(central_longitude=-50.0,
                                    central_latitude=75.0, false_easting=0.0,
                                    false_northing=0.0,)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-75, -15, 58, 88]

        if args.area == 'neu':
            self.projection = ccrs.AlbersEqualArea(central_longitude=11.0, central_latitude=0.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-2.1, 20, 51, 65]

        if args.area == 'disko':
            self.projection = ccrs.AlbersEqualArea(central_longitude=-51.7, central_latitude=70,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-56, -46.5, 66.33, 71.52]

        if args.area == 'sjalland':
            self.projection = ccrs.AlbersEqualArea(central_longitude=11.0, central_latitude=0.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [10, 13.4, 54.3, 56.45]

        if args.area == 'europe':
            self.projection = ccrs.AlbersEqualArea(central_longitude=11.0, central_latitude=0.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-20, 30, 40, 75]

        if args.area == 'faroes':
            self.projection = ccrs.AlbersEqualArea(central_longitude=353.0, central_latitude=62.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-8, -6, 61.3, 62.5]

        if args.area == 'scoresbysund':
            self.projection = ccrs.AlbersEqualArea(central_longitude=-25., central_latitude=71.5,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-30.5, -20.5, 69.7, 73.7]

        if args.area == 'sgl':
            self.projection = ccrs.AlbersEqualArea(central_longitude=-42.5, central_latitude=61,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-52, -37, 59, 63]

        if args.area == 'tas':
            self.projection = ccrs.AlbersEqualArea(central_longitude=-37.7, central_latitude=65.9,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-43.0, -32.4, 63.9, 67.9]

        if args.area == 'nkb':
            self.projection = ccrs.AlbersEqualArea(central_longitude=-52.0, central_latitude=63.2,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [-57.5, -46.5, 61.5, 68.5]

        return


    def plotdiff_t2m(self, args:argparse.Namespace, data_list:list) -> None:
        """Plot 2m temperature difference

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        data_list : list of xr.Dataset
            Data to plot
        """

        exps = args.experiments.split(",")

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data_list[0]['lon'], data_list[0]['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data_list[0])

        colors = levels_and_colors.t2m.colors
        levels = [k for k in levels_and_colors.t2m.levels]

        cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend='both')

        analysis = data_list[0]['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data_list[0]['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if ( (self.check_for_empty_array(data_list[0]['t2m'][k,:,:])) | 
                 (self.check_for_empty_array(data_list[1]['t2m'][k,:,:])) ): continue

            self.add_title(axes,valid_time,analysis,'2m Temperature Difference: {}-{}'.format(exps[0],exps[1]))

            cs = plt.pcolormesh(plons, plats, data_list[0]['t2m'][k,:,:]-data_list[1]['t2m'][k,:,:],
                                cmap=cmap,
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels, extend='both')
            cb.set_label(r"$^\circ C$", rotation=270)

            fig.canvas.draw()

            figure_name = "{}/T2Mdiff_{}_{}_{}-{}.png".format(args.output_dir,
                                                    exps[0], exps[1],
                                                    analysis.strftime('%Y%m%d_%H%M'),
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)

        return


    def plotdiff_td2m(self, args:argparse.Namespace, data_list:list) -> None:
        """Plot 2m dew point temperature difference

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        data_list : list of xr.Dataset
            Data to plot
        """
        exps = args.experiments.split(",")

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data_list[0]['lon'], data_list[0]['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data_list[0])

        colors = levels_and_colors.t2m.colors
        levels = [k for k in levels_and_colors.t2m.levels]

        cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend='both')

        analysis = data_list[0]['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data_list[0]['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if ( (self.check_for_empty_array(data_list[0]['td2m'][k,:,:])) | 
                 (self.check_for_empty_array(data_list[1]['td2m'][k,:,:])) ): continue

            self.add_title(axes,valid_time,analysis,'2m Dew Point Temperature Difference: {}-{}'.format(exps[0],exps[1]))

            cs = plt.pcolormesh(plons, plats, data_list[0]['td2m'][k,:,:]-data_list[1]['td2m'][k,:,:],
                                cmap=cmap,
                                norm=norm,
                                transform=self.data_crs)


            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels, extend='both')
            cb.set_label(r"$^\circ C$", rotation=270)

            fig.canvas.draw()

            figure_name = "{}/TD2Mdiff_{}_{}_{}-{}.png".format(args.output_dir,
                                                    exps[0], exps[1],
                                                    analysis.strftime('%Y%m%d_%H%M'),
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)

        return


    def plotdiff_w10m(self, args:argparse.Namespace, data_list:list) -> None:
        """Plot 10m wind speed

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        data_list : list of xr.Dataset
            Data to plot
        """
        exps = args.experiments.split(",")

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = np.array(data_list[0]['lon']), np.array(data_list[0]['lat'])
        plons, plats = self.get_pcolormesh_center_coordinates(data_list[0])

        colors = levels_and_colors.w10m.colors
        levels = [k for k in levels_and_colors.w10m.levels]

        cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend='both')

        analysis = data_list[0]['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data_list[0]['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if ( (self.check_for_empty_array(data_list[0]['u10m'][k,:,:])) | 
                 (self.check_for_empty_array(data_list[1]['u10m'][k,:,:])) ): continue
            if ( (self.check_for_empty_array(data_list[0]['v10m'][k,:,:])) | 
                 (self.check_for_empty_array(data_list[1]['v10m'][k,:,:])) ): continue

            self.add_title(axes,valid_time,analysis,'10m Wind Difference: {}-{}'.format(exps[0],exps[1]))

            w10m = []
            for i in range(len(data_list)):
                u = data_list[i]['u10m'][k,:,:].values
                v = data_list[i]['v10m'][k,:,:].values
    
                w10m.append(np.sqrt(u*u+v*v))

            cs = plt.pcolormesh(plons, plats, w10m[0]-w10m[1],
                                cmap=cmap,
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels, extend='both')
            cb.set_label(r"$m/s$", rotation=270)


            fig.canvas.draw()

            figure_name = "{}/W10Mdiff_{}_{}_{}-{}.png".format(args.output_dir,
                                                    exps[0], exps[1],
                                                    analysis.strftime('%Y%m%d_%H%M'),
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)

        return


    def plotdiff_precip(self, args:argparse.Namespace, data_list:list) -> None:
        """Plot precipitation

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        data_list : list of xr.Dataset
            Data to plot
        """
        exps = args.experiments.split(",")

        # Fix that pcolormesh uses cell lower left corners
        plons, plats = self.get_pcolormesh_center_coordinates(data_list[0])

        colors = levels_and_colors.precip.colors
        levels = [k for k in levels_and_colors.precip.levels]

        cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend='both')

        analysis = data_list[0]['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data_list[0]['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if ( (self.check_for_empty_array(data_list[0]['precip'][k,:,:])) | 
                 (self.check_for_empty_array(data_list[1]['precip'][k,:,:])) ): continue
            

            self.add_title(axes,valid_time,analysis,'Precipitation Difference: {}-{}'.format(exps[0],exps[1]))

            cs = plt.pcolormesh(plons, plats, data_list[0]['precip'][k,:,:].values-data_list[1]['precip'][k,:,:].values,
                                cmap=cmap,
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels, extend='both')
            cb.set_label(r"$mm$", rotation=270)

            fig.canvas.draw()

            figure_name = "{}/PRECIPdiff_{}_{}_{}-{}.png".format(args.output_dir,
                                                    exps[0], exps[1],
                                                    analysis.strftime('%Y%m%d_%H%M'),
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)

        return


    def plotdiff_z(self, args:argparse.Namespace, data_list:list) -> None:
        """Plot orography (surface geopotential)

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments from command line
        data_list : list of xr.Dataset
            Data to plot
        """
        exps = args.experiments.split(",")

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data_list[0]['lon'], data_list[0]['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data_list[0])

        colors = levels_and_colors.z.colors
        levels = [k for k in levels_and_colors.z.levels]

        cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend='both')

        analysis = data_list[0]['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data_list[0]['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if ( (self.check_for_empty_array(data_list[0]['z'][k,:,:])) |
                 (self.check_for_empty_array(data_list[1]['z'][k,:,:])) ): continue

            self.add_title(axes,valid_time,analysis,'Orography Difference: {}-{}'.format(exps[0],exps[1]))

            cs = plt.pcolormesh(plons, plats, data_list[0]['z'][k,:,:]-data_list[1]['z'][k,:,:],
                                cmap=cmap,
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels, extend='both')
            cb.set_label(r"$m$", rotation=270)

            fig.canvas.draw()

            figure_name = "{}/zdiff_{}_{}_{}-{}.png".format(args.output_dir,
                                                    exps[0], exps[1],
                                                    analysis.strftime('%Y%m%d_%H%M'),
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)


        return



    def fig_ax(self, w:int, h:int, **kwargs:dict) -> tuple:
        """Get figure and axes

        Parameters
        ----------
        w : int
            Width of figure
        h : int
            Height of figure

        Returns
        -------
        tuple
            Figure and axes
        """
        fig = plt.figure(figsize=(w, h))
        ax = fig.subplots(1,1, **kwargs)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        return fig, ax


    def add_coastlines(self, ax:plt.subplots, **kwargs:dict) -> tuple:
        """Add coastlines to plot

        Parameters
        ----------
        ax : plt.subplots
            Axes to add coastlines to

        Returns
        -------
        tuple
            extent and coastlines
        """
        # set_extent causes segmentation faults on some old cartopy installations
        extent = ax.set_extent(self.extent, self.data_crs)
        coastline = ax.coastlines(resolution='10m', color=(0.2,0.2,0.2), linewidth=0.7)
        return extent, coastline


    def add_title(self, ax:plt.subplots, validtime:dt.datetime,
                  analysis:dt.datetime, headline:str, **kwargs:dict) -> tuple:
        """Add title to plot

        Parameters
        ----------
        ax : plt.subplots
            axes to add title to
        validtime : dt.datetime
            Valid Time of data
        analysis : dt.datetime
            Analysis Time of data
        headline : str
            Headline to add to title

        Returns
        -------
        tuple
            Title_left, Title_center, Title_right
        """
        title_left = ax.set_title(validtime.strftime('Valid: %Y-%m-%d %H:%Mz'),
                                  fontsize=10, loc='left', y=-0.07)
        title_center = ax.set_title(headline, fontsize=9, loc='center')
        title_right = ax.set_title(analysis.strftime('Analysis: %Y-%m-%d %H:%Mz'),
                                   fontsize=10, loc='right', y=-0.07)
        return title_left, title_center, title_right


    def get_pcolormesh_center_coordinates(self, data:xr.Dataset) -> tuple:
        """Get pcolormesh center coordinates

        Parameters
        ----------
        data : xr.Dataset
            Data holding lon and lat

        Returns
        -------
        tuple
            lons, lats of center of each cell
        """
        # Subtract 1/2 the grid size from both lon and lat arrays
        dlon = (data['lon'][0,1] - data['lon'][0,0]).values
        dlat = (data['lat'][1,0] - data['lat'][0,0]).values
        lons = data['lon'].values - dlon/2
        lats = data['lat'].values - dlat/2
        # Add 1 grid spacing to the right column of lon array and concatenate it
        # as an additional column to the right
        lons = np.c_[ lons, lons[:,-1]+dlon ]
        # Duplicate the bottom row of the lon array and concatenate it to the bottom
        lons = np.r_[ lons, [lons[-1,:]] ]
        # Duplicate the right-most column of lats array and concatenate it on the right
        lats = np.c_[ lats, lats[:,-1] ]
        # Add 1 grid spacing to the bottom row of lat array and concatenate it
        # as an additional row below
        lats = np.r_[ lats, [lats[-1,:]+dlat] ]
        return lons, lats


    def check_for_empty_array(self, da:xr.DataArray) -> bool:
        """Check if array is empty

        Parameters
        ----------
        da : xr.DataArray
            Array of data

        Returns
        -------
        bool
            If empty return True, False otherwise
        """
        empty = False
        flat = da.values.ravel()
        notna = flat[~np.isnan(flat)]
        if len(notna) == 0:
            empty = True
        return empty


class levels_and_colors:
    """Class for levels and colors
    """

    class t2m:
        """Class for t2m levels and colors
        """
        # 21 levels 
        levels=[-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-0.5,0.,0.5,2.,3.,4.,5.,6.,7.,8.,9.,10.0]

        # 22 colors
        colors = [(0.00, 0.00, 0.30), (0.00, 0.00, 0.43), (0.00, 0.00, 0.56), (0.00, 0.00, 0.70),
                  (0.00, 0.00, 0.83), (0.00, 0.00, 0.96), (0.13, 0.13, 1.00),
                  (0.32, 0.32, 1.00), (0.51, 0.51, 1.00), (0.69, 0.69, 1.00),
                  (0.88, 0.88, 1.00), (1.00, 0.93, 0.93),
                  (1.00, 0.74, 0.74), (1.00, 0.55, 0.55), (1.00, 0.36, 0.36),
                  (1.00, 0.18, 0.18), (0.99, 0.00, 0.00), (0.90, 0.00, 0.00),
                  (0.81, 0.00, 0.00), (0.71, 0.00, 0.00), (0.62, 0.00, 0.00), (0.52, 0.00, 0.00)]

    class w10m:
        """Class for w10m levels and colors
        """
        # 21 levels 
        levels=[-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.0]

        # 22 colors
        colors = [(0.00, 0.00, 0.30), (0.00, 0.00, 0.43), (0.00, 0.00, 0.56), (0.00, 0.00, 0.70),
                  (0.00, 0.00, 0.83), (0.00, 0.00, 0.96), (0.13, 0.13, 1.00),
                  (0.32, 0.32, 1.00), (0.51, 0.51, 1.00), (0.69, 0.69, 1.00),
                  (0.88, 0.88, 1.00), (1.00, 0.93, 0.93),
                  (1.00, 0.74, 0.74), (1.00, 0.55, 0.55), (1.00, 0.36, 0.36),
                  (1.00, 0.18, 0.18), (0.99, 0.00, 0.00), (0.90, 0.00, 0.00),
                  (0.81, 0.00, 0.00), (0.71, 0.00, 0.00), (0.62, 0.00, 0.00), (0.52, 0.00, 0.00)]

    class precip:
        """Class for precip levels and colors
        """
        # 21 levels
        levels = [-20.0, -10.0, -7.0, -5.0, -4.0, -3.0, -2.0, -1.0, -0.5, -0.1, 0.0,
                  0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0]

        # 22 colors
        colors = [(0.00, 0.00, 0.30), (0.00, 0.00, 0.43), (0.00, 0.00, 0.56), (0.00, 0.00, 0.70),
                  (0.00, 0.00, 0.83), (0.00, 0.00, 0.96), (0.13, 0.13, 1.00),
                  (0.32, 0.32, 1.00), (0.51, 0.51, 1.00), (0.69, 0.69, 1.00),
                  (0.88, 0.88, 1.00), (1.00, 0.93, 0.93),
                  (1.00, 0.74, 0.74), (1.00, 0.55, 0.55), (1.00, 0.36, 0.36),
                  (1.00, 0.18, 0.18), (0.99, 0.00, 0.00), (0.90, 0.00, 0.00),
                  (0.81, 0.00, 0.00), (0.71, 0.00, 0.00), (0.62, 0.00, 0.00), (0.52, 0.00, 0.00)]

    class z:
        """Class for z levels and colors
        """
        # 21 levels
        levels = [-50.0, -20.0, -10.0, -7.5, -5.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0,
                  0.5, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5, 10.0, 20.0, 50.0]

        # 22 colors
        colors = [(0.00, 0.00, 0.30), (0.00, 0.00, 0.43), (0.00, 0.00, 0.56), (0.00, 0.00, 0.70),
                  (0.00, 0.00, 0.83), (0.00, 0.00, 0.96), (0.13, 0.13, 1.00),
                  (0.32, 0.32, 1.00), (0.51, 0.51, 1.00), (0.69, 0.69, 1.00),
                  (0.88, 0.88, 1.00), (1.00, 0.93, 0.93),
                  (1.00, 0.74, 0.74), (1.00, 0.55, 0.55), (1.00, 0.36, 0.36),
                  (1.00, 0.18, 0.18), (0.99, 0.00, 0.00), (0.90, 0.00, 0.00),
                  (0.81, 0.00, 0.00), (0.71, 0.00, 0.00), (0.62, 0.00, 0.00), (0.52, 0.00, 0.00)]
 
