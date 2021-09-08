#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.plot
Called from dlotter.__main__
"""
import sys
import argparse
import xarray as xr
import numpy as np
import datetime as dt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import cartopy.crs as ccrs
import cartopy.feature as cfeature


class plot:

    def __init__(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        self.projections(args)
        
        parameters = args.parameters
        avail_parameters = list(data.data_vars)

        self.nt = data.dims['time']
        
        if 't2m' in parameters:
            if 't2m' not in avail_parameters:
                print('t2m was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_t2m(args, data)

        return


    def projections(self,args:argparse.Namespace) -> None:
        if args.area == 'dk':
            self.projection = ccrs.AlbersEqualArea(central_longitude=11.0, central_latitude=0.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [6, 16, 53, 59]

        return


    def plot_t2m(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        print('- Plotting', flush=True)
        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data['lon'], data['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        colors = ListedColormap(levels_and_colors.t2m.colors)
        levels = [k for k in levels_and_colors.t2m.levels]
        contour_levels = [k for k in levels_and_colors.t2m.contour_levels]
        
        norm = self.color_norm(levels)

        data = data.sortby('time')
        
        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['t2m'][k,:,:]): continue
       
            self.add_title(axes,valid_time,analysis,'2m Temperature')

            cs = plt.pcolormesh(plons, plats, data['t2m'][k,:,:], 
                                cmap=colors, 
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
            cb.set_label(r"$^\circ C$", rotation=270)
            
            # cc = plt.contour(clons, clats, data['t2m'][k,:,:], 
            #                  levels=contour_levels, 
            #                  transform=self.data_crs, 
            #                  colors='black',
            #                  linewidths=0.5)

            fig.canvas.draw()       
            
            figure_name = "{}/T2M_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            #cc.remove()
            print("-- {}".format(figure_name), flush=True)
        
        return


    def fig_ax(self, w:int, h:int, **kwargs:dict) -> tuple:
        fig = plt.figure(figsize=(w, h))
        ax = fig.subplots(1,1, **kwargs)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        return fig, ax


    def add_contour(self, ax, X, Y, data, levels, **kwargs) -> None:
        cc = plt.contour(X, Y, data, levels, **kwargs)
        ax.clabel(cc, fmt='%2.0f', inline=True, fontsize=10)
        return cc


    def add_coastlines(self, ax:plt.subplots, **kwargs:dict) -> tuple:
        extent = ax.set_extent(self.extent, self.data_crs)
        coastline = ax.coastlines(resolution='10m', color=(0.2,0.2,0.2), linewidth=0.7)
        return extent, coastline

    def add_title(self, ax:plt.subplots, validtime:dt.datetime, 
                  analysis:dt.datetime, headline:str, **kwargs:dict) -> tuple:
        title_left = ax.set_title(validtime.strftime('Valid: %Y-%m-%d %H:%M'), fontsize=10, loc='left')
        title_center = ax.set_title(headline, fontsize=9, loc='center')
        title_right = ax.set_title(analysis.strftime('Analysis: %Y-%m-%d %H:%M'), fontsize=10, loc='right')
        return title_left, title_center, title_right


    def get_pcolormesh_center_coordinates(self, data:xr.Dataset) -> tuple:
        # Subtract 1/2 the grid size from both lon and lat arrays
        dlon = (data['lon'][0,1] - data['lon'][0,0]).values
        dlat = (data['lat'][1,0] - data['lat'][0,0]).values
        lons = data['lon'].values - dlon/2
        lats = data['lat'].values - dlat/2
        # Add 1 grid spacing to the right column of lon array and concatenate it as an additional column to the right
        lons = np.c_[ lons, lons[:,-1]+dlon ]
        # Duplicate the bottom row of the lon array and concatenate it to the bottom
        lons = np.r_[ lons, [lons[-1,:]] ]
        # Duplicate the right-most column of lats array and concatenate it on the right
        lats = np.c_[ lats, lats[:,-1] ]
        # Add 1 grid spacing to the bottom row of lat array and concatenate it as an additional row below
        lats = np.r_[ lats, [lats[-1,:]+dlat] ]
        return lons, lats


    def color_norm(self, levels:list) -> matplotlib.colors.BoundaryNorm:
        len_lab = len(levels)
        norm_bins = np.sort([*levels])
        norm = matplotlib.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
        return norm


    def check_for_empty_array(self, da:xr.DataArray) -> bool:

        empty = False
        flat = da.values.ravel()
        notna = flat[~np.isnan(flat)]
        if len(notna) == 0:
            empty = True

        return empty


class levels_and_colors:

    class t2m:
        levels=[-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,22,24,26,28,30,32,34,36,38,40,42]

        colors = [(0.14, 0.00, 0.15),(0.31, 0.00, 0.33), (0.49, 0.00, 0.54), (0.71, 0.00, 0.77),
                      (0.93, 0.00, 1.00),(0.70, 0.00, 1.00),(0.49, 0.00, 1.00),(0.30, 0.00, 1.00),
                      (0.13, 0.00, 1.00),(0.00, 0.00, 1.00),(0.05, 0.25, 1.00),(0.11, 0.49, 1.00),
                      (0.16, 0.73, 1.00),(0.22, 1.00, 1.00),(0.24, 1.00, 0.71),(0.24, 1.00, 0.45),
                      (0.24, 0.94, 0.19),(0.16, 0.85, 0.00),(0.08, 0.78, 0.00),(1.00, 1.00, 0.00),
                      (0.98, 0.88, 0.00),(0.98, 0.76, 0.00),(0.97, 0.64, 0.00),(0.97, 0.53, 0.00),
                      (0.96, 0.42, 0.00),(0.96, 0.31, 0.00),(0.96, 0.22, 0.00),(0.96, 0.12, 0.00),
                      (0.96, 0.00, 0.00),(0.95, 0.00, 0.00),(0.95, 0.00, 0.13),(0.93, 0.00, 0.33),
                      (0.96, 0.00, 0.54),(0.96, 0.00, 0.77)]

        contour_levels = [-10, -5, 0, 5, 10, 15, 20, 25]