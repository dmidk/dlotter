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

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

class plot:

    def __init__(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        self.projections(args)
        
        parameters = args.parameters
        avail_parameters = list(data.data_vars)
        
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

        print('Plotting')

        fig, axes = self.fig_ax('T2m', 10, 8, subplot_kw={'projection': self.projection})
        ax.set_extent(self.extent, self.data_crs)
        ax.coastlines(resolution='10m', color=(0.2,0.2,0.2))
        
        fig.canvas.draw()


        plt.show()
        
        return
            # fig, ax = plt.subplots(1, 1, figsize=(10,8),
            #                        subplot_kw={'projection': projection})

            # ax = plt.axes(projection=projection)
            # #ax = plt.axes(projection=osm_img.crs)

    def fig_ax(self, title:str, w:int, h:int, **kwargs:dict) -> tuple:
        fig = plt.figure(figsize=(w, h))
        ax = fig.subplots(1,1, **kwargs)
        ax.set_title(title)
        return fig, ax


    def add_contour(ax, X, Y, data, levels, **kwargs) -> None:
        cs = ax.contour(X, Y, data, levels, **kwargs)
        ax.clabel(cs, fmt='%2.0f', inline=True, fontsize=10)
        return