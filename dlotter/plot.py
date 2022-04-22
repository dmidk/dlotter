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
        
        print('- Plotting', flush=True)

        if 't2m' in parameters:
            if 't2m' not in avail_parameters:
                print('t2m was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_t2m(args, data)

        if 'td2m' in parameters:
            if 'td2m' not in avail_parameters:
                print('td2m was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_td2m(args, data)
        
        if 'w10m' in parameters:
            if 'u10m' not in avail_parameters or 'v10m' not in avail_parameters:
                print('u10m/v10m was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_w10m(args, data)

        if 'precip' in parameters:
            if 'precip' not in avail_parameters:
                print('precip was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_precip(args, data)

        if 'slp' in parameters:
            if 'slp' not in avail_parameters:
                print('slp was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_slp(args, data)

        if 'tcc' in parameters:
            if 'tcc' not in avail_parameters:
                print('tcc was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_tcc(args, data)

        if 'lmhc' in parameters:
            if 'lcc' not in avail_parameters or 'mcc' not in avail_parameters or 'hcc' not in avail_parameters:
                print('lmhc was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_lmhc(args, data)

        if 'snow' in parameters:
            if 'snow' not in avail_parameters:
                print('snow was not found in available parameters: {}, cannot plot'.format(avail_parameters), flush=True)
                sys.exit(1)
            self.plot_snow(args, data)

        return


    def projections(self,args:argparse.Namespace) -> None:
        if args.area == 'dk':
            self.projection = ccrs.AlbersEqualArea(central_longitude=11.0, central_latitude=0.0,
                                    false_easting=0.0, false_northing=0.0,
                                    standard_parallels=(20.0, 50.0), globe=None)
            self.data_crs = ccrs.PlateCarree()
            self.extent = [6, 16, 53, 59]

        if args.area == 'gl':
            # self.projection = ccrs.AlbersEqualArea(central_longitude=-50.0, central_latitude=75.0,
            #                         false_easting=0.0, false_northing=0.0,
            #                         standard_parallels=(-50.0, 70.0), globe=None)
            self.projection = ccrs.TransverseMercator(central_longitude=-50.0, central_latitude=75.0,
                                     false_easting=0.0, false_northing=0.0,)
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

        return


    def plot_t2m(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data['lon'], data['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        if args.dataset == 'DANRAmm':
            colors = ListedColormap(levels_and_colors.t2mdanra.colors)
            levels = [k for k in levels_and_colors.t2mdanra.levels]
            contour_levels = [k for k in levels_and_colors.t2mdanra.contour_levels]
        else:
            colors = ListedColormap(levels_and_colors.t2m.colors)
            levels = [k for k in levels_and_colors.t2m.levels]
            contour_levels = [k for k in levels_and_colors.t2m.contour_levels]
        
        norm = self.color_norm(levels)

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

            cl = plt.contour(clons,clats,data['t2m'][k,:,:].values,
                             colors='black',
                             levels=contour_levels,
                             linewidths=0.7,
                             transform=self.data_crs,
                             transform_first=True)
            
            fig.canvas.draw()       
            
            figure_name = "{}/T2M_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            [cl.collections[k].remove() for k in range(len(cl.collections))]
            print("-- {}".format(figure_name), flush=True)
        
        return


    def plot_td2m(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data['lon'], data['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        colors = ListedColormap(levels_and_colors.t2m.colors)
        levels = [k for k in levels_and_colors.t2m.levels]
        contour_levels = [k for k in levels_and_colors.t2m.contour_levels]
        
        norm = self.color_norm(levels)

        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)
        
        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})
        
        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['td2m'][k,:,:]): continue
       
            self.add_title(axes,valid_time,analysis,'2m Dew Point Temperature')

            cs = plt.pcolormesh(plons, plats, data['td2m'][k,:,:], 
                                cmap=colors, 
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
            cb.set_label(r"$^\circ C$", rotation=270)

            cl = plt.contour(clons,clats,data['td2m'][k,:,:].values,
                             colors='black',
                             levels=contour_levels,
                             linewidths=0.7,
                             transform=self.data_crs)
            
            fig.canvas.draw()       
            
            figure_name = "{}/TD2M_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            [cl.collections[k].remove() for k in range(len(cl.collections))]
            print("-- {}".format(figure_name), flush=True)
        
        return

    
    def plot_w10m(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = np.array(data['lon']), np.array(data['lat'])
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        colors = ListedColormap(levels_and_colors.w10m.colors)
        levels = [k for k in levels_and_colors.w10m.levels]
        
        norm = self.color_norm(levels)

        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        barbs_opts = {'height': 0.7, 'spacing': 0.2}
        bt = self.barbs_thin(clons)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['u10m'][k,:,:]): continue
            if self.check_for_empty_array(data['v10m'][k,:,:]): continue
       
            self.add_title(axes,valid_time,analysis,'10m Wind')

            u = data['u10m'][k,:,:].values
            v = data['v10m'][k,:,:].values

            w10m = np.sqrt(u*u+v*v)

            cs = plt.pcolormesh(plons, plats, w10m, 
                                cmap=colors, 
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
            cb.set_label(r"$m/s$", rotation=270)
            
            
            axes.barbs(clons[::bt,::bt], clats[::bt,::bt], 
                       u[::bt,::bt], v[::bt,::bt], 
                       length=5,
                       sizes=barbs_opts,
                       linewidth=0.95, 
                       transform=self.data_crs,
                       color='#404040')

            fig.canvas.draw()       
            
            figure_name = "{}/W10M_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)
        
        return

    
    def plot_precip(self, args:argparse.Namespace, data:xr.Dataset) -> None:
        # Fix that pcolormesh uses cell lower left corners
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        colors = ListedColormap(levels_and_colors.precip.colors)
        levels = [k for k in levels_and_colors.precip.levels]
        
        norm = self.color_norm(levels)
        
        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['precip'][k,:,:]): continue
            precip = data['precip'][k,:,:].values
       
            self.add_title(axes,valid_time,analysis,'Precipitation')

            cs = plt.pcolormesh(plons, plats, precip, 
                                cmap=colors, 
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
            cb.set_label(r"$mm$", rotation=270)

            fig.canvas.draw()       
            
            figure_name = "{}/PRECIP_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)

        return

    
    def plot_slp(self, args:argparse.Namespace, data:xr.Dataset) -> None:
        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data['lon'], data['lat']

        levels = [k for k in levels_and_colors.slp.levels]
        
        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})

        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['slp'][k,:,:]): continue
            slp = data['slp'][k,:,:].values
       
            self.add_title(axes,valid_time,analysis,'Sea Level Pressure')

            cl = plt.contour(clons,clats,slp,
                             colors='black',
                             levels=levels,
                             linewidths=0.8,
                             transform=self.data_crs)

            fig.canvas.draw()       
            
            figure_name = "{}/SLP_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            [cl.collections[k].remove() for k in range(len(cl.collections))]
            print("-- {}".format(figure_name), flush=True)

        return

    def plot_tcc(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        # Fix that pcolormesh uses cell lower left corners
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        colors = ListedColormap(levels_and_colors.tcc.colors)
        levels = [k for k in levels_and_colors.tcc.levels]        
        norm = self.color_norm(levels)

        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)
        
        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})
        
        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['tcc'][k,:,:]): continue
       
            self.add_title(axes,valid_time,analysis,'Total Cloud Cover')

            cs = plt.pcolormesh(plons, plats, data['tcc'][k,:,:], 
                                cmap=colors, 
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
            cb.set_label(r"$\%$", rotation=270)
            
            fig.canvas.draw()       
            
            figure_name = "{}/TCC_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
            print("-- {}".format(figure_name), flush=True)
        
        return


    def plot_lmhc(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        # Fix that pcolormesh uses cell lower left corners
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        lcc_colors = ListedColormap(levels_and_colors.lcc.colors)
        lcc_levels = [k for k in levels_and_colors.lcc.levels]        
        lcc_norm = self.color_norm(lcc_levels)

        mcc_colors = ListedColormap(levels_and_colors.mcc.colors)
        mcc_levels = [k for k in levels_and_colors.mcc.levels]        
        mcc_norm = self.color_norm(mcc_levels)

        hcc_colors = ListedColormap(levels_and_colors.hcc.colors)
        hcc_levels = [k for k in levels_and_colors.hcc.levels]        
        hcc_norm = self.color_norm(hcc_levels)

        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)
        
        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})
        
        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['lcc'][k,:,:]): continue
            if self.check_for_empty_array(data['mcc'][k,:,:]): continue
            if self.check_for_empty_array(data['hcc'][k,:,:]): continue
       
            self.add_title(axes,valid_time,analysis,'Low Medium High clouds')


            #TODO: pcolormesh with multiple layers
            hcc = data['hcc'][k,:,:].values
            hcc[np.where(hcc<hcc_levels[0])] = np.nan
            cs_hcc = plt.pcolormesh(plons, plats, hcc, 
                                cmap=hcc_colors, 
                                norm=hcc_norm,
                                transform=self.data_crs)

            mcc = data['mcc'][k,:,:].values
            mcc[np.where(mcc<mcc_levels[0])] = np.nan
            cs_mcc = plt.pcolormesh(plons, plats, data['mcc'][k,:,:], 
                                cmap=mcc_colors, 
                                norm=mcc_norm,
                                transform=self.data_crs)

            lcc = data['lcc'][k,:,:].values
            lcc[np.where(lcc<lcc_levels[0])] = np.nan
            cs_lcc = plt.pcolormesh(plons, plats, lcc, 
                                cmap=lcc_colors, 
                                norm=lcc_norm,
                                transform=self.data_crs)

            #cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
            #cb.set_label(r"$\%$", rotation=270)
            
            fig.canvas.draw()       
            
            figure_name = "{}/LMHC_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            #cb.remove()
            
            cs_lcc.remove()
            cs_mcc.remove()
            cs_hcc.remove()

            print("-- {}".format(figure_name), flush=True)
        
        return

    def plot_snow(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        # Fix that pcolormesh uses cell lower left corners
        clons, clats = data['lon'], data['lat']
        plons, plats = self.get_pcolormesh_center_coordinates(data)

        colors = ListedColormap(levels_and_colors.precip.colors)
        levels = [k for k in levels_and_colors.precip.levels]
        
        norm = self.color_norm(levels)

        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)
        
        fig, axes = self.fig_ax(10, 8, subplot_kw={'projection': self.projection})
        
        self.add_coastlines(axes)

        for k in range(self.nt):
            valid_time = data['time'][k].values
            valid_time = dt.datetime.utcfromtimestamp(valid_time.astype(int) * 1e-9)

            if self.check_for_empty_array(data['snow'][k,:,:]): continue
       
            self.add_title(axes,valid_time,analysis,'snow')

            cs = plt.pcolormesh(plons, plats, data['snow'][k,:,:], 
                                cmap=colors,
                                norm=norm,
                                transform=self.data_crs)

            cb = plt.colorbar(cs, fraction=0.046, pad=0.04, ticks=levels)
#            cb = plt.colorbar(cs, fraction=0.046, pad=0.04)
            cb.set_label(r"$mm$", rotation=270)
            
            fig.canvas.draw()       
            
            figure_name = "{}/SNOW_{}-{}.png".format(args.output_dir,
                                                    analysis.strftime('%Y%m%d_%H%M'), 
                                                    valid_time.strftime('%Y%m%d_%H%M'))
            plt.savefig(figure_name)
            cb.remove()
            cs.remove()
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
        # set_extent causes segmentation faults on some old cartopy installations
        extent = ax.set_extent(self.extent, self.data_crs)
        coastline = ax.coastlines(resolution='10m', color=(0.2,0.2,0.2), linewidth=0.7)
        return extent, coastline

    def add_title(self, ax:plt.subplots, validtime:dt.datetime, 
                  analysis:dt.datetime, headline:str, **kwargs:dict) -> tuple:
        title_left = ax.set_title(validtime.strftime('Valid: %Y-%m-%d %H:%Mz'), fontsize=10, loc='left')
        title_center = ax.set_title(headline, fontsize=9, loc='center')
        title_right = ax.set_title(analysis.strftime('Analysis: %Y-%m-%d %H:%Mz'), fontsize=10, loc='right')
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


    def barbs_thin(self, clons:np.array) -> int:
        dx = abs(clons.flatten()[0] - clons.flatten()[1])
        thinner = int(0.5/dx)
        return thinner


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

    class t2mdanra:

        levels=[-10,-6,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,25]

        colors=[(0.5, 0.0, 1.0, 1.0), (0.41919191919191917, 0.12659245357374926, 0.9979866764718844, 1.0), (0.33838383838383834, 0.2511479871810792, 0.9919548128307953, 1.0), (0.25757575757575757, 0.3716624556603275, 0.9819286972627067, 1.0), (0.17676767676767674, 0.4861967361004687, 0.9679487013963562, 1.0), (0.09595959595959591, 0.5929079290546405, 0.9500711177409453, 1.0), (0.015151515151515138, 0.6900790114821119, 0.9283679330160726, 1.0), (0.10606060606060608, 0.8145759520503357, 0.8888354486549235, 1.0), (0.20707070707070718, 0.8959937742913359, 0.8497254299495144, 1.0), (0.3080808080808082, 0.954902241444074, 0.8052702575310585, 1.0), (0.4292929292929294, 0.9938384644612541, 0.7452644496757547, 1.0), (0.5505050505050506, 0.9968547759519424, 0.6785094115571322, 1.0), (0.6515151515151516, 0.9718115683235417, 0.6181589862206052, 1.0), (0.792929292929293, 0.8959937742913359, 0.5272254676105024, 1.0), (0.8737373737373739, 0.8325698546347712, 0.47227107477268254, 1.0), (0.9545454545454546, 0.7557495743542583, 0.41541501300188644, 1.0), (1.0, 0.6900790114821119, 0.3716624556603275, 1.0), (1.0, 0.6181589862206051, 0.32706796331742155, 1.0), (1.0, 0.5406408174555974, 0.2817325568414296, 1.0), (1.0, 0.4582265217274105, 0.23575893550942728, 1.0), (1.0, 0.31203344569848696, 0.1580013959733498, 1.0), (1.0, 0.22031053278654036, 0.11083819990101088, 1.0), (1.0, 0.15800139597334945, 0.07924995685678822, 1.0), (1.0, 0.12659245357374938, 0.06342391965656456, 1.0), (1.0, 0.09505604330418288, 0.0475819158237424, 1.0), (1.0, 0.0634239196565644, 0.031727933498067594, 1.0), (1.0, 1.2246467991473532e-16, 6.123233995736766e-17, 1.0)]



#        colors = [(0.00, 0.00, 0.40),(0.00, 0.20, 0.60), (0.20, 0.40, 0.80), (0.00, 0.60, 1.00),
#                      (0.00, 0.80, 1.00),(0.00, 1.00, 1.00),(0.60, 1.00, 0.60),(0.60, 1.00, 0.40),
#                      (0.60, 1.00, 0.20),(0.00, 1.00, 0.00),(0.16, 0.85, 0.00),(0.08, 0.78, 0.00),
#                      (0.80, 0.80, 0.00),(0.90, 0.90, 0.00),(1.00, 1.00, 0.00),(1.00, 0.80, 0.00),
#                      (1.00, 0.70, 0.00),(1.00, 0.60, 0.00),(1.00, 0.50, 0.00),(1.00, 0.40, 0.00),
#                      (0.98, 0.3, 0.00),(0.98, 0.2, 0.00),(0.4, 0.00, 0.00),(0.80, 0.40, 0.00),
#                      (0.80, 0.60, 0.00),(0.96, 0.31, 0.00),(0.96, 0.22, 0.00),(0.96, 0.12, 0.00),
#                      (0.96, 0.00, 0.00),(0.95, 0.00, 0.00)]

        contour_levels = [-10, -5, 0, 5, 10, 15, 20, 25]


    class w10m:
        levels=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]

        colors =  [(0.00,0.00,1.00),(0.00,0.00,1.00),(0.04,0.25,1.00),(0.08,0.37,1.00),
                       (0.10,0.49,1.00),(0.14,0.60,1.00),(0.16,0.73,1.00),(0.20,0.88,0.93),
                       (0.21,0.90,0.53),(0.21,0.91,0.33),(0.21,0.91,0.10),(0.22,0.96,0.00),
                       (0.24,1.00,0.00),(1.00,1.00,0.00),(0.99,0.87,0.00),(0.98,0.76,0.00),
                       (0.97,0.64,0.00),(0.97,0.53,0.00),(0.85,0.42,0.00),(0.85,0.32,0.00),
                       (0.80,0.22,0.03),(0.70,0.13,0.06),(0.62,0.05,0.07),(0.58,0.00,0.08),
                       (0.62,0.00,0.15),(0.71,0.00,0.33),(0.96,0.00,0.54),(0.96,0.00,0.77),
                       (0.96,0.00,1.00),(0.82,0.00,0.90),(0.70,0.00,0.77),(0.60,0.00,0.65),
                       (0.49,0.00,0.54),(0.39,0.00,0.44),(0.31,0.00,0.33)]

    class precip:
        # 17 levels
        levels = [0.0, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
        
        # 16 colors
        colors = [(1.00, 1.00, 1.00), (0.84, 0.69, 0.63), (0.72, 0.56, 0.51), (0.80, 0.86, 1.00),
                  (0.65, 0.73, 1.00), (0.49, 0.63, 1.00), (0.25, 0.52, 0.97), (0.00, 0.00, 0.98),
                  (0.13, 0.57, 0.09), (0.20, 0.71, 0.18), (0.41, 0.80, 0.36), (0.68, 0.99, 0.36),
                  (1.00, 1.00, 0.08), (0.99, 0.57, 0.05), (0.86, 0.00, 0.02), (0.68, 0.00, 0.02)]


    class slp:
        # 17 levels
        levels = np.arange(950,1050,2)


    class tcc:
        # 9 levels
        levels=[0.0, 1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 8/8]

        # 8 Colors
        colors = [(1.00, 1.00, 1.00), (0.95, 0.95, 0.95), (0.90, 0.90, 0.90), (0.80, 0.80, 0.80),
                      (0.70, 0.70, 0.70), (0.50, 0.50, 0.50), (0.40, 0.40, 0.40), (0.30, 0.30, 0.30)]
        
    class lcc:
        # 9 levels
        levels=[2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 8/8]

        # 8 Colors
        colors = [(0.70, 0.70, 0.70), (0.60, 0.60, 0.60),
                      (0.50, 0.50, 0.50), (0.40, 0.40, 0.40), (0.30, 0.30, 0.30), (0.20, 0.20, 0.20)]

    class mcc:
        # 9 levels
        levels=[2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 8/8]

        # 8 Colors
        # colors = [(0.00, 0.40, 0.00), (0.00, 0.50, 0.00),
        #               (0.00, 0.60, 0.00), (0.00, 0.70, 0.00), (0.00, 0.80, 0.00), (0.00, 0.90, 0.00)]
        colors = [(0.00, 0.90, 0.00), (0.00, 0.80, 0.00),
                      (0.00, 0.70, 0.00), (0.00, 0.60, 0.00), (0.00, 0.50, 0.00), (0.00, 0.40, 0.00)]

    class hcc:
        # 9 levels
        levels=[2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 8/8]

        # 8 Colors
        # colors = [(0.00, 0.00, 0.40), (0.00, 0.00, 0.50),
        #               (0.00, 0.00, 0.60), (0.00, 0.00, 0.70), (0.00, 0.00, 0.80), (0.00, 0.00, 0.90)]
        colors = [(0.00, 0.00, 0.90), (0.00, 0.00, 0.80),
                      (0.00, 0.00, 0.70), (0.00, 0.00, 0.60), (0.00, 0.00, 0.50), (0.00, 0.00, 0.40)]
