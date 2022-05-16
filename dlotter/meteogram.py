#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.meteogram
Called from dlotter.__main__
Plots meteograms
"""
import sys
import argparse
import xarray as xr
import numpy as np
import datetime as dt

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


class meteogram:


    def __init__(self, args:argparse.Namespace, data:xr.Dataset) -> None:

        print('- Plotting', flush=True)

        self.define_icons()

        self.plot_epsmeteogram(args, data)

        return


    def define_icons(self) -> None:
        #Day
        self.sun          = 'dlotter/icons/0.png'
        self.partcloud    = 'dlotter/icons/1.png'
        self.cloudy       = 'dlotter/icons/2.png'
        self.haze         = 'dlotter/icons/3.png'

        self.rain_shower        = 'dlotter/icons/101.png'
        self.rain               = 'dlotter/icons/102.png'
        self.heavy_rain_shower  = 'dlotter/icons/111.png'
        self.heavy_rain         = 'dlotter/icons/112.png'

        self.sleet_shower       = 'dlotter/icons/201.png'
        self.sleet              = 'dlotter/icons/202.png'
        self.heavy_sleet_shower = 'dlotter/icons/211.png'
        self.heavy_sleet        = 'dlotter/icons/212.png'

        self.snow_shower        = 'dlotter/icons/301.png'
        self.snow               = 'dlotter/icons/302.png'
        self.heavy_snow_shower  = 'dlotter/icons/311.png'
        self.heavy_snow         = 'dlotter/icons/312.png'

        # Night
        self.moon               = 'dlotter/icons/1000.png'
        self.partcloud_night    = 'dlotter/icons/1001.png'
        self.cloudy_night       = 'dlotter/icons/1002.png'
        self.haze_night         = 'dlotter/icons/1003.png'

        self.rain_shower_night        = 'dlotter/icons/1010.png'
        self.rain_night               = 'dlotter/icons/1020.png'
        self.heavy_rain_shower_night  = 'dlotter/icons/1110.png'
        self.heavy_rain_night         = 'dlotter/icons/1120.png'

        self.sleet_shower_night       = 'dlotter/icons/2010.png'
        self.sleet_night              = 'dlotter/icons/2020.png'
        self.heavy_sleet_shower_night = 'dlotter/icons/2110.png'
        self.heavy_sleet_night        = 'dlotter/icons/2120.png'

        self.snow_shower_night        = 'dlotter/icons/3010.png'
        self.snow_night               = 'dlotter/icons/3020.png'
        self.heavy_snow_shower_night  = 'dlotter/icons/3110.png'
        self.heavy_snow_night         = 'dlotter/icons/3120.png'

        return
    

    def plot_epsmeteogram(self, args:argparse.Namespace, data:xr.Dataset) -> None:
        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)
        
        x = np.arange(args.files_per_member)
        y = np.arange(args.members)
        xx, yy = np.meshgrid(x,y)

        locations = args.latlon.split(':')

        l=0
        for loc in locations:

            fig, ax = plt.subplots(figsize=(14,8))
            self.add_title(ax, analysis, "location: {}".format(loc))
            
            k = 0
            for x0, y0 in zip(xx.flatten(), yy.flatten()):
                #filenumber = x0 #aka time-dimension
                #member = y0
                
                cloud = data['tcc'][x0,y0,l]
                precip = data['precip'][x0,y0,l]
                #rain = data['rain'][x0,y0,:,:].values.flatten()[location_idx]
                visibility = data['visibility'][x0,y0,l]

                precip_type = 'ra' # TODO fix this
                
                symbol = self.get_weather_symbol(cloud, precip, precip_type, visibility)
                #path = self.sun

                path = symbol
                k+=1

                ab = AnnotationBbox(self.getImage(path), (x0, y0), frameon=False)
                ax.add_artist(ab)

            ylabel = ['mbr{}'.format(k) for k in y]
            plt.yticks(y, ylabel)

            xlabel = [(analysis+dt.timedelta(hours=int(k))).strftime('%a\n%Hz') for k in x]
            plt.xticks(x[::2], xlabel[::2])
            plt.margins(x=0, y=10)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='both', which='major', pad=25)

            plt.show()
            l+=1

        return


    def getImage(self, path):
        return OffsetImage(plt.imread(path), zoom=0.1)

    def fig_ax(self, w:int, h:int, **kwargs:dict) -> tuple:
        fig = plt.figure(figsize=(w, h))
        ax = fig.subplots(1,1, **kwargs)
        return fig, ax

    def add_title(self, ax:plt.subplots, analysis:dt.datetime, headline:str, **kwargs:dict) -> tuple:
        title_center = ax.set_title(headline, fontsize=9, loc='center')
        title_right = ax.set_title(analysis.strftime('Analysis: %Y-%m-%d %H:%M'), fontsize=10, loc='right', pad=20)
        return title_center, title_right



    def get_weather_symbol(self, cloud:float, precip:float, precip_type:str, visibility:float) -> str:

        symbol = 0
        night = False

        if precip >= 2.0: # Heavy Rain
            if precip_type == 'ra':
                symbol+=110
            elif precip_type == 'sl':
                symbol+=210
            elif precip_type == 'sn':
                symbol+=310
        elif precip > 0.1: # Light Rain
            if precip_type == 'ra':
                symbol+=100
            elif precip_type == 'sl':
                symbol+=200
            elif precip_type == 'sn':
                symbol+=300
        elif precip <= 0.1: # No Rain
            symbol += 0


        if visibility <= 200:
            symbol+=3
            if precip > 0.1: # If Rain + Low visibility
                symbol -= 1
        elif cloud > 0.8: # Cloudy
            symbol+=2
        elif cloud > 0.2: # Broken
            symbol+=1
        elif cloud <= 0.2: # Clear
            symbol+=0

        icon = 'dlotter/icons/{}.png'.format(symbol)

        return icon