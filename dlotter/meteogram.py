#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.meteogram
Called from dlotter.__main__
Plots meteograms
"""
import argparse
import xarray as xr
import numpy as np
import datetime as dt

import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


class meteogram:
    """Class for plotting meteograms
    """


    def __init__(self, args:argparse.Namespace, data:xr.Dataset) -> None:
        """Constructor for meteogram class

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments to dlotter from command line
        data : xr.Dataset
            Data to plot
        """

        print('- Plotting', flush=True)

        self.define_icons()
        self.plot_epsmeteogram(args, data)

        print('- Done Plotting', flush=True)

        return


    def define_icons(self) -> None:
        """Define relative paths to icons
        """
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
        """Plot EPS meteogram with weather icons

        Parameters
        ----------
        args : argparse.Namespace
            Input to dlotter from command line
        data : xr.Dataset
            Data to plot
        """
        analysis = data['time'][0].values
        analysis = dt.datetime.utcfromtimestamp(analysis.astype(int) * 1e-9)

        x = np.arange(args.files_per_member)
        y = np.arange(args.members)
        xx, yy = np.meshgrid(x,y)

        locations = args.latlon.split(':')

        l=0
        for loc in locations:

            _, ax = plt.subplots(figsize=(14,8))
            self.add_title(ax, analysis, "location: {}".format(loc))

            k = 0
            for x0, y0 in zip(xx.flatten(), yy.flatten()):
                #filenumber = x0 #aka time-dimension
                #member = y0

                cloud = data['tcc'][x0,y0,l]
                precip = data['precip'][x0,y0,l]
                precip_solid = data['precip_solid'][x0,y0,l]
                visibility = data['visibility'][x0,y0,l]

                precip_type = self.get_precip_type(precip, precip_solid)

                night = data['night'][x0,l]

                symbol = self.get_weather_symbol(cloud, precip, precip_type, visibility, night)

                path = symbol
                k+=1

                ab = AnnotationBbox(self.getImage(path), (x0, y0), frameon=False)
                ax.add_artist(ab)

            ylabel = ['mbr{:02d}'.format(k) for k in y]
            plt.yticks(y, ylabel)

            xlabel = [(analysis+dt.timedelta(hours=int(k))).strftime('%a\n%Hz') for k in x]
            plt.xticks(x[::2], xlabel[::2])
            #plt.margins(x=0, y=10)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='both', which='major', pad=25)

            plt.savefig('{}/epsplot_{}.png'.format(args.output_dir, l))
            l+=1

        return


    def getImage(self, path:str) -> OffsetImage:
        """Read image from file and return it as matplotlib.offsetbox.OffsetImage

        Parameters
        ----------
        path : str
            Path to image file

        Returns
        -------
        OffsetImage
            Image object
        """
        return OffsetImage(plt.imread(path), zoom=0.1)

    def fig_ax(self, w:int, h:int, **kwargs:dict) -> tuple:
        """Create figure and axes

        Parameters
        ----------
        w : int
            Width of figure
        h : int
            Height of figure

        Returns
        -------
        tuple
            Tuple of figure and axis
        """
        fig = plt.figure(figsize=(w, h))
        ax = fig.subplots(1,1, **kwargs)
        return fig, ax


    def add_title(self, ax:plt.subplots, analysis:dt.datetime, headline:str,
                  **kwargs:dict) -> tuple:
        """Add a title to the plot

        Parameters
        ----------
        ax : plt.subplots
            Subplot axes
        analysis : dt.datetime
            Analysis time of the forecast
        headline : str
            Headline to add to plot

        Returns
        -------
        tuple
            Title for center position and right position
        """
        title_center = ax.set_title(headline, fontsize=9, loc='center')
        title_right = ax.set_title(analysis.strftime('Analysis: %Y-%m-%d %H:%M'),
                                   fontsize=10, loc='right', pad=20)
        return title_center, title_right



    def get_weather_symbol(self, cloud:float, precip:float, precip_type:str,
                           visibility:float, night:float) -> str:
        """Get weather symbol

        Parameters
        ----------
        cloud : float
            Cloud fraction
        precip : float
            Amount of precipitation
        precip_type : str
            Type of precipitation
        visibility : float
            Visibility in meters
        night : float
            1 if night, 0 otherwise

        Returns
        -------
        str
            Path to icon
        """

        symbol = 0

        is_night = False
        if night==1: is_night = True

        bool_precip = False


        if precip >= 2.0: # Heavy Rain
            if precip_type == 'ra':
                symbol+=110
            elif precip_type == 'sl':
                symbol+=210
            elif precip_type == 'sn':
                symbol+=310
            bool_precip = True
        elif precip > 0.1: # Light Rain
            if precip_type == 'ra':
                symbol+=100
            elif precip_type == 'sl':
                symbol+=200
            elif precip_type == 'sn':
                symbol+=300
            bool_precip = True
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
        elif cloud <= 0.2 and not bool_precip: # Clear
            symbol+=0
        elif cloud <= 0.2 and bool_precip: # Clear but heavy shower
            symbol+=1

        if is_night: symbol+=1000

        icon = 'dlotter/icons/{}.png'.format(symbol)
        return icon


    def get_precip_type(self, precip:float, precip_solid:float) -> str:
        """Get type of precipitation

        Parameters
        ----------
        precip : float
            Amount of precipitation
        precip_solid : float
            Amount of solid precipitation

        Returns
        -------
        str
            Type of precipitation (ra, sn, sl)
        """

        if precip == 0:
            precip_type = 'ra' # Default value, in this case it will not be used
        elif precip != 0 and precip_solid == precip: # Everything is snow
            precip_type = 'sn'
        elif precip != 0 and precip_solid > 0 and precip > precip_solid:
            # Both rain and snow, choose sleet
            precip_type = 'sl'
        else:
            precip_type = 'ra'

        return precip_type
