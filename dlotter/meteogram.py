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

#from pywaffle import Waffle

class meteogram:


    #def __init__(self, args:argparse.Namespace, data:xr.Dataset) -> None:
    def __init__(self, args:argparse.Namespace) -> None:

        print('- Plotting', flush=True)

        self.icon_paths = [
                        'dlotter/icons/0.png',
                        'dlotter/icons/1.png',
                        'dlotter/icons/2.png',
                        'dlotter/icons/102.png']

        self.plot_epsmeteogram(args)

        

        return

    def plot_epsmeteogram(self, args:argparse.Namespace) -> None:

        # # read image file
        # im='dlotter/icons/0.png' #use absolute path instead
        # img=plt.imread(im)

        # fig, axes = self.fig_ax(10, 8)
        # analysis = dt.datetime(2021,9,15,12,0)
        # self.add_title(axes,analysis,'DMI Lyngbyvej')

        # x = np.arange(0,24)
        # y = np.arange(0,20)
        # xx, yy = np.meshgrid(x,y)

        # frame_height=1
        # for ix in x:
        #     x_start=ix-0.5
        #     y_start=-.5
        #     axes.imshow(img,extent=[x_start,x_start+frame_height,y_start,y_start+frame_height])
        # #     print(ix,iy)
        # # x_start=-.5
        # # y_start=-.5
        # # axes.imshow(img,extent=[x_start,x_start+frame_height,y_start,y_start+frame_height])


        # plt.plot(xx, yy, linestyle="None", marker="None")

        # ylabel = ['mbr{}'.format(k) for k in y]
        # plt.yticks(y, ylabel)

        # xlabel = [(analysis+dt.timedelta(hours=int(k))).strftime('%a\n%Hz') for k in x]
        # plt.xticks(x[::2], xlabel[::2])

        # plt.show()



        # fig = plt.figure(
        # FigureClass=Waffle,
        # rows=5,
        # values=[30, 16, 3],
        # colors=["#ff00bf", "#4384FF", "#C0C0C0"],
        # icons=['sun', 'cloud-showers-heavy', 'snowflake'], #custom icons instead???????
        # font_size=20,
        # icon_style='solid',
        # icon_legend=False,
        # )
       # plt.show()

        # x = [0,1,2,3,4]
        # y = [0,1,2,3,4]

        analysis = dt.datetime(2021,9,15,12,0)

        x = np.arange(0,24)
        y = np.arange(0,20)

        xx, yy = np.meshgrid(x,y)

        fig, ax = plt.subplots(figsize=(10,8))
        self.add_title(ax,analysis,'DMI Lyngbyvej')
        ax.scatter(xx, yy) 

        k = 0
        for x0, y0 in zip(xx.flatten(), yy.flatten()):

            if k%2==0:
                path = self.icon_paths[0]
            elif k%5==0:
                path = self.icon_paths[1]
            elif k%3==0:
                path = self.icon_paths[3]
            else:
                path = self.icon_paths[2]
            k+=1

            #print(x0,y0,k, path)
            ab = AnnotationBbox(self.getImage(path), (x0, y0), frameon=False)
            ax.add_artist(ab)

        ylabel = ['mbr{}'.format(k) for k in y]
        plt.yticks(y, ylabel)

        xlabel = [(analysis+dt.timedelta(hours=int(k))).strftime('%a\n%Hz') for k in x]
        plt.xticks(x[::2], xlabel[::2])

        plt.show()

    

        return

    def getImage(self, path):
        return OffsetImage(plt.imread(path), zoom=0.1)

    def fig_ax(self, w:int, h:int, **kwargs:dict) -> tuple:
        fig = plt.figure(figsize=(w, h))
        ax = fig.subplots(1,1, **kwargs)
        return fig, ax

    def add_title(self, ax:plt.subplots, analysis:dt.datetime, headline:str, **kwargs:dict) -> tuple:
        title_center = ax.set_title(headline, fontsize=9, loc='center')
        title_right = ax.set_title(analysis.strftime('Analysis: %Y-%m-%d %H:%M'), fontsize=10, loc='right')
        return title_center, title_right