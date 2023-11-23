#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter
"""

__author__ = "K. Hintz"
__copyright__ = "Danish Meteorological Institute"

__license__ = "MIT"
__version__ = "0.0.2.2"
__maintainer__ = "K. Hintz"
__email__ = "kah@dmi.dk"
__status__ = "Development"

import sys
import os
sys.path.insert(0, os.path.abspath('./dlotter/'))
import argparse

from .prepare import prepare
from .read import grib2Read, netcdf2read
from .plot import plot
from .plotdiff import plotdiff
from .arguments import arguments
from .meteogram import meteogram
from .comeps import comeps
from .utils import s3util
from dmit import ostools

class MyParser(argparse.ArgumentParser):
    """Parser for dlotter arguments

    Parameters
    ----------
    argparse : argparse.ArgumentParser
        Parser for dlotter arguments
    """
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


if __name__ == '__main__':

    modargs = arguments()
    args = modargs.get_args(sys.argv)

    if args.verbose:
        print('---- Input Arguments ----', flush=True)
        for p in args._get_kwargs():
                print("{}: {}".format(p[0], p[1]), flush=True)
        print('---- --------------- ----', flush=True)

    # Check if output directory exist. If not try to create it.
    if not ostools.does_dir_exist(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)


    if args.cmd == 'plot':
        prepwork = prepare(args)
        files_to_read = prepwork.files_to_read

        if args.verbose:
            print('Found the following files:', flush=True)
            print(files_to_read, flush=True)

        if args.filetype == 'grib2':
            datareader = grib2Read(args, files_to_read)
            data = datareader.data
        elif args.filetype == 'nc':
            datareader = netcdf2read(args, files_to_read)
            data = datareader.data
        else:
            print('Filetype: "{}", not supported.'.format(args.filetype), flush=True)
            sys.exit(1)

        plotwork = plot(args, data)

        if args.s3:
            plot_files = ostools.find_files(args.output_dir, postfix='.png', recursive=False, fullpath=True)
            s3tool = s3util()
            s3tool.send_files_to_s3(plot_files, args.s3path) # args.s3path: 'ig', 'dini' etc.


    if args.cmd == 'epsmeteogram':
        eps = comeps(args)
        data = eps.data

        plotwork = meteogram(args, data)


    if args.cmd == 'plotdiff':
        if len(args.directory.split(",")) != 2:
            print('Number of specified directories are {:d}, which is not supported. \
                   Must be 2.'.format(len(args.directory.split(","))), flush=True)
            sys.exit(1)

        prepwork = prepare(args)
        files_to_read = prepwork.files_to_read

        if args.verbose:
            print('Found the following files:', flush=True)
            print(files_to_read, flush=True)

        if args.filetype == 'grib2':
            datareader = grib2Read(args, files_to_read)
            data = datareader.data
        elif args.filetype == 'nc':
            datareader = netcdf2read(args, files_to_read)
            data = datareader.data
        else:
            print('Filetype: "{}", not supported.'.format(args.filetype), flush=True)
            sys.exit(1)

        plotwork = plotdiff(args, data)

