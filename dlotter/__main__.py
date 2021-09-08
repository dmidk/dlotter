#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter
"""

__author__ = "K. Hintz"
__copyright__ = "Danish Meteorological Institute"

__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "K. Hintz"
__email__ = "kah@dmi.dk"
__status__ = "Development"

import sys
import os
sys.path.insert(0, os.path.abspath('./dlotter/'))
import argparse
from argparse import ArgumentDefaultsHelpFormatter

from .prepare import prepare
from .read import grib2Read
from .plot import plot
from .arguments import arguments


class MyParser(argparse.ArgumentParser):
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

    if args.cmd == 'plot':
        prepwork = prepare(args)
        files_to_read = prepwork.files_to_read
        
        if args.filetype == 'grib2':
            datareader = grib2Read(args, files_to_read)
            data = datareader.data
        else:
            print('Filetype: "{}", not supported.'.format(args.filetype), flush=True)
            sys.exit(1)

        plotwork = plot(args, data)