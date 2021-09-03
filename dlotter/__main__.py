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


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


if __name__ == '__main__':
    parent_parser = MyParser(
        description='Plot data quick and dirty from NWP output',
        formatter_class=ArgumentDefaultsHelpFormatter)
    
    subparser = parent_parser.add_subparsers(dest="cmd")

    parser_prepare = subparser.add_parser('plot', help='Plot NWP output')

    parser_prepare.add_argument('-p',
                    '--parameters',
                    metavar='PARAMETERS',
                    type=str,
                    help='Parameters to plot. Seperate with ":", eg: "t2m:w10m".',
                    required=True)
    
    parser_prepare.add_argument('-f',
                    '--filetype',
                    metavar='FILETYPE',
                    type=str,
                    help='What filetype are we using? (Options are: grib2)',
                    default='grib2',
                    required=False)
    
    parser_prepare.add_argument('-d',
                    '--directory',
                    type=str,
                    help='directory to read data from',
                    default='.')
    
    parser_prepare.add_argument('--prefix',
                    type=str,
                    help='Set to prefix of files if any',
                    default='',
                    required=False)

    parser_prepare.add_argument('--postfix',
                    type=str,
                    help='Set to postfix of files if any',
                    default='',
                    required=False)
    
    parser_prepare.add_argument('-o',
                    '--output-dir',
                    metavar='OUTDIR',
                    type=str,
                    help='Directory to place output into',
                    default='.',
                    required=False)

    parser_prepare.add_argument('-l',
                    '--limit-files',
                    metavar='LIMIT',
                    type=int,
                    help='Only use the first LIMIT files. If set to 0, not limit is used. If Limit > 0, files will be sorted by name first',
                    default=0,
                    required=False)

    parser_prepare.add_argument('--verbose',
                    action='store_true',
                    help='Verbose output', 
                    default=False)


    if len(sys.argv)==1:
        parent_parser.print_help()
        sys.exit(2)

    args = parent_parser.parse_args()

    if args.verbose:
        print('---- Input Arguments ----', flush=True)
        for p in args._get_kwargs():
                print("{}: {}".format(p[0], p[1]), flush=True)
        print('---- --------------- ----', flush=True)

    if args.cmd == 'plot':
        prepwork = prepare(args)
        files_to_read = prepwork.files_to_read
        
        if args.filetype == 'grib2':
            grib2Read(args, files_to_read)
        else:
            print('Filetype: "{}", not supported.'.format(args.filetype), flush=True)
            sys.exit(1)