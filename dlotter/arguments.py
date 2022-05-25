#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.arguments
Called from dlotter.__main__
"""
import sys
import argparse
from argparse import ArgumentDefaultsHelpFormatter


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


class arguments:
    """Holds all the arguments for dlotter
    """

    def __init__(self) -> None:
        """Constructor for arguments class
        """

        return


    def get_args(self, sysargs):
        """Get arguments from command line

        Parameters
        ----------
        sysargs : sys.argv
            Arguments given to dlotter on command line

        Returns
        -------
        argparse.Namespace
            Namespace with all the arguments
        """

        parent_parser = MyParser(
            description='Plot data quick and dirty from NWP output',
            formatter_class=ArgumentDefaultsHelpFormatter)

        subparser = parent_parser.add_subparsers(dest="cmd")

        # Parser for NWP Deterministic
        parser_plot = subparser.add_parser('plot', help='Plot NWP output')

        parser_plot.add_argument('-p',
                        '--parameters',
                        metavar='PARAMETERS',
                        type=str,
                        help='Parameters to plot. Seperate with ":",\
                            eg: "t2m:w10m:precip:slp:td2m:tcc:lmhc".',
                        required=True)

        parser_plot.add_argument('-f',
                        '--filetype',
                        metavar='FILETYPE',
                        type=str,
                        help='What filetype are we using? (Options are: grib2, nc)',
                        default='grib2',
                        required=False)

        parser_plot.add_argument('-d',
                        '--directory',
                        type=str,
                        help='directory to read data from',
                        default='.')

        parser_plot.add_argument('--prefix',
                        type=str,
                        help='Set to prefix of files if any',
                        default='',
                        required=False)

        parser_plot.add_argument('--postfix',
                        type=str,
                        help='Set to postfix of files if any',
                        default='',
                        required=False)

        parser_plot.add_argument('-o',
                        '--output-dir',
                        metavar='OUTDIR',
                        type=str,
                        help='Directory to place output into',
                        default='.',
                        required=False)

        parser_plot.add_argument('-l',
                        '--limit-files',
                        metavar='LIMIT',
                        type=int,
                        help='Only use the first LIMIT files. If set to 0, not limit is used. \
                            If Limit > 0, files will be sorted by name first',
                        default=0,
                        required=False)

        parser_plot.add_argument('-a',
                        '--area',
                        metavar='AREA',
                        type=str,
                        help='Over which area to plot (Options are: dk, neu, sjalland, disko)',
                        default="dk",
                        required=False)

        parser_plot.add_argument('--verbose',
                        action='store_true',
                        help='Verbose output',
                        default=False)

        # Parser for EPS Meteogram
        parser_epsmeteogram = subparser.add_parser('epsmeteogram', help='Plot EPS meteogram')

        parser_epsmeteogram.add_argument('-d',
                            '--directory',
                            type=str,
                            help='directory to read data from',
                            default='.')

        parser_epsmeteogram.add_argument('-m',
                            '--members',
                            type=int,
                            help='Number of members',
                            required=True)

        parser_epsmeteogram.add_argument('-f',
                            '--files-per-member',
                            type=int,
                            help='Number of gribfiles per member',
                            required=True)

        parser_epsmeteogram.add_argument('--latlon',
                            type=str,
                            help='Coordinates to plot. Seperate with ":", eg: "55,12:60,14"',
                            required=True)

        parser_epsmeteogram.add_argument('--prefix',
                            type=str,
                            help='Set to prefix of files if any',
                            default='',
                            required=False)

        parser_epsmeteogram.add_argument('--postfix',
                            type=str,
                            help='Set to postfix of files if any',
                            default='',
                            required=False)

        parser_epsmeteogram.add_argument('-o',
                            '--output-dir',
                            metavar='OUTDIR',
                            type=str,
                            help='Directory to place output into',
                            default='.',
                            required=False)

        parser_epsmeteogram.add_argument('--verbose',
                            action='store_true',
                            help='Verbose output',
                            default=False)


        if len(sysargs)==1:
            parent_parser.print_help()
            sys.exit(2)

        args = parent_parser.parse_args()

        return args