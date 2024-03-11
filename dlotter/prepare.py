#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.prepare
Called from dlotter.__main__
"""
import sys
import argparse
from dmit import ostools

class prepare:
    """Class for preparing input
    """

    def __init__(self, args:argparse.Namespace) -> None:
        """Constructor for prepare class

        Parameters
        ----------
        args : argparse.Namespace
            Input to dlotter from command line
        """

        if args.cmd == 'plot':
            self.valid_parameters = ['t2m', 'w10m', 'precip', 'slp', 'td2m',
                                     'tcc', 'lmhc', 'snow', 'ws', 'cape', 'z']
            self.check_meta(args)
            self.files_to_read = self.find_files_to_read(args)

        if args.cmd == 'plotdiff':
            self.valid_parameters = ['t2m', 'w10m', 'precip', 'td2m', 'z']
            self.check_meta(args)
            self.files_to_read = self.find_files_to_read(args)

        # if args.cmd == 'epsmeteogram':
        #     self.files_to_read = self.find_files_to_read(args, epsmode=True)

        return


    def check_meta(self, args:argparse.Namespace) -> None:
        """Checks if it safe to proceed or user needs to give change inputs

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments to dlotter.__main__
        """

        directory = args.directory.split(',')
        for d in directory:
            dir_state = ostools.does_dir_exist(d)
            if not dir_state:
                print("Input directory: {}, does not exist".format(args.directory), flush=True)
                sys.exit(1)

        allowed_found = False
        parameters = args.parameters.split(':')
        for p in parameters:
            if p in self.valid_parameters:
                allowed_found = True
        if not allowed_found:
            print("No valid parameters found in input argument", flush=True)
            sys.exit(1)

        return


    def find_files_to_read(self, args:argparse.Namespace, epsmode:bool=False) -> list:
        """Finds the file(s) to read and plot from

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments to dlotter.__main__

        Returns
        -------
        list
            List with path(s) to file(s)
        """

        directory = args.directory.split(",")

        if not epsmode:
            if args.limit_files > 0:
                inorder = True
            else:
                inorder = False

        if epsmode:
            recursive = True
            inorder = True
        else:
            recursive = False

        files = []
        for d in directory:
            f = ostools.find_files(d,
                                prefix=args.prefix,
                                postfix=args.postfix,
                                recursive=recursive,
                                onlyfiles=True,
                                fullpath=True,
                                olderthan=None,
                                inorder=inorder)

#       Use this with care
        if args.harmonie_hourly:
            new_files = []
            for f in files:
                filename = f.split('/')[-1]
                if len(filename)<4: new_files.append(f)
            files = new_files


        if not epsmode:
            if args.limit_files > 0:
                if args.limit_files >= len(files):
                    limit = len(files)
                else:
                    limit = args.limit_files

                    f = f[0:limit]

        if epsmode:
            epsfiles = []
            for member in range(args.members):
                for fm in range(args.files_per_member):
                    #Follows: /something/mbr000/000
                    epsfiles.append("{}/mbr{:03d}/{:03d}".format(directory,member,fm))

            if len(epsfiles) != len(f):
                print('Number of automatic found files was different than specified', \
                        flush=True)
            f = epsfiles

        files.append(f)


        return files
