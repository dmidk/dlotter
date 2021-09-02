#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.prepare
Called from dlotter.__main__
"""
import sys
import argparse
from dmit import ostools

class prepare:

    def __init__(self, args:argparse.Namespace) -> None:

        self.check_meta(args)

        self.files_to_read = self.find_files_to_read(args)
        
        return


    def check_meta(self, args:argparse.Namespace) -> None:
        """Checks if it safe to proceed or user needs to give change inputs

        Parameters
        ----------
        args : argparse.Namespace
            Input arguments to dlotter.__main__
        """

        dir_state = ostools.does_dir_exist(args.directory)
        if not dir_state:
            print("Input directory: {}, does not exist".format(args.directory), flush=True)
            sys.exit(1)

        return


    def find_files_to_read(self, args:argparse.Namespace) -> list:
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

        directory = args.directory
        files = ostools.find_files(directory, 
                                   prefix=args.prefix, 
                                   postfix=args.postfix,
                                   recursive=False, 
                                   onlyfiles=True,
                                   fullpath=True,
                                   olderthan=None,
                                   inorder=False)

        return files