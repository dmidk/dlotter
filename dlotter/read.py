#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for dlotter.read
Called from dlotter.__main__
"""
import sys
import argparse
import pandas as pd
import eccodes as ec

class grib2Read:

    def __init__(self, args:argparse.Namespace, files_to_read:list) -> None:
        if args.verbose:
            print("Reading GRIB2", flush=True)

        data = self.read(args, files_to_read)
        
        return


    def read(self, args:argparse.Namespace, files_to_read:list) -> pd.DataFrame:

        

        return