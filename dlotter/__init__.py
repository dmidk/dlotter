"""__init__ of dlotter
"""
import os
import sys
sys.path.insert(0, os.path.abspath('./dlotter/'))
from .prepare import prepare
from .read import grib2Read
from .plot import plot
from .arguments import arguments
from .meteogram import meteogram