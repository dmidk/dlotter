#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Master module for setup dlotter
"""
import setuptools
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='dlotter',
     version='0.1.0',
     author="Kasper Hintz",
     author_email="kah@dmi.dk",
     description="Quick and dirty static plots of NWP output",
     long_description=long_description,
     long_description_content_type="text/markdown",
     packages=setuptools.find_packages(),
     setup_requires=[
         'wheel',
         ],
     install_requires=[
         'wheel',
         'dmit',
         'xarray',
         'numpy',
         'cartopy',
         'eccodes',
         'pygrib',
         'scipy',
         ],
     url="https://dmidk.github.io/dlotter/",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
     include_package_data=True,
     package_data={'dlotter': glob('dlotter/icons/*.png')},
 )
