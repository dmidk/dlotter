# dlotter
Plot data from NWP output (BETA - in development).\
The "word" *dlotter* is a contraction of *DMI* and *Plotter*.

# Install
Installation should be relatively simple. 

Install using pip: `pip install dlotter`

*One dependency is "cartopy". pip installation of cartopy seems to be currently broken, thus leads to a failed installation of dlotter. If this happens, either build cartopy from source or with conda (see section "Prerequisites").*

# Usage
Submodules can be imported, but this software is designed to be used as a CLI (Command-Line-Interface), by calling the `__main__` module.

`python -m dlotter`: Prints a help message

`python -m dlotter plot`: Calls the plot submodule and prints a help message (as flags are needed for this command to work)

`python -m dlotter plot -p t2m:w10m -d /path/to/grib2/files/ --filetype grib2 --verbose`: Plots t2m and w10m from grib2 files in the directory specified.

`python -m dlotter plot -p t2m:w10m -d /path/to/grib2/files/ --prefix 0 --filetype grib2 --limit-files 4 --verbose`: Limits to use the first 4 files sorted by name and only use files starting with "0" in their name.

## Flags and commands
'-a' *or* '--area' `dk,neu`\
'-d' *or* '--directory' `/path/to/files/`\
'-p' *or* '--parameters' `t2m,w10m,precip`


# Important assumptions
- All files should be in the same projection. This is to speed up the plotting.
- All files will be plotted in the same domain. This is also to speed up plotting.
- All data at one time are in the same file

# Prerequisites
*dlotter* is a script tool based on python. Dependencies is defined in `dlotter.yml` and a python environment can be made using (mini)conda:

`conda env create -f dlotter.yml`


## Build package
To build from the source code:
```sh
python setup.py sdist bdist_wheel
```
To install locally:
```sh
pip install dist/dlotter-$VERSION-py3-none-any.whl --upgrade
```
To deploy to PyPi (Only users with access can do this):
```sh
python -m twine upload dist/*
```