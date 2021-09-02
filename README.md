# dlotter
Plot data from NWP output (BETA - in development).\
The "word" *dlotter* is a contraction of *DMI* and *Plotter*.

# Usage
Submodules can be imported, but this software is designed to be used as a CLI (Command-Line-Interface), by calling the `__main__` module.

`python -m dlotter`: Prints a help message

`python -m dlotter plot`: Calls the plot submodule and prints a help message (as flags are needed for this command to work)

`python -m dlotter plot -p t2m:w10m -d /path/to/grib2/files/ --filetype grib2 --verbose`: Plots t2m and w10m from grib2 files in the directory specified.

# Important assumptions
- All files should be in the same projection. This is to speed up the plotting.
- All files will be plotted in the same domain. This is also to speed up plotting.
- All data at one time are in the same file

# Prerequisites
*dlotter* is a script tool based on python. Dependencies is defined in `dlotter.yml` and a python environment can be made using (mini)conda:

`conda env create -f dlotter.yml`