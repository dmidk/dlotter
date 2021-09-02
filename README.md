# dlotter
Plot data from NWP output (BETA - in development).\
The "word" *dlotter* is a contraction of *DMI* and *Plotter*.

# Usage
Submodules can be imported, but this software is designed to be used as a CLI (Command-Line-Interface), by calling the `__main__` module.

`python -m dlotter`: Prints a help message

`python -m dlotter plot`: Calls the plot submodule and prints a help message (as flags are needed for this command to work)

# Prerequisites
*dlotter* is a script tool based on python. Dependencies is defined in `dlotter.yml` and a python environment can be made using (mini)conda:

`conda env create -f dlotter.yml`