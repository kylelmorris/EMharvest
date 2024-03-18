## EMharvest

# Background

These scripts will parse the xml and mrc files in an EPU session directory and produce a file containing information for deposition.

# Installation

Clone the repository to your local machine

In a terminal, change directory into the cloned repository directory

Run the following command to set up a vritual environment for EMharvest

$ python -m venv python

$ source python/bin/activate

$ pip install .

Each time you want to use EMharvest, remember to activate that virtual environment

$ source python/bin/activate

# Usage instructions

Navigate to the directory you want to send outputs to (optional)

$ cd /path/to/working/output/area

Run a harvest on a single visit diretory

$ emharvest.harvest.py --help

$ emharvest.harvest.py

Use a GUI to achieve the same as above

$ emharvest.py