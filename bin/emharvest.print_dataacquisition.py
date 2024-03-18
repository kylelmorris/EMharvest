#!/usr/bin/env python
#
# On local system you can use
#! /usr/bin/env python
# On DLS system you need to use a specific python installation
#! /dls_sw/apps/EM/ebic.scripts/0.0/ebic.scripts_python/bin/python

############################################################################
#
# Author: "Kyle L. Morris"
# eBIC Diamond Light Source 2022
#
############################################################################

# You will need to load dials on the DLS system to get xmltodict
# module load dials
# libtbx.python

# Changed xml parsing method to follow that of DJHatton (DLS)
# https://github.com/d-j-hatton/python-smartem/blob/main/src/smartem/parsing/epu.py

import os
import argparse
parser = argparse.ArgumentParser()
from pathlib import Path
import json
from rich.pretty import pprint
from typing import Any, Dict, Tuple
import xmltodict

# Argument parsing
# https://stackoverflow.com/questions/11604653/how-to-add-command-line-arguments-with-flags-in-python3
parser.add_argument("-i", "--input", help="Required: Input directory")
args = parser.parse_args()
print( "xml {}".format(
        args.input
        ))

def main():
    # Default behaviour for finding timestamps of exposures is fast method
    if not args.input:
        print('Have to define xml input')
        exit(1)

    print_epu_xml(args.input)

def print_epu_xml(xml_path: Path) -> Dict[str, Any]:
    # Use this function for troubleshooting/viewing the raw xml to find data structure
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    data = json.loads(json.dumps(data))
    pprint(data)

main()
