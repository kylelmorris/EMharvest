#!/usr/bin/env python3

import os, re
import fnmatch
import math
import argparse
import datetime
import dateutil.parser
import glob
from pathlib import Path
import pandas as pd
import numpy as np

import json
from rich.pretty import pprint
from typing import Any, Dict, Tuple
import xmltodict

import hashlib

import gemmi

from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxWriter import PdbxWriter

prog = "EM HARVEST"
usage = """
        Harvesting microscopy data for automatic depostion.
        Example:
        For single particle data usage is:
        python emharvest.harvest.py -m SPA -c epu -e ../_repo_data/Supervisor_20230919_140141_84_bi23047-106_grid1/EpuSession.dm -a ../_repo_data/atlas/Supervisor_20230919_115905_Atlas_bi23047-106/ScreeningSession.dm  -o ../_repo_data/
        or for non-standard epu files the command is like:
        python emharvest.harvest.py -m SPA -c non-epu -i ../_repo_data/cm33879-3/raw5/metadata/Images-Disc1/GridSquare_30454884/GridSquare_20230531_130838.xml  -o ../_repo_data/dep_cm33879-3
        or for tomogram data usage is:
        python emharvest.harvest.py -m TOMO -t ../_repo_data/tomo_data/SearchMaps/Overview.xml -o ../_repo_data/ -m TOMO -d ../_repo_data/tomo_data/Position_1_33.mdoc
        """

parser = argparse.ArgumentParser(description="Microscopy Data Harvest Script")
parser.add_argument("-m", "--mode", choices=["SPA", "TOMO"], required=True,
                    help="Mode: SPA for Single Particle Analysis or TOMO for Tomography")
parser.add_argument("-c", "--category", choices=["epu", "non-epu"], help="Kind of microscopy input files that needs to be harvested (Required for SPA mode)")
parser.add_argument("-i", "--input_file", help="Any input SPA file which is in xml format (Required for SPA mode)")
parser.add_argument("-e", "--epu", help="EPU session file: Session.dm (Required for SPA mode)")
parser.add_argument("-a", "--atlas", help="Atlas session file: ScreeningSession.dm (Required for SPA mode)")
parser.add_argument("-o", "--output", help="Output directory for generated reports")
parser.add_argument("-p", "--print", help="Optional: Y = Only print xml and exit")
parser.add_argument("-t", "--tomogram_file", help="Tomogram file (Required for TOMO mode)")
parser.add_argument("-d", "--mdoc_file", help="Tomography mdoc file (Required for TOMO mode)")
args = parser.parse_args()


def main():
    if args.mode == "SPA" and args.category == "epu":
        if not args.epu or not args.atlas:
            parser.error("SPA mode requires both --epu and --atlas files.")

        main.epu_xml = args.epu
        main.epu_directory = os.path.dirname(args.epu)

        if args.print:
            print_epu_xml(main.epu_xml)
            exit(1)

        main.atlas_xml = args.atlas
        main.atlas_directory = os.path.dirname(args.atlas)

        output_dir = os.path.join(os.getcwd(), "emharvest")
        main.dep_dir = os.path.join(output_dir, "dep")

        if not os.path.exists(main.dep_dir):
            os.makedirs(main.dep_dir)

        perform_minimal_harvest_epu(main.epu_xml, output_dir)

    if args.mode == "SPA" and args.category == "non-epu":
        main.input_xml = args.input_file
        main.input_directory = os.path.dirname(args.input_file)

        if args.print:
            print_epu_xml(main.input_xml)
            exit(1)

        output_dir = os.path.join(os.getcwd(), "emharvest")
        main.dep_dir = os.path.join(output_dir, "dep")

        if not os.path.exists(main.dep_dir):
            os.makedirs(main.dep_dir)

        perform_spa_harvest_nonepu(main.input_xml, output_dir)

    elif args.mode == "TOMO":
        tomogram_file = args.tomogram_file
        mdoc_file = args.mdoc_file
        if not tomogram_file or not mdoc_file:
            parser.error("TOMO mode requires both --tomogram_file and a --mdoc file.")

        output_dir = os.path.join(os.getcwd(), "emharvest")
        main.dep_dir = os.path.join(output_dir, "dep_tomo")

        if not os.path.exists(main.dep_dir):
            os.makedirs(main.dep_dir)

        perform_tomogram_harvest(tomogram_file, mdoc_file, output_dir)


def perform_tomogram_harvest(tomogram_file, mdoc_file, output_dir):
    print(f"Processing tomogram data from file: {tomogram_file} and {mdoc_file}")
    print(f"Output will be saved to: {output_dir}/dep_tomo")

    FoilDataDict = FoilHoleData(tomogram_file)
    FoilDataDict['tiltAngleMax'] = "?"
    FoilDataDict['tiltAngleMin'] = "?"
    main_sessionName = FoilDataDict["sessionName"]

    EpuDataDict = dict(main_sessionName=main_sessionName, grid_topology="?", grid_material="?",
                       nominal_defocus_min_microns="?", nominal_defocus_max_microns="?",
                       collection="?", number_of_images="?", spot_size="?", C2_micron="?", Objective_micron="?",
                       Beam_diameter_micron="?")

    TomoOverViewDataDict = {**FoilDataDict, **EpuDataDict}

    OverViewDataDict = AnyXMLDataFile(tomogram_file)
    TomoDataDict = {**TomoOverViewDataDict, **OverViewDataDict}

    TomoMdocDataDict = TomoMdocData(mdoc_file)

    TomoDataDict['xmlMag'] = TomoMdocDataDict['Magnification']
    CompleteTomoDataDict = {**TomoDataDict, **TomoMdocDataDict}

    save_deposition_file(CompleteTomoDataDict)

def perform_spa_harvest_nonepu(input_spa_file, output_dir):
    print(f"Processing tomogram data from file: {input_spa_file}")
    print(f"Output will be saved to: {output_dir}/dep_tomo")

    SPADataDict = FoilHoleData(input_spa_file)
    main_sessionName = SPADataDict["sessionName"]

    EpuDataDict = dict(main_sessionName=main_sessionName, grid_topology="?", grid_material="?",
                       nominal_defocus_min_microns="?", nominal_defocus_max_microns="?",
                       collection="?", number_of_images="?", spot_size="?", C2_micron="?", Objective_micron="?",
                       Beam_diameter_micron="?")

    NonEpuDataDict = AnyXMLDataFile(input_spa_file)

    SPATotalDataDict = {**SPADataDict, **EpuDataDict}

    CompleteSPADataDict = {**SPATotalDataDict, **NonEpuDataDict}

    save_deposition_file(CompleteSPADataDict)


def findpattern(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    # result = glob.glob('**/'+str(pattern), recursive=True)
    return result


def roundup(n, decimals=0):
    # https://realpython.com/python-rounding/
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier


def checksum(path, out):
    # https://www.quickprogrammingtips.com/python/how-to-calculate-sha256-hash-of-a-file-in-python.html
    filename = path
    sha256_hash = hashlib.sha256()
    with open(filename, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
        checksum = sha256_hash.hexdigest()

    # Open file for writing
    checkfile = open(out, "w")
    # Write checksum string
    n = checkfile.write(checksum)
    # Close file
    checkfile.close()

    print('Created checksum')
    print()


def print_epu_xml(xml_path: Path) -> Dict[str, Any]:
    # Use this function for troubleshooting/viewing the raw xml to find data structure
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    data = json.loads(json.dumps(data))
    pprint(data)


def searchSupervisorAtlas(path):
    print('Searching Supervisor Atlas directory for xmls, mrc, and jpg')
    print()

    xmlAtlasList = findpattern('Atlas*.xml', path)  # You then need to remove any item containing *Data*
    try:
        xmlAtlas = xmlAtlasList[0]
        # print(xml)
    except:
        xmlAtlas = 'None'

    searchSupervisorAtlas.xmlAtlasList = xmlAtlasList
    searchSupervisorAtlas.xmlAtlas = xmlAtlas

    xmlAtlasTileList = findpattern('Tile*.xml', path)  # You then need to remove any item containing *Data*
    try:
        xmlAtlasTile = xmlAtlasTileList[0]
        # print(xml)
    except:
        xmlAtlasTile = 'None'

    searchSupervisorAtlas.xmlAtlasTileList = xmlAtlasTileList
    searchSupervisorAtlas.xmlAtlasTile = xmlAtlasTile

    print('Found representative xml file for pulling meta data about Atlas session')
    print('Atlas: ' + searchSupervisorAtlas.xmlAtlas)
    print('Atlas tile: ' + searchSupervisorAtlas.xmlAtlasTile)
    print()

    # Store representative xml as global dictionary for reference anywhere in script (reduce I/O)
    with open(xmlAtlas, "r") as xml:
        for_parsing = xml.read()
        searchSupervisorAtlas.xmlAtlasDict = xmltodict.parse(for_parsing)

    with open(xmlAtlasTile, "r") as xml:
        for_parsing = xml.read()
        searchSupervisorAtlas.xmlAtlasTileDict = xmltodict.parse(for_parsing)


def searchSupervisorData(path):
    print('Searching Supervisor Data directory for xmls, mrc, and jpg')
    print()

    print('Finding GridSquare xml')
    xmlSquareList = findpattern('GridSquare*.xml', path)
    try:
        xmlSquare = xmlSquareList[0]
        print('Done')
        # print(xml)
    except:
        xmlSquare = 'None'
        print('None found')

    print('Finding FoilHole xml')
    xmlHoleList = findpattern('FoilHole*.xml', path)
    xmlHoleList = [x for x in xmlHoleList if
                   "Data" not in x]  # This will remove items in list containing *Data*, i.e. DataAcquisition xml files
    try:
        xmlHole = xmlHoleList[0]
        print('Done')
        # print(xml)
    except:
        xmlHole = 'None'
        print('None found')

    print('Finding AcquisitionData xml')
    xmlDataList = findpattern('FoilHole*Data*.xml', path)
    try:
        xmlData = xmlDataList[0]
        print('Done')
        # print(xml)
    except:
        xmlData = 'None'
        print('None found')

    print('Finding AquisitionData mrc')
    mrcDataList = findpattern('FoilHole*Data*.mrc', path)
    try:
        mrc = mrcDataList[0]
        print('Done')
        # print(mrc)
    except:
        mrc = 'None'
        print('None found')

    print('Finding AquisitionData jpg')
    jpgDataList = findpattern('FoilHole*Data*.jp*g', path)
    try:
        jpg = jpgDataList[0]
        print('Done')
        # print(jpg)
    except:
        jpg = 'None'
        print('None found')

    searchSupervisorData.xmlSquareList = xmlSquareList
    searchSupervisorData.xmlHoleList = xmlHoleList
    searchSupervisorData.xmlDataList = xmlDataList

    searchSupervisorData.xmlSquare = xmlSquare
    searchSupervisorData.xmlHole = xmlHole
    searchSupervisorData.xmlData = xmlData

    print('Found representative xml file for pulling meta data about EPU session')
    print('Square: ' + searchSupervisorData.xmlSquare)
    print('Hole: ' + searchSupervisorData.xmlHole)
    print('Acquisition: ' + searchSupervisorData.xmlData)
    print()

    # Store representative xml as global dictionary for reference anywhere in script (reduce I/O)
    try:
        with open(xmlSquare, "r") as xml:
            for_parsing = xml.read()
            searchSupervisorData.xmlSquareDict = xmltodict.parse(for_parsing)
    except:
        print('searchSupervisorData error')

    try:
        with open(xmlHole, "r") as xml:
            for_parsing = xml.read()
            searchSupervisorData.xmlHoleDict = xmltodict.parse(for_parsing)
    except:
        print('searchSupervisorData error')

    try:
        with open(xmlData, "r") as xml:
            for_parsing = xml.read()
            searchSupervisorData.xmlDataDict = xmltodict.parse(for_parsing)
    except:
        print('searchSupervisorData error')

    searchSupervisorData.mrcDataList = mrcDataList

    searchSupervisorData.jpgDataList = jpgDataList


def getXmlMag(xml_path: Path) -> Dict[str, Any]:
    try:
        with open(xml_path, "r") as xml:
            for_parsing = xml.read()
            data = xmltodict.parse(for_parsing)
        data = data["MicroscopeImage"]
    except:
        xmlMag = 0
        xmlAPix = 0
    else:
        xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
        xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]

        xmlAPix = float(xmlMetrePix) * 1e10
        xmlAPix = roundup(xmlAPix, 1)

    return xmlMag, str(xmlAPix)


def xml_presets(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    ## Presets
    # Loop through the presets in the Microscope Settings list
    presets = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
        "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"]
    length = len(presets)
    camera = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
        "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][0]["value"]["b:Acquisition"]["c:camera"][
        "c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"]
    lengthCam = len(camera)

    # Create list for gathering preset conditions for reporting
    # DEV DEV DEV might be more flexible to have these going into dataframe
    namePresetList = []
    probePresetList = []
    magPresetList = []
    apixPresetList = []
    spotPresetList = []
    c2PresetList = []
    beamDPresetList = []
    defocusPresetList = []
    timePresetList = []
    binPresetList = []

    # Loop to gather all microscope presets used for session
    for x in range(0, length):
        name = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["key"]

        # Get magnifications from image xml, they are not stored in the epu session file
        if name == 'Atlas':
            mag = getXmlMag(searchSupervisorAtlas.xmlAtlas)[0]
            apix = getXmlMag(searchSupervisorAtlas.xmlAtlas)[1]
        elif name == 'GridSquare':
            mag = getXmlMag(searchSupervisorData.xmlSquare)[0]
            apix = getXmlMag(searchSupervisorData.xmlSquare)[1]
        elif name == 'Hole':
            mag = getXmlMag(searchSupervisorData.xmlHole)[0]
            apix = getXmlMag(searchSupervisorData.xmlHole)[1]
        elif name == 'Acquisition':
            mag = getXmlMag(searchSupervisorData.xmlData)[0]
            apix = getXmlMag(searchSupervisorData.xmlData)[1]
        else:
            mag = 0
            apix = 0

        probeMode = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:ProbeMode"]
        spot = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:SpotIndex"]
        c2 = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:Apertures"][
            "c:C2Aperture"]["c:Diameter"]
        beamD = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:BeamDiameter"]
        # Deal with two condensor lens systems that don't know beam diameter
        if isinstance(beamD, dict):
            beamD = 0
        else:
            beamD = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
                "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"][
                "c:BeamDiameter"]
        beamDmicron = float(beamD) * 1e6
        DF = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:Defocus"]
        DFmicron = float(DF) * 1e6
        time = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"][
            "c:ExposureTime"]
        epuBin = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"][
            "KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"][
            "c:Binning"]["d:x"]

        # Here we face data in lists and not always in the same list position so need to loop to find position
        # for y in range(0, lengthCam):
        # listKeyValue = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"][y]["key"]["#text"]
        # if listKeyValue == 'SuperResolutionFactor':
        # superRes = listKeyValue
        # superResOn = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"][y]["value"]["#text"]

        # Rounding
        spot = round(float(spot))
        c2 = round(float(c2))
        beamDmicron = roundup(float(beamDmicron), 1)
        DFmicron = roundup(float(DFmicron), 1)
        time = roundup(float(time), 2)
        mag = round(float(mag))
        apix = roundup(float(apix), 3)

        # Report presets or continue silentily
        print(name)
        print('Nominal magnification: ' + str(mag) + ' X')
        print('Pixel size: ' + str(apix) + ' apix')
        print('Probe mode: ' + str(probeMode))
        print('Spot: ' + str(spot))
        print('C2 apeture: ' + str(c2) + ' microns')
        print('Beam diameter: ' + str(beamDmicron) + ' microns')
        print('Defocus: ' + str(DFmicron) + ' microns')
        print('Exposure time: ' + str(time) + ' seconds')
        print('')

        # Append presets to the preset lists for reporting
        namePresetList.append(name)
        magPresetList.append(mag)
        apixPresetList.append(apix)
        probePresetList.append(probeMode)
        spotPresetList.append(spot)
        c2PresetList.append(c2)
        beamDPresetList.append(beamDmicron)
        defocusPresetList.append(DFmicron)
        timePresetList.append(time)
        binPresetList.append(epuBin)

        # Gather main params for reporting
        if name == 'Acquisition':
            xml_presets.time = time
            xml_presets.beamD = beamDmicron
            xml_presets.probe = probeMode
            xml_presets.C2 = c2
            xml_presets.spot = spot
            xml_presets.epuBin = epuBin
            xml_presets.mag = mag
            xml_presets.apix = apix
            # xml_presets.superRes = superRes
            # xml_presets.superResOn = superResOn
        if name == 'AutoFocus':
            xml_presets.beamDAutoFocus = beamDmicron

        # Gather all presets for mass reporting
        xml_presets.namePresetList = namePresetList
        xml_presets.magPresetList = magPresetList
        xml_presets.apixPresetList = apixPresetList
        xml_presets.probePresetList = probePresetList
        xml_presets.spotPresetList = spotPresetList
        xml_presets.c2PresetList = c2PresetList
        xml_presets.beamDPresetList = beamDPresetList
        xml_presets.defocusPresetList = defocusPresetList
        xml_presets.timePresetList = timePresetList
        xml_presets.binPresetList = binPresetList

    # report complete
    print('Finished gathering all microscope presets')


def xml_presets_data(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    ## SuperResolutionBinning Factor
    # The SuperResolutionFactor is not always in the same list position in a:KeyValueOfstringanyType
    keyValueList = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"]

    # Loop through the list to find the SuperResolutionFactor list position
    i = 0
    for value in keyValueList:
        key = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i][
            "a:Key"]
        if key == "SuperResolutionFactor":
            j = i
        i = i + 1

    try:
        superResBin = \
        data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][j][
            "a:Value"]["#text"]
    except:
        superResBin = 'Unknown'

    ## Energy filter
    # Known error in nt29493-49 - glacios
    try:
        filterSlit = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitInserted"]
    except:
        filterSlit = 'None'
    else:
        filterSlit = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitInserted"]
    try:
        filterSlitWidth = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitWidth"]
    except:
        filterSlitWidth = 'None'
    else:
        filterSlitWidth = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitWidth"]

    # Aperture(s)
    # Loop through the list to find the Objective aperture list position
    i = 0
    for value in keyValueList:
        key = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Key"]
        if key == "Aperture[OBJ].Name":
            keyvalue = i
            objectiveAperture = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Value"]["#text"]
        i = i + 1

    # Stage tilt
    stageAlphaRad = getStageTilt(micpath)[0]
    stageBetaRad = getStageTilt(micpath)[1]
    stageAlpha = roundup(math.degrees(float(stageAlphaRad)), 2)
    stageBeta = roundup(math.degrees(float(stageBetaRad)), 2)

    # Report
    xml_presets_data.superResBin = superResBin
    xml_presets_data.filterSlitWidth = filterSlitWidth
    xml_presets_data.filterSlit = filterSlit
    xml_presets_data.stageAlpha = roundup(float(stageAlpha), 1)
    xml_presets_data.stageBeta = roundup(float(stageBeta), 1)
    xml_presets_data.objective = objectiveAperture


def getStageTilt(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # Find the stage Alpha (DEV DEV think the units might be 1/100th)
    stageAlpha = data["microscopeData"]["stage"]["Position"]["A"]
    stageBeta = data["microscopeData"]["stage"]["Position"]["B"]

    return [stageAlpha, stageBeta]


def xml_session(xml_path: Path) -> pd.DataFrame:
    data_dict = {}
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    # Location of EPU session directory on which this script was ran
    data_dict['realPath'] = os.path.realpath(xml_path)

    # EPU version
    epuId = data["Version"]["@z:Id"]
    epuBuild = data["Version"]["a:_Build"]
    epuMajor = data["Version"]["a:_Major"]
    epuMinor = data["Version"]["a:_Minor"]
    epuRevision = data["Version"]["a:_Revision"]
    data_dict['epuVersion'] = str(epuMajor) + '.' + str(epuMinor) + '.' + str(epuRevision) + '-' + str(
        epuId) + '.' + str(epuBuild)

    # Output format
    data_dict['doseFractionOutputFormat'] = data["DoseFractionsOutputFormat"]["#text"]

    # Autoloader slot (starts at 0)
    autoSlot = data["AutoloaderSlot"]
    data_dict['autoSlot'] = float(autoSlot) + 1

    # SampleXml is a list denoted inside [], pull 0 out
    # Then the text of AtlasId is in a key pair in the dictionary denoted by {}, so call that key pair
    # Use the print_epu_xml def to see full [] and {} formatted data structure
    data_dict['atlasDir'] = data["Samples"]["_items"]["SampleXml"][0]["AtlasId"]["#text"]

    # Session name and creation time
    sessionName = xml_sessionName(xml_path)
    sessionDate = data["Samples"]["_items"]["SampleXml"][0]["StartDateTime"]
    sessionDateFormat = formatEPUDate(sessionDate)
    data_dict['sessionName'] = sessionName
    data_dict['sessionDate'] = sessionDateFormat

    # Grid type - lacey or holeycarbon or holeygold
    try:
        data_dict['gridType'] = data["Samples"]["_items"]["SampleXml"][0]["GridType"]
    except:
        data_dict['gridType'] = 'Unknown'

    # The I0 filter settings may hint at what grid type is being used
    try:
        data_dict['I0set'] = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["IsCalibrated"]
    except:
        data_dict['I0set'] = 'Unknown'

    try:
        data_dict['I0MaxInt'] = round(
            float(data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["MaximumIntensity"]))
    except:
        data_dict['I0MaxInt'] = 'Unknown'

    try:
        data_dict['I0MinInt'] = round(
            float(data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["MinimumIntensity"]))
    except:
        data_dict['I0MinInt'] = 'Unknown'

    # Clustering method
    data_dict['clustering'] = data["ClusteringMode"]
    data_dict['clusteringRadius'] = float(data["ClusteringRadius"]) * 1e6 if data[
                                                                                 "ClusteringMode"] == 'ClusteringWithImageBeamShift' else np.nan

    data_dict['focusWith'] = \
    data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["FocusWith"]["#text"]
    data_dict['focusRecurrence'] = \
    data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["Recurrence"]["#text"]

    data_dict['delayImageShift'] = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"][
        "DelayAfterImageShift"]
    data_dict['delayStageShift'] = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"][
        "DelayAfterStageShift"]

    # AFIS or accurate
    if data["ClusteringMode"] == 'ClusteringWithImageBeamShift':
        data_dict['afisMode'] = 'AFIS'
        data_dict['afisRadius'] = data_dict['clusteringRadius']
    else:
        data_dict['afisMode'] = 'Accrt'
        data_dict['afisRadius'] = np.nan

    # Send xml dict over to function to get defocus range
    defocusRange = getDefocusRange(data)

    # In some cases the defocus list is a single value, test to deal with
    # Also find max and min defocus values
    if isinstance(defocusRange, str):
        # data_dict['defocusList'] = "xml read error"
        data_dict['defocusMax'] = "xml read error"
        data_dict['defocusMin'] = "xml read error"
    elif isinstance(defocusRange, list):
        defocusRangeMicron = [float(element) * 1 for element in defocusRange]
        defocusRangeRound = [round(num, 2) for num in defocusRangeMicron]
        # data_dict['defocusList'] = defocusRangeRound
        data_dict['defocusMax'] = min(defocusRangeRound)
        data_dict['defocusMin'] = max(defocusRangeRound)

    # print(data_dict)
    df = pd.DataFrame(data_dict, index=[0])

    # Print to terminal
    print('EPU version:', data_dict['epuVersion'])
    print('Dose fraction output:', data_dict['doseFractionOutputFormat'])
    print('EPU Session:', data_dict['sessionName'])
    print('EPU Session date:', data_dict['sessionDate'])
    print('EPU Session path:', data_dict['realPath'])
    print()
    print('Clustering mode:', data_dict['clustering'])
    print('Clustering radius:', data_dict['clusteringRadius'])
    print('Focus mode:', data_dict['focusWith'])
    print('Focus recurrance:', data_dict['focusRecurrence'])
    print()
    print('Autoloader slot:', data_dict['autoSlot'])
    print('Atlas directory:', data_dict['atlasDir'])
    print()
    # print('Defocus list:', data_dict['defocus'])
    print('Defocus max:', data_dict['defocusMax'])
    print('Defocus min:', data_dict['defocusMin'])
    print()
    print('Image shift delay:', data_dict['delayImageShift'])
    print('Stage shift delay:', data_dict['delayStageShift'])
    print('Grid type:', data_dict['gridType'])
    print('I0 set:', data_dict['I0set'])
    print('I0 max:', data_dict['I0MaxInt'])
    print('I0 min:', data_dict['I0MinInt'])
    print()
    print('\033[1m' + 'Finished gathering metadata from main EPU session file' + '\033[0m')
    print()

    return df


def xml_sessionName(xml_path):
    # It is necessary to have a function for getting xml session name elsewhere in script
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    sessionName = data["Name"]["#text"]

    return sessionName


def formatEPUDate(d):
    # Returns formatted in datetime, needs to be string for printing

    # https://stackoverflow.com/questions/17594298/date-time-formats-in-python
    # https://www.w3schools.com/python/python_datetime.asp
    # https://www.tutorialexample.com/python-detect-datetime-string-format-and-convert-to-different-string-format-python-datetime-tutorial/amp/
    # Read in EPU formatted date and time - remember input is a string
    epuDate = dateutil.parser.parse(d)
    epuDate = epuDate.strftime("%Y-%m-%d %H:%M:%S")
    # new_date = datetime.datetime.strptime(d,"%Y-%m-%dT%H:%M:%S.%fZ")
    # Reformat date into YY-MM-DD HH-MM-SS
    # return new_date.strftime("%Y-%m-%d %H:%M:%S")
    return datetime.datetime.strptime(epuDate, "%Y-%m-%d %H:%M:%S")


def getDefocusRange(data):
    # Need to deal with cases where it's multi shot or single shot
    templates = \
    data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"][
        "b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]

    if isinstance(templates, list):
        shotType = 'Multishot'
        # DEV DEV This defocus code might need a revisit if you get lots of errors
        d = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"][
            "a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][0]["b:value"][
            "ImageAcquisitionSettingXml"]["Defocus"]
        if d.get("a:double"):
            try:
                df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"][
                    "a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][0]["b:value"][
                    "ImageAcquisitionSettingXml"]["Defocus"]["a:double"]
                # Sometimes the values contain unicode en-dash and not ASCII hyphen
                # df.replace('\U00002013', '-')
            except:
                print('Warning, could not find defocus range in xml file')
                df = ['xml read error']
        else:
            try:
                df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"][
                    "a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"][
                    "ImageAcquisitionSettingXml"]["Defocus"]["a:_items"]["a:double"]
                # Sometimes the values contain unicode en-dash and not ASCII hyphen
            # df.replace('\U00002013', '-')
            except:
                print('Warning, could not find defocus range in xml file')
                df = ['xml read error']
    else:
        shotType = 'Single'
        # There is sometimes a data structure change I think in single shot acqusition, cause currently unknown, check for it
        d = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"][
            "a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"][
            "ImageAcquisitionSettingXml"]["Defocus"]
        if d.get("a:double"):
            try:
                df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"][
                    "a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"][
                    "ImageAcquisitionSettingXml"]["Defocus"]["a:double"]
                # Sometimes the values contain unicode en-dash and not ASCII hyphen
                # df.replace('\U00002013', '-')
            except:
                print('Warning, could not find defocus range in xml file')
                df = ['xml read error']
        else:
            try:
                df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"][
                    "a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"][
                    "ImageAcquisitionSettingXml"]["Defocus"]["a:_items"]["a:double"]
                # Sometimes the values contain unicode en-dash and not ASCII hyphen
                # df.replace('\U00002013', '-')
            except:
                print('Warning, could not find defocus range in xml file')
                df = ['xml read error']

    getDefocusRange.shotType = shotType

    # Remember df in this case means defocus, not dataframe!!
    # Sometimes there is a single value in the defocus list and then this gets stored as a single string
    if isinstance(df, str):
        # Convert to list is str and thus single value
        df = df.split(sep=None, maxsplit=-1)

    # Check for error, which is stored as single item list
    read = df[0]
    if df[0] == "xml read error":
        return df
    # Otherwise convert metres into microns
    else:
        dfMicron = [float(item) * 1e6 for item in df]
        return dfMicron


def AnyXMLDataFile(xmlpath: Path) -> Dict[str, Any]:
    """Reading the Overview.xml file information and storing in a dictionary."""
    with open(xmlpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    acqusition_date = data["microscopeData"]["acquisition"]["acquisitionDateTime"]
    date = acqusition_date.split("T", 1)[0]
    model = data["microscopeData"]["instrument"]["InstrumentModel"]
    microscope_mode = data["microscopeData"]["optics"]["ColumnOperatingTemSubMode"]
    eV = data["microscopeData"]["gun"]["AccelerationVoltage"]
    xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]
    xmlAPix = float(xmlMetrePix) * 1e10
    xmlAPix = roundup(xmlAPix, 1)
    soft_name = data["microscopeData"]["core"]["ApplicationSoftware"]
    if soft_name == "Tomography":
        software_name = "FEI tomography"
    else:
        software_name = soft_name
    software_version = data["microscopeData"]["core"]["ApplicationSoftwareVersion"]
    illumination = data["microscopeData"]["optics"]["IlluminationMode"]

    objectiveAperture, C2_micron = "", ""
    keyValueList = data["CustomData"]["a:KeyValueOfstringanyType"]
    for i, value in enumerate(keyValueList):
        key = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Key"]
        if key == "Aperture[OBJ].Name":
            objectiveAperture = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Value"]["#text"]
        if key == "Aperture[C2].Name":
            C2_micron = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Value"]["#text"]

    OverViewDataDict = dict(date=date, model=model, microscope_mode=microscope_mode, eV=eV, xmlMag=xmlMag,
                            xmlMetrePix=xmlMetrePix, xmlAPix=xmlAPix, objectiveAperture=objectiveAperture,
                            C2_micron=C2_micron, software_name=software_name, software_version=software_version,
                            illumination=illumination)

    return OverViewDataDict


def unique_values(existing_list, new_values):
    """Add unique values to the list, handling NaN separately."""
    for value in new_values:
        if isinstance(value, float) and math.isnan(value):
            if not any(isinstance(v, float) and math.isnan(v) for v in existing_list):
                existing_list.append(value)
        elif value not in existing_list:
            existing_list.append(value)
    return existing_list


def TomoMdocData(mdocpath: Path) -> Dict[str, Any]:
    """Reading the mdoc file information and storing in a dictionary."""
    mdoc_data = {}

    with open(mdocpath, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                line = line[1:-1].strip()
            if not line:
                continue

            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()

                # Handle multiple values in one line
                if " " in value:
                    values = value.split()
                    try:
                        values = [float(v) for v in values]
                    except ValueError:
                        pass
                else:
                    try:
                        values = [float(value)]
                    except ValueError:
                        values = [value]

                # Update the mdoc_data dictionary
                if key in mdoc_data:
                    mdoc_data[key] = unique_values(mdoc_data[key], values)
                else:
                    mdoc_data[key] = values

    for key, value in mdoc_data.items():
        if isinstance(value, list):
            if len(value) == 1:
                mdoc_data[key] = value[0]
            else:
                unique_val = list(set(value))
                if len(unique_val) == 1:
                    mdoc_data[key] = unique_val[0]

    return mdoc_data


def FoilHoleData(xmlpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(xmlpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    sessionName = data["uniqueID"]

    # The values are not always in the same list position in a:KeyValueOfstringanyType
    keyValueList = data["CustomData"]["a:KeyValueOfstringanyType"]

    keyMicroscopeList = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"][
        "a:KeyValueOfstringanyType"]

    # Loop through the list to find the DoseRate list position
    keyvalue = 0
    detectorName, detectorMode, counting, superResolution, objectiveAperture = "", "", "", "", ""
    detector_keys = [
        "Detectors[EF-CCD].CommercialName",
        "Detectors[EF-Falcon].CommercialName"
    ]

    for i, value in enumerate(keyValueList):
        key = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Key"]

        if key == "Detectors[BM-Falcon].DoseRate" or key == "Detectors[EF-Falcon].DoseRate":
            keyvalue = i

        if key in detector_keys:
            detectorName = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Value"]["#text"]
            if detectorName == "BioQuantum K3":
                detectorName = "GATAN K3 BIOQUANTUM (6k x 4k)"
            elif detectorName == "Falcon 4i":
                detectorName = "TFS FALCON 4i (4k x 4k)"

        if key == "Aperture[OBJ].Name":
            objectiveAperture = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Value"]["#text"]

    for i, value in enumerate(keyMicroscopeList):
        keyMicroscopeData = \
        data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i]["a:Key"]
        if keyMicroscopeData == "ElectronCountingEnabled":
            counting = \
            data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i][
                "a:Value"]["#text"]

        if keyMicroscopeData == "SuperResolutionFactor":
            superResolution = \
            data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i][
                "a:Value"]["#text"]

    if counting == "true":
        if superResolution == "1":
            detectorMode = "SUPER-RESOULTION"
        elif superResolution == "2":
            detectorMode = "COUNTING"

    # Retrieve the values
    # xmlDoseRate = data["CustomData"]["a:KeyValueOfstringanyType"][keyvalue]["a:Value"]["#text"]
    xmlDoseRate = "?"  # the data file has only electron_dose on camera and not the dose used on the specimen
    avgExposureTime = data["microscopeData"]["acquisition"]["camera"]["ExposureTime"]
    slitWid = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitWidth"]
    slitInserted = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitInserted"]
    if slitInserted == "true":
        slitWidth = slitWid
    else:
        slitWidth = "?"
    electronSource = data["microscopeData"]["gun"]["Sourcetype"]
    tiltAngleMin = round(float(data["microscopeData"]["stage"]["Position"]["A"]) * (180 / math.pi), 5)
    tiltAngleMax = round(float(data["microscopeData"]["stage"]["Position"]["B"]) * (180 / math.pi), 5)

    FoilHoleDataDict = dict(sessionName=sessionName, xmlDoseRate=xmlDoseRate, detectorName=detectorName,
                            avgExposureTime=avgExposureTime, detectorMode=detectorMode, slitWidth=slitWidth,
                            electronSource=electronSource, tiltAngleMin=tiltAngleMin, tiltAngleMax=tiltAngleMax)

    return FoilHoleDataDict


def find_mics(path, search):
    # Need to have an independent function to find the mics, then move into search_mics to sort them out
    # So find mics can be used independently

    print('Looking for micrograph data in EPU directory using extension: ' + search)
    print('')
    # Just get file names
    # files = glob.glob("./Images-Disc1/GridSquare*/Data/*xml")
    # files.sort(key=os.path.getmtime)
    # print("\n".join(files))

    # Old method of finding xml files
    searchedFiles = glob.glob(path + "/**/GridSquare*/Data/*" + search + '*')
    # searchedFiles = glob.iglob(main.epu+"/Images-Disc1/GridSquare*/Data/*"+search)
    if searchedFiles:
        print('Found micrograph data: ' + str(len(searchedFiles)))
    else:
        print('No micrographs found with search term: ' + search)
        searchedFiles = 'exit'

    # New method of finding xml files
    # searchedFiles = searchSupervisorData.xmlList

    return searchedFiles


def deposition_file(xml):
    # Get EPU session name from main EPU xml file, this is a function
    main_sessionName = xml_sessionName(xml)

    # This is the data xml metadata file already in a dictionary
    data = searchSupervisorData.xmlDataDict["MicroscopeImage"]
    software_version = df_lookup(main.masterdf, 'epuVersion')
    date = df_lookup(main.masterdf, 'sessionDate').strftime("%Y-%m-%d %H:%M:%S")
    nominal_defocus_min_microns = df_lookup(main.masterdf, 'defocusMin')
    nominal_defocus_max_microns = df_lookup(main.masterdf, 'defocusMax')
    collection = df_lookup(main.masterdf, 'afisMode')
    number_of_images = main.mic_count
    spot_size = xml_presets.spot
    C2_micron = xml_presets.C2
    Objective_micron = str(xml_presets_data.objective)
    Beam_diameter_micron = xml_presets.beamD

    # Get mag
    xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]
    xmlAPix = float(xmlMetrePix) * 1e10
    xmlAPix = roundup(xmlAPix, 1)

    # Get scope and kV
    model = data["microscopeData"]["instrument"]["InstrumentModel"]
    eV = data["microscopeData"]["gun"]["AccelerationVoltage"]

    microscope_mode = data["microscopeData"]["optics"]["ColumnOperatingTemSubMode"]
    # illumination = data["microscopeData"]["optics"]["IlluminationMode"]

    grid_type = df_lookup(main.masterdf, 'gridType')
    grid_parts = re.findall(r'[A-Z][a-z]*', grid_type)

    # Now, parts will be ['Holey', 'Carbon']
    grid_topology = grid_parts[0]
    grid_material = grid_parts[1]

    EpuDataDict = dict(main_sessionName=main_sessionName, xmlMag=xmlMag, xmlMetrePix=xmlMetrePix, xmlAPix=xmlAPix,
                       model=model, eV=eV, microscope_mode=microscope_mode, grid_topology=grid_topology,
                       grid_material=grid_material,
                       software_name="EPU", software_version=software_version, date=date,
                       nominal_defocus_min_microns=nominal_defocus_min_microns,
                       nominal_defocus_max_microns=nominal_defocus_max_microns,
                       collection=collection, number_of_images=number_of_images, spot_size=spot_size,
                       C2_micron=C2_micron, Objective_micron=Objective_micron,
                       Beam_diameter_micron=Beam_diameter_micron, illumination="?",
                       PixelSpacing="?", SubFramePath="?")

    FoilHoleDataDict = FoilHoleData(searchSupervisorData.xmlData)
    CompleteDataDict = {**EpuDataDict, **FoilHoleDataDict}
    save_deposition_file(CompleteDataDict)


def SubFramePath(CompleteDataDict, n):
    """ To get the information from SubFramPath value in mdoc file"""
    SubFramePath = CompleteDataDict['SubFramePath'][n]
    file_name = SubFramePath.split('\\')[-1]
    parts = file_name.split('.')[0]
    angle = parts.split('_')[-1]
    return angle


def save_deposition_file(CompleteDataDict):
    # Save doppio deposition csv file
    dictHorizontal1 = {
        'Microscope': CompleteDataDict['model'],
        'software_version': CompleteDataDict['software_version'],
        'date': CompleteDataDict['date'],
        'eV': CompleteDataDict['eV'],
        'mag': CompleteDataDict['xmlMag'],
        'apix': CompleteDataDict['xmlAPix'],
        'nominal_defocus_min_microns': CompleteDataDict['nominal_defocus_min_microns'],
        'grid_topology': CompleteDataDict['grid_topology'],
        'grid_material': CompleteDataDict['grid_material'],
        'nominal_defocus_max_microns': CompleteDataDict['nominal_defocus_max_microns'],
        'spot_size': CompleteDataDict['spot_size'],
        'C2_micron': CompleteDataDict['C2_micron'],
        'Objective_micron': CompleteDataDict['Objective_micron'],
        'Beam_diameter_micron': CompleteDataDict['Beam_diameter_micron'],
        'collection': CompleteDataDict['collection'],
        'number_of_images': CompleteDataDict['number_of_images'],
        'software_name': CompleteDataDict["software_name"],
        'software_category': "IMAGE ACQUISITION",
        'microscope_mode': CompleteDataDict['microscope_mode'],
        'detector_name': CompleteDataDict['detectorName'],
        'dose_rate': CompleteDataDict['xmlDoseRate'],
        'avg_exposure_time': CompleteDataDict['avgExposureTime'],
        'detector_mode': CompleteDataDict['detectorMode'],
        'illumination_mode': CompleteDataDict['illumination'].upper(),
        'slit_width': CompleteDataDict['slitWidth'],
        'electron_source': CompleteDataDict['electronSource'],
        'tilt_angle_min': CompleteDataDict['tiltAngleMin'],
        'tilt_angle_max': CompleteDataDict['tiltAngleMax']
    }
    if args.mode == "TOMO":
        dictHorizontal1.update({'pixel_spacing_x': CompleteDataDict['PixelSpacing'],
                                'pixel_spacing_y': CompleteDataDict['PixelSpacing'],
                                'pixel_spacing_z': CompleteDataDict['PixelSpacing'],
                                'angle_increment': float(SubFramePath(CompleteDataDict, 2)) - float(
                                    SubFramePath(CompleteDataDict, 1)),
                                'rotation_axis': CompleteDataDict['RotationAngle'],
                                'max_angle': SubFramePath(CompleteDataDict, -1),
                                'min_angle': SubFramePath(CompleteDataDict, -2),
                                'angle2_increment': '?',
                                'max_angle2': '?',
                                'min_angle2': '?'
                                })

    df1 = pd.DataFrame([dictHorizontal1])

    # Sample data for the second row
    dictHorizontal2 = {
        'Microscope': 'em_imaging.microscope_model',
        'software_name': 'em_software.name',
        'software_version': 'em_software.version',
        'software_category': 'em_software.category',
        'date': 'em_imaging.date',
        'eV': 'em_imaging.accelerating_voltage',
        'mag': 'em_imaging.nominal_magnification',
        'apix': '?',
        'nominal_defocus_min_microns': 'em_imaging.nominal_defocus_min',
        'nominal_defocus_max_microns': 'em_imaging.nominal_defocus_max',
        'spot_size': '?',
        'C2_micron': 'em_imaging.c2_aperture_diameter',
        'Objective_micron': '?',
        'Beam_diameter_micron': '?',
        'collection': '?',
        'number_of_images': '?',
        "microscope_mode": 'em_imaging.mode',
        "grid_material": 'em_support_film.material',
        "grid_topology": 'em_support_film.topology',
        "detector_name": "em_image_recording.film_or_detector_model",
        "dose_rate": "em_image_recording.avg_electron_dose_per_image",
        "avg_exposure_time": "em_image_recording.average_exposure_time",
        "detector_mode": "em_image_recording.detector_mode",
        "illumination_mode": "em_imaging.illumination_mode",
        "slit_width": "em_imaging_optics.energyfilter_slit_width",
        "electron_source": "em_imaging.electron_source",
        "tilt_angle_min": "em_imaging.tilt_angle_min",
        "tilt_angle_max": "em_imaging.tilt_angle_max"
    }
    if args.mode == "TOMO":
        dictHorizontal2.update({"pixel_spacing_x": "em_map.pixel_spacing_x",
                                "pixel_spacing_y": "em_map.pixel_spacing_y",
                                "pixel_spacing_z": "em_map.pixel_spacing_z",
                                "angle_increment": "em_tomography.axis1_angle_increment",
                                "rotation_axis": "em_tomography.dual_tilt_axis_rotation",
                                "max_angle": "em_tomography.axis1_max_angle",
                                "min_angle": "em_tomography.axis1_min_angle",
                                "angle2_increment": "em_tomography.axis2_angle_increment",
                                "max_angle2": "em_tomography.axis2_max_angle",
                                "min_angle2": "em_tomography.axis2_min_angle"})

    df2 = pd.DataFrame([dictHorizontal2])

    # Append the second row to the DataFrame
    df = pd.concat([df1, df2], ignore_index=True, sort=False)
    # df = df1.merge(df2, left_index=0, right_index=0)

    ## Deposition file
    depfilepath = main.dep_dir + '/' + CompleteDataDict['main_sessionName'] + '_dep.json'
    checksumpath = main.dep_dir + '/' + CompleteDataDict['main_sessionName'] + '_dep.checksum'

    # Human readable deposition file
    # df.to_csv (main.dep_dir+'/'+sessionName+'.dep', index = False, header=True)
    # Manual headings
    headings = ['Value', 'JSON', 'mmCIF']
    # Duplicate the last row
    mmcif_name = df.iloc[-1]
    # Append to DF
    df = df.append(mmcif_name, ignore_index=True)

    tfs_xml_path_list = [
        '[MicroscopeImage][microscopeData][instruments][InstrumentModel]',
        '[MicroscopeImage][microscopeData][core][ApplicationSoftware]',
        '[MicroscopeImage][microscopeData][acquisitionDateTime]',
        '[MicroscopeImage][microscopeData][gun][AccelerationVoltage]',
        '[MicroscopeImage][microscopeData][optics][TemMagnification][NominalMagnification]',
        '[MicroscopeImage][SpatialScale][pixelSize][x][numericValue]',
        '[MicroscopeImage][microscopeData][optics][TemMagnification][NominalMagnification]',
        '[Samples][_items][SampleXml][0][GridType]',
        '[Samples][_items][SampleXml][0][GridType]',
        '[MicroscopeImage][microscopeData][optics][Defocus]',
        '[MicroscopeImage][microscopeData][optics][SpotIndex]',
        '[MicroscopyImage][CustomData][a:KeyValueOfstringanyType][a:Key] is Aperture[C2].Name then extract [<a:Value>]',
        '[MicroscopyImage][CustomData][a:KeyValueOfstringanyType][a:Key] is Aperture[OBJ].Name then extract [<a:Value>]',
        '[MicroscopeImage][microscopeData][optics][BeamDiameter]',
        '?',
        '?',
        '[MicroscopeImage][microscopeData][core][ApplicationSoftwareVersion]',
        'PREDEFINED VALUE',
        '[MicroscopeImage][microscopeData][optics][ColumnOperatingTemSubMode]',
        '[MicroscopyImage][CustomData][a:KeyValueOfstringanyType][a:Key] isDetectorCommercialName then extract [<a:Value>]',
        '?',
        '[MicroscopeImage][microscopeData][acquisition][camera][ExposureTime]',
        '[MicroscopeImage][microscopeData][acquisition][camera][CameraSpecificInput][a:KeyValueOfstringanyType][a:Key] is ElectronCountingEnabled and [<a:Vallue>] is true then COUNTING',
        '[MicroscopeImage][microscopeData][optics][IlluminationMode]',
        '[MicroscopeImage][microscopeData][optics][EnergyFilter][EnergySelectionSlitWidth]',
        '[MicroscopeImage][microscopeData][gun][Sourcetype]',
        '[MicroscopeImage][microscopeData][stage][Position][A]',
        '[MicroscopeImage][microscopeData][stage][Position][B]'
    ]
    if args.mode == "TOMO":
        tfs_xml_path_list.extend(['[PixelSpacing]',
                                  '[PixelSpacing]',
                                  '[PixelSpacing]',
                                  '[SubFramPath]',
                                  '[RotationAngle]',
                                  '[SubFramePath]',
                                  '[SubFramePath]',
                                  '[CryoTomo is usually single axis tilt]',
                                  '[CryoTomo is usually single axis tilt]',
                                  '[CryoTomo is usually single axis tilt]'])

    emdb_xml_path_list = [
        '[emd][structure_determination_list][structure_determination][microscopy_list]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][software_list][software][version]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][date]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][acceleration_voltage]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][nominal_magnification]',
        '?',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][nominal_defocus_min]',
        '[emd][structure_determination_list][structure_determination][specimen_preparation_list][single_particle_preparation][grid][support_film][film_topolgy]',
        '[emd][structure_determination_list][structure_determination][specimen_preparation_list][single_particle_preparation][grid][support_film][film_material]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][nominal_defocus_max]',
        '?',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][c2_aperture_diameter]',
        '?',
        '?',
        '?',
        '?',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][software_list][software][name]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][imaging_mode]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][image_recording_list][image_recording][film_or_detector_model]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][image_recording_list][image_recording][average_electron_dose_per_image]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][image_recording_list][image_recording][average_exposure_time]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][image_recording_list][image_recording][detector_mode]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][microscopy][illumination_mode]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][specialist_optics][energyfilter][slith_width]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][electron_source]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][tilt_angle_min]',
        '[emd][structure_determination_list][structure_determination][microscopy_list][single_particle_microscopy][tilt_angle_max]'
    ]
    if args.mode == "TOMO":
        emdb_xml_path_list.extend(['[emd][map][pixel_spacing][x]',
                                   '[emd][map][pixel_spacing][y]',
                                   '[emd][map][pixel_spacing][z]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis1][angle_increment]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis_rotation]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis1][max_angle]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis1][min_angle]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis2][angle_increment]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis2][max_angle]',
                                   '[emd][structure_determination_list][structure_determination][microscopy_list][tomgraphy_microscopy][tilt_series][axis2][min_angle]'
                                   ])

    # Transpose
    df_transpose = df.T
    # Set headings
    df_transpose.columns = headings
    # Add XML paths
    if len(tfs_xml_path_list) != len(df_transpose):
        raise ValueError("Length of tfs_xml_path_list must match the number of rows in the DataFrame")
    df_transpose['TFS XML Path'] = tfs_xml_path_list
    if len(emdb_xml_path_list) != len(df_transpose):
        raise ValueError("Length of emdb_xml_path_list must match the number of rows in the DataFrame")
    df_transpose['EMDB XML Path'] = emdb_xml_path_list
    # Set the name of the index
    df_transpose.index.name = 'Items'
    # add square brackets around keys
    df_transpose['JSON'] = df_transpose['JSON'].apply(lambda x: '[' + x.replace('.', '][') + ']')
    # Save to CSV
    df_transpose.to_csv(main.dep_dir + '/' + CompleteDataDict['main_sessionName'] + '_dep.csv', index=True, header=True)

    df1_selected = df1.applymap(lambda x: x.item() if isinstance(x, (np.generic, np.ndarray)) else x)

    # Create nested dictionary from the DataFrame
    nested_dict = {}
    for col, new_col in dictHorizontal2.items():
        if '?' not in new_col and col in df1_selected.columns:
            top_key, sub_key = new_col.split('.')
            if top_key not in nested_dict:
                nested_dict[top_key] = {}
            nested_dict[top_key][sub_key] = df1.at[0, col]

    # Convert nested dictionary to JSON
    json_output = json.dumps(nested_dict, indent=4, default=str)

    # write to json file
    with open(depfilepath, 'w') as f:
        f.write(json_output)

    # This can be run before doing full analysis of the session directories
    print("Created deposition file")

    # Deposition Checksum
    checksum(depfilepath, checksumpath)

    # Input data as cif dictionary
    cif_dict = {}

    for key in dictHorizontal2:
        cif_dict[dictHorizontal2[key]] = dictHorizontal1[key]

    # transalating and writting to cif file
    print("CIF_DICTIONARY", cif_dict)
    translate_xml_to_cif(cif_dict, CompleteDataDict['main_sessionName'])


def df_lookup(df, column):
    """
    Perform a lookup in a DataFrame based on a given column and value.

    Parameters:
        df (pd.DataFrame): The DataFrame to search.
        column (str): The name of the column to search.
        value: The value to look for in the specified column.

    Returns:
        pd.DataFrame: Subset of the original DataFrame containing rows where the specified column matches the given value.
    """
    return df[column][0]


def perform_minimal_harvest_epu(xml_path, output_dir):
    # Before running full eminsight analysis, look for all image files, via xml, mrc or jpg
    # This will declare global searchSupervisorData variables with the file lists, for xml, mrc and jpg
    searchSupervisorAtlas(main.atlas_directory)
    searchSupervisorData(main.epu_directory)

    # Get presets for EPU session xml
    print('')
    print('\033[1m' + 'Finding all presets from EPU session:' + '\033[0m')
    print('')

    xml_presets(xml_path)

    # Get presets specific to acqusition magnifcation which are only contained in an acqusition image xml
    xml_presets_data(searchSupervisorData.xmlData)

    # Get main set up parameters from EPU session xml
    print('')
    print('\033[1m' + 'Finding main EPU session parameters:' + '\033[0m')
    print('')

    # xml_session(xml_path)
    # Call the function with the path to your XML file
    main.masterdf = xml_session(xml_path)

    # Find mics via xml for counting
    searchedFiles = find_mics(main.epu_directory, 'xml')
    if searchedFiles == 'exit':
        print("exiting due to not finding any image xml data")
        exit()
    main.mic_count = len(searchedFiles)

    # Create a deposition file
    deposition_file(xml_path)

def write_mmcif_file(data_list, sessionName):
    """
    pdbx writer is used to write data stored in self.__dataList
    :return written: a boolean; True when pdf writer is finished
    """
    written = False
    depfilepath = main.dep_dir + '/' + sessionName
    if depfilepath:
        mmcif_filename = depfilepath + '_dep.cif'
        with open(mmcif_filename, "w") as cfile:
            pdbx_writer = PdbxWriter(cfile)
            pdbx_writer.write(data_list)
        written = True
    return written


def add_container(data_list, container_id):
    """
    Adds a container with the specified container_id to the data_list.
    """
    container = DataContainer(container_id)
    data_list.append(container)
    return container


def add_category(container, category_id, items):
    """
    Adds a category with the specified category_id and items to the container.
    """
    category = DataCategory(category_id)
    for item in items:
        category.appendAttribute(item)
    container.append(category)


def insert_data(container, category_id, data_list):
    """
    Inserts data_list into the specified category_id within the container.
    """
    cat_obj = container.getObj(category_id)
    if cat_obj is None:
        return

    if all(isinstance(i, list) for i in data_list):
        list_values = [list(t) for t in zip(*data_list)]
        cat_obj.extend(list_values)
    else:
        cat_obj.append(data_list)


def translate_xml_to_cif(input_data, sessionName):
    """
    Translates input XML data into a CIF file.
    """
    if not input_data:
        return False

    input_data_filtered = {key: value for key, value in input_data.items() if key != "?"}
    cif_data_list = []

    # Create a dictionary to store container_id with corresponding categories and values
    container_dict = {}

    # Iterate over the input data and accumulate category names and values
    for key, value in input_data_filtered.items():
        if isinstance(key, str) and key != "?":
            container_id, category = key.split(".")
            if container_id not in container_dict:
                container_dict[container_id] = {}
            if category not in container_dict[container_id]:
                container_dict[container_id][category] = [value]
            else:
                # Handle duplicate categories if needed
                pass

    # Add all accumulated categories and values to the CIF container
    container = add_container(cif_data_list, sessionName)

    for category_name, categories_values in container_dict.items():
        category_list = []
        cif_values_list = []
        for category, cif_values in categories_values.items():
            if category == "date":
                cif_values = [cif_values[0].split(" ")[0]]
            elif category == "accelerating_voltage":
                if cif_values[0] != "?":
                    cif_values = [int(float(cif_values[0]) / 1000)]
            elif category == "nominal_defocus_min":
                if cif_values[0] != "?":
                    cif_values = [int((cif_values[0]) * -1000)]
            elif category == "nominal_defocus_max":
                if cif_values[0] != "?":
                    cif_values = [int((cif_values[0]) * -1000)]
            elif category == "mode":
                if cif_values[0] == "BrightField":
                    cif_values = ["BRIGHT FIELD"]
            elif category == "topology" or category == "material":
                cif_values = [cif_values[0].upper()]
            elif category == "electron_source":
                if cif_values[0] == "FieldEmission":
                    cif_values = ["FIELD EMISSION GUN"]
            elif category == "illumination_mode":
                if cif_values[0] == "PARALLEL":
                    cif_values = ["FLOOD BEAM"]

            category_list.append(category)
            cif_values_list.append(cif_values)
        add_category(container, category_name, category_list)
        insert_data(container, category_name, cif_values_list)

    # Write the modified CIF data to a file
    return write_mmcif_file(cif_data_list, sessionName)


main()