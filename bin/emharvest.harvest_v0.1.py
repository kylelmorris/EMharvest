#!/usr/bin/env python3

import os
import fnmatch
import math
import argparse
import datetime
import dateutil.parser
import glob
from pathlib import Path
import pandas as pd

import json
from rich.pretty import pprint
from typing import Any, Dict, Tuple
import xmltodict

import hashlib

parser = argparse.ArgumentParser(description="Microscopy Data harvest Script")
parser.add_argument("-e", "--epu", help="EPU session file: Session.dm")
parser.add_argument("-a", "--atlas", help="Atlas session file: ScreeningSession.dm")
parser.add_argument("-o", "--output", help="Output directory for generated reports")
parser.add_argument("-p", "--print", help="Optional: Y = Only print xml and exit")
args = parser.parse_args()

def main():
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

    perform_minimal_harvest(main.epu_xml, output_dir)

def findpattern(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    #result = glob.glob('**/'+str(pattern), recursive=True)
    return result

def roundup(n, decimals=0):
    # https://realpython.com/python-rounding/
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

def checksum(path, out):
    # https://www.quickprogrammingtips.com/python/how-to-calculate-sha256-hash-of-a-file-in-python.html
    filename = path
    sha256_hash = hashlib.sha256()
    with open(filename,"rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
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

    xmlAtlasList = findpattern('Atlas*.xml', path) # You then need to remove any item containing *Data*
    try:
      xmlAtlas = xmlAtlasList[0]
      #print(xml)
    except:
      xmlAtlas = 'None'

    searchSupervisorAtlas.xmlAtlasList = xmlAtlasList
    searchSupervisorAtlas.xmlAtlas = xmlAtlas

    xmlAtlasTileList = findpattern('Tile*.xml', path) # You then need to remove any item containing *Data*
    try:
      xmlAtlasTile = xmlAtlasTileList[0]
      #print(xml)
    except:
      xmlAtlasTile = 'None'

    searchSupervisorAtlas.xmlAtlasTileList = xmlAtlasTileList
    searchSupervisorAtlas.xmlAtlasTile = xmlAtlasTile

    print('Found representative xml file for pulling meta data about Atlas session')
    print('Atlas: '+searchSupervisorAtlas.xmlAtlas)
    print('Atlas tile: '+searchSupervisorAtlas.xmlAtlasTile)
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
      #print(xml)
    except:
      xmlSquare = 'None'
      print('None found')

    print('Finding FoilHole xml')
    xmlHoleList = findpattern('FoilHole*.xml', path)
    xmlHoleList = [ x for x in xmlHoleList if "Data" not in x ] # This will remove items in list containing *Data*, i.e. DataAcquisition xml files
    try:
      xmlHole = xmlHoleList[0]
      print('Done')
      #print(xml)
    except:
      xmlHole = 'None'
      print('None found')

    print('Finding AcquisitionData xml')
    xmlDataList = findpattern('FoilHole*Data*.xml', path)
    try:
      xmlData = xmlDataList[0]
      print('Done')
      #print(xml)
    except:
      xmlData = 'None'
      print('None found')

    print('Finding AquisitionData mrc')
    mrcDataList = findpattern('FoilHole*Data*.mrc', path)
    try:
      mrc = mrcDataList[0]
      print('Done')
      #print(mrc)
    except:
      mrc = 'None'
      print('None found')

    print('Finding AquisitionData jpg')
    jpgDataList = findpattern('FoilHole*Data*.jp*g', path)
    try:
      jpg = jpgDataList[0]
      print('Done')
      #print(jpg)
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
    print('Square: '+searchSupervisorData.xmlSquare)
    print('Hole: '+searchSupervisorData.xmlHole)
    print('Acquisition: '+searchSupervisorData.xmlData)
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
    presets = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"]
    length=len(presets)
    camera = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][0]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"]
    lengthCam=len(camera)

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
        name = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["key"]

        # Get magnifications from image xml, they are not stored in the epu session file
        if  name == 'Atlas':
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

        probeMode = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:ProbeMode"]
        spot = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:SpotIndex"]
        c2 = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:Apertures"]["c:C2Aperture"]["c:Diameter"]
        beamD = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:BeamDiameter"]
        # Deal with two condensor lens systems that don't know beam diameter
        if isinstance(beamD, dict):
            beamD = 0
        else:
            beamD = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:BeamDiameter"]
        beamDmicron = float(beamD)*1e6
        DF = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:Defocus"]
        DFmicron = float(DF)*1e6
        time = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:ExposureTime"]
        epuBin = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:Binning"]["d:x"]

        # Here we face data in lists and not always in the same list position so need to loop to find position
        #for y in range(0, lengthCam):
            #listKeyValue = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"][y]["key"]["#text"]
            #if listKeyValue == 'SuperResolutionFactor':
                #superRes = listKeyValue
                #superResOn = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"][y]["value"]["#text"]

        # Rounding
        spot = round(float(spot))
        c2 = round(float(c2))
        beamDmicron = roundup(float(beamDmicron),1)
        DFmicron = roundup(float(DFmicron),1)
        time = roundup(float(time),2)
        mag = round(float(mag))
        apix = roundup(float(apix),3)

        # Report presets or continue silentily
        print(name)
        print('Nominal magnification: '+str(mag)+' X')
        print('Pixel size: '+str(apix)+' apix')
        print('Probe mode: '+str(probeMode))
        print('Spot: '+str(spot))
        print('C2 apeture: '+str(c2)+' microns')
        print('Beam diameter: '+str(beamDmicron)+' microns')
        print('Defocus: '+str(DFmicron)+' microns')
        print('Exposure time: '+str(time)+' seconds')
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
        if  name == 'Acquisition':
            xml_presets.time = time
            xml_presets.beamD = beamDmicron
            xml_presets.probe = probeMode
            xml_presets.C2 = c2
            xml_presets.spot = spot
            xml_presets.epuBin = epuBin
            xml_presets.mag = mag
            xml_presets.apix = apix
            #xml_presets.superRes = superRes
            #xml_presets.superResOn = superResOn
        if  name == 'AutoFocus':
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
       key = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i]["a:Key"]
       if key == "SuperResolutionFactor":
          j = i
       i = i+1

    try:
       superResBin = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][j]["a:Value"]["#text"]
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
       i = i+1

    # Stage tilt
    stageAlphaRad = getStageTilt(micpath)[0]
    stageBetaRad = getStageTilt(micpath)[1]
    stageAlpha = roundup(math.degrees(float(stageAlphaRad)),2)
    stageBeta = roundup(math.degrees(float(stageBetaRad)),2)

    # Report
    xml_presets_data.superResBin = superResBin
    xml_presets_data.filterSlitWidth = filterSlitWidth
    xml_presets_data.filterSlit = filterSlit
    xml_presets_data.stageAlpha = roundup(float(stageAlpha),1)
    xml_presets_data.stageBeta = roundup(float(stageBeta),1)
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

def xml_session(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    # Location of EPU session directory on which this script was ran
    realPath = os.path.realpath(xml_path)

    # EPU version
    epuId = data["Version"]["@z:Id"]
    epuBuild = data["Version"]["a:_Build"]
    epuMajor = data["Version"]["a:_Major"]
    epuMinor = data["Version"]["a:_Minor"]
    epuRevision = data["Version"]["a:_Revision"]
    epuVersion = str(epuMajor)+'.'+str(epuMinor)+'.'+str(epuRevision)+'-'+str(epuId)+'.'+str(epuBuild)

    # Output format
    doseFractionOutputFormat = data["DoseFractionsOutputFormat"]["#text"]

    # Autoloader slot (starts at 0)
    autoSlot = data["AutoloaderSlot"]
    autoSlotReal = float(autoSlot)+1
    #print('Autoloader position: '+str(autoSlotReal))

    # SampleXml is a list denoed inside [], pull 0 out
    # Then the text of AtlasId is in a key pair in the dictionary denoted by {}, so call that key pair
    # Use the print_epu_xml def to see full [] and {} formatted data structure
    atlasDir = data["Samples"]["_items"]["SampleXml"][0]["AtlasId"]["#text"]
    #print('Atlas location: '+atlasDir)
    #print('')

    # Session name and creation time
    sessionName = xml_sessionName(xml_path)
    #sessionName = data["Name"]["#text"]
    #sessionDate = data["StartDateTime"]
    sessionDate = data["Samples"]["_items"]["SampleXml"][0]["StartDateTime"]
    sessionDateFormat = formatEPUDate(sessionDate)
    #print('EPU session name: '+sessionName)
    #print('EPU session created: '+str(sessionDateFormat))
    #print('EPU session created: '+str(sessionDate))
    #print('')

    # Grid type - lacey or holeycarbon or holeygold
    try:
       gridType = data["Samples"]["_items"]["SampleXml"][0]["GridType"]
    except:
       gridType = 'Unknown'

    # The I0 filter settings may hint at what grid type is being used
    try:
       I0set = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["IsCalibrated"]
    except:
       I0set = 'Unknown'

    try:
       I0MaxInt = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["MaximumIntensity"]
    except:
       I0MaxInt = 'Unknown'

    try:
       I0MinInt = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["MinimumIntensity"]
    except:
       I0MinInt = 'Unknown'

    # Clustering method
    clustering = data["ClusteringMode"]
    #print('Clustering mode: '+clustering)
    clusteringRadius = data["ClusteringRadius"]
    clusteringRadiusMicron = float(clusteringRadius)*1e6
    #print('Clustering radius: '+str(clusteringRadiusMicron)+' microns')
    #print('')
    if clustering == 'ClusteringWithImageBeamShift':
        afisMode = 'AFIS'
        afisRadius = str(clusteringRadiusMicron)
    else:
        afisMode = 'Accrt'
        afisRadius = np.nan

    focusWith = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["FocusWith"]["#text"]
    focusRecurrence = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["Recurrence"]["#text"]
    #print('Focus with: '+focusWith)
    #print('Focus recurrence: '+focusRecurrence)
    #print('')

    delayAfterImageShift = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DelayAfterImageShift"]
    delayAfterStageShift = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DelayAfterStageShift"]
    #print('Delay after Image Shift: '+delayAfterImageShift+' seconds')
    #print('Delay after Stage Shift: '+delayAfterStageShift+' seconds')
    #print('')

    # Send xml dict over to function to get defocus range
    defocusRange=getDefocusRange(data)

    # In some cases the defocus list is a single value, test to deal with
    # Also find max and min defocus values
    if isinstance(defocusRange, str):
        if defocusRange == "xml read error":
           defocusRangeMicron = "xml read error"
           defocusRangeRound = "xml read error"
           defocusMax = "xml read error"
           defocusMin = "xml read error"
        else:
           #defocusRangeMicron = float(defocusRange) * 1e6
           defocusRangeMicron = float(defocusRange) * 1e6
           defocusRangeRound = roundup(defocusRangeMicron, 2)
           defocusMax = min(defocusRangeRound)
           defocusMin = max(defocusRangeRound)
    elif isinstance(defocusRange, list):
        if defocusRange[0] == "xml read error":
           defocusRangeMicron = "xml read error"
           defocusRangeRound = "xml read error"
           defocusMax = "xml read error"
           defocusMin = "xml read error"
        else:
           #defocusRangeMicron = [float(element) * 1e6 for element in defocusRange]
           defocusRangeMicron = [float(element) * 1 for element in defocusRange]
           defocusRangeRound = [roundup(num, 2) for num in defocusRangeMicron]
           defocusMax = min(defocusRangeRound)
           defocusMin = max(defocusRangeRound)
    #print('Defocus range (um): ')
    #print(defocusRangeRound)
    #print('Defocus max: '+str(defocusMax))
    #print('Defocus min: '+str(defocusMin))

    # Report
    xml_session.epuVersion = epuVersion
    xml_session.doseFractionOutputFormat = doseFractionOutputFormat
    xml_session.sessionName = sessionName
    xml_session.sessionDate = sessionDateFormat
    #xml_session.sessionDate = sessionDate
    xml_session.clustering = clustering
    xml_session.clusteringRadius = clusteringRadiusMicron
    xml_session.focusWith = focusWith
    xml_session.focusRecurrence = focusRecurrence
    xml_session.xmlPath = realPath
    # This should get the directory name two levels up from EpuSession.dm and on diamond systems be the DLS Visit
    xml_session.autoSlot = autoSlotReal
    xml_session.atlasDirOrig = atlasDir

    xml_session.defocus = defocusRangeRound
    xml_session.defocusMax = defocusMax
    xml_session.defocusMin = defocusMin

    xml_session.delayImageShift = delayAfterImageShift
    xml_session.delayStageShift = delayAfterStageShift
    xml_session.afisMode = afisMode
    xml_session.afisRadius = afisRadius
    xml_session.gridType = gridType
    xml_session.I0set = I0set
    xml_session.I0MaxInt = round(float(I0MaxInt))
    xml_session.I0MinInt = round(float(I0MinInt))

    print('EPU version: '+ xml_session.epuVersion)
    print('Dose fraction output: '+ xml_session.doseFractionOutputFormat)
    print('EPU Session: '+ xml_session.sessionName)
    print('EPU Session date: '+ str(xml_session.sessionDate))
    print('EPU Session path '+ xml_session.xmlPath)
    print()
    
    print('Clustering mode: '+ xml_session.clustering)
    print('Clustering radius: '+ str(xml_session.clusteringRadius))
    print('Focus mode: '+ xml_session.focusWith)
    print('Focus recurrance: '+ xml_session.focusRecurrence)
    print()
    
    print('Autoloader slot: '+ str(xml_session.autoSlot))
    print('Atlas directory: '+ xml_session.atlasDirOrig)
    print()

    print('Defocus list: '+ str(xml_session.defocus))
    print('Defocus max: '+ str(xml_session.defocusMax))
    print('Defocus min: '+ str(xml_session.defocusMin))
    print()

    print('Image shift delay: '+ xml_session.delayImageShift)
    print('Stage shift delay: '+ xml_session.delayStageShift)
    print('AFIS mode: '+ xml_session.afisMode)
    print('AFIS radius '+ xml_session.afisRadius)
    print('Grid type: '+ xml_session.gridType)
    print('I0 set: '+ xml_session.I0set)
    print('I0 max: '+ str(xml_session.I0MaxInt))
    print('I0 min '+ str(xml_session.I0MinInt))
    print()

    print('\033[1m'+'Finished gathering metadata from main EPU session file'+'\033[0m')
    print('')

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
    #new_date = datetime.datetime.strptime(d,"%Y-%m-%dT%H:%M:%S.%fZ")
    # Reformat date into YY-MM-DD HH-MM-SS
    #return new_date.strftime("%Y-%m-%d %H:%M:%S")
    return datetime.datetime.strptime(epuDate,"%Y-%m-%d %H:%M:%S")

def getDefocusRange(data):
    # Need to deal with cases where it's multi shot or single shot
    templates = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]

    if isinstance(templates, list):
        shotType = 'Multishot'
        # DEV DEV This defocus code might need a revisit if you get lots of errors
        d = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][0]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]
        if d.get("a:double"):
           try:
              df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][0]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:double"]
              #Sometimes the values contain unicode en-dash and not ASCII hyphen
              #df.replace('\U00002013', '-')
           except:
              print('Warning, could not find defocus range in xml file')
              df = ['xml read error']
        else:
            try:
               df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:_items"]["a:double"]
               #Sometimes the values contain unicode en-dash and not ASCII hyphen
            #df.replace('\U00002013', '-')
            except:
               print('Warning, could not find defocus range in xml file')
               df = ['xml read error']
    else:
        shotType = 'Single'
        # There is sometimes a data structure change I think in single shot acqusition, cause currently unknown, check for it
        d = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]
        if d.get("a:double"):
            try:
               df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:double"]
               #Sometimes the values contain unicode en-dash and not ASCII hyphen
               #df.replace('\U00002013', '-')
            except:
               print('Warning, could not find defocus range in xml file')
               df = ['xml read error']
        else:
            try:
               df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:_items"]["a:double"]
               #Sometimes the values contain unicode en-dash and not ASCII hyphen
               #df.replace('\U00002013', '-')
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

def xmlDoseRate(xmlpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(xmlpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # The values are not always in the same list position in a:KeyValueOfstringanyType
    keyValueList = data["CustomData"]["a:KeyValueOfstringanyType"]

    # Loop through the list to find the DoseRate list position
    i = 0
    for value in keyValueList:
       key = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Key"]
       # Krios III has EF selectris falcon, Glacios II has standard Falcon no energy filter, BM is base model?
       if key == "Detectors[BM-Falcon].DoseRate" or key == "Detectors[EF-Falcon].DoseRate":
          keyvalue = i
       i = i+1

    xmlDoseRate = data["CustomData"]["a:KeyValueOfstringanyType"][keyvalue]["a:Value"]["#text"]

    return xmlDoseRate

def find_mics(path, search):
    # Need to have an independent function to find the mics, then move into search_mics to sort them out
    # So find mics can be used independently

    print('Looking for micrograph data in EPU directory using extension: '+search)
    print('')
    # Just get file names
    #files = glob.glob("./Images-Disc1/GridSquare*/Data/*xml")
    #files.sort(key=os.path.getmtime)
    #print("\n".join(files))

    # Old method of finding xml files
    searchedFiles = glob.glob(path+"/**/GridSquare*/Data/*"+search+'*')
    #searchedFiles = glob.iglob(main.epu+"/Images-Disc1/GridSquare*/Data/*"+search)
    if searchedFiles:
       print('Found micrograph data: '+str(len(searchedFiles)))
    else:
       print('No micrographs found with search term: '+search)
       searchedFiles = 'exit'

    # New method of finding xml files
    #searchedFiles = searchSupervisorData.xmlList

    return searchedFiles

def deposition_file(xml):

    # Get EPU session name from main EPU xml file, this is a function
    sessionName = xml_sessionName(xml)

    # This is the data xml metadata file already in a dictionary
    data = searchSupervisorData.xmlDataDict["MicroscopeImage"]

    # Get mag
    xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]
    xmlAPix = float(xmlMetrePix) * 1e10
    xmlAPix = roundup(xmlAPix, 1)

    # Get scope and kV
    model = data["microscopeData"]["instrument"]["InstrumentModel"]
    eV = data["microscopeData"]["gun"]["AccelerationVoltage"]

    # Save doppio deposition csv file
    dictHorizontal1 = {
    'Microscope': model,
    'epuversion': xml_session.epuVersion,
    'date': xml_session.sessionDate,
    'eV': eV,
    'mag': xmlMag,
    'apix': str(xmlAPix),
    'nominal_defocus_min_microns': xml_session.defocusMin,
    'nominal_defocus_max_microns': xml_session.defocusMax,
    'spot_size': xml_presets.spot,
    'C2_micron': xml_presets.C2,
    'Objective_micron': str(xml_presets_data.objective),
    'Beam_diameter_micron': xml_presets.beamD,
    'collection': xml_session.afisMode,
    'number_of_images': main.mic_count
    }
    df1 = pd.DataFrame([dictHorizontal1])

    # Sample data for the second row
    dictHorizontal2 = {
    'Microscope': '_em_imaging.microscope_model',
    'epuversion': '_em_software.name',
    'date': '_em_imaging.date',
    'eV': '_em_imaging.accelerating_voltage',
    'mag': '_em_imaging.nominal_magnification',
    'apix': '?',
    'nominal_defocus_min_microns': '_em_imaging.nominal_defocus_min',
    'nominal_defocus_max_microns': '_em_imaging.nominal_defocus_max',
    'spot_size': '?',
    'C2_micron': '_em_imaging.c2_aperture_diameter',
    'Objective_micron': '?',
    'Beam_diameter_micron': '?',
    'collection': '?',
    'number_of_images': '?'
    }
    df2 = pd.DataFrame([dictHorizontal2])

    # Append the second row to the DataFrame
    df = pd.concat([df1, df2], ignore_index=True, sort=False)
    #df = df1.merge(df2, left_index=0, right_index=0)
    print(df)

    ## Deposition file
    depfilepath = main.dep_dir+'/'+sessionName+'_dep.json'
    checksumpath = main.dep_dir+'/'+sessionName+'_dep.checksum'

    # Human readable deposition file
    #df.to_csv (main.dep_dir+'/'+sessionName+'.dep', index = False, header=True)
    # Manual headings
    headings = ['Value', 'mmCIF']
    df_transpose = df.T
    # Set headings
    df_transpose.columns = headings
    df_transpose.to_csv(main.dep_dir+'/'+sessionName+'_dep.csv', index = True, header=True)

    # Computer readable deposition file
    df_transpose.to_json(depfilepath, index = True)

    # This can be run before doing full analysis of the session directories
    print("Created deposition file")

    # Deposition Checksum
    checksum(depfilepath, checksumpath)

def deposition_file_append_mmcif(df):
   
   return df

def perform_minimal_harvest(xml_path, output_dir):

    # Before running full eminsight analysis, look for all image files, via xml, mrc or jpg
    # This will declare global searchSupervisorData variables with the file lists, for xml, mrc and jpg
    searchSupervisorAtlas(main.atlas_directory)
    searchSupervisorData(main.epu_directory)

    # Get presets for EPU session xml
    print('')
    print('\033[1m'+'Finding all presets from EPU session:'+'\033[0m')
    print('')

    xml_presets(xml_path)

    # Get presets specific to acqusition magnifcation which are only contained in an acqusition image xml
    xml_presets_data(searchSupervisorData.xmlData)

    # Get main set up parameters from EPU session xml
    print('')
    print('\033[1m'+'Finding main EPU session parameters:'+'\033[0m')
    print('')

    xml_session(xml_path)

    # Find mics via xml for counting
    searchedFiles = find_mics(main.epu_directory, 'xml')
    if searchedFiles == 'exit':
        print("exiting due to not finding any image xml data")
        exit()
    main.mic_count = len(searchedFiles)

    # Create a deposition file
    deposition_file(xml_path)

    # Create an mmCIF file


main()