#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
From manual measurements to opensim scale config file

Ideally in a mocap system with 3d coordinates of markers, we would project some of the markers to define the 3 axes of each body (as defined in 'material-scaling-Rosti.pdf', which follows best practices).
    Then, in the scaling tool of opensim we would define which pairs of markers define which axes of each body ('measurement' method as described in https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scale+Setup+File).
    Opensim would calculate the distances, work out the scaling factors and scale the model.

Here, because we only have a limited set of manual measurements (distance between two anatomical landmarks that aren't perfectly aligning with the axes),
    we need to use the 'manualScale' scaling in opensim (where only the "ScaleSet" tag matters) instead of the 'measurements' (which is marker based) one.
    We first need to approximate the sizes by computing the euclidean distances between the landmarks in the generic opensim model,
    and then compare them to the manual measurements in order to find the scaling factors.

Finally, this script will generate the scaling config file and scale the model using opensim.

"""

###############################################################################
# Imports
###############################################################################

import sys, os, copy
CURR_DIR = os.path.dirname(os.path.realpath(__file__))
import pandas, pandas as pd
import numpy as np
import xml.etree.ElementTree as xml

###############################################################################
# Manual lengths from the excel file
###############################################################################

# load the excel file
fileExcel = CURR_DIR+'/INPUTFILE.xlsx'
df = pd.read_excel(fileExcel, engine='openpyxl')

excel_thorax_z = df.iloc[8, 3] # IJ to C7 : 130
excel_thorax_y = df.iloc[7, 3] # IJ to PX : 145
excel_thorax_x = df.iloc[9, 3] # IJ to AC : 165

excel_scapula_x = df.iloc[18, 3] # TS to AA : 125
excel_scapula_y = df.iloc[17, 3] # AI to TS : 130
excel_scapula_z = df.iloc[16, 3] # AA to AC : 150

excel_clavicle_z = df.iloc[19, 3] # SC to AC : 185

excel_humerus_y = df.iloc[13, 3] # GH Centre to mid EpL-Epm (length of humerus) : 255
excel_humerus_z = df.iloc[12, 3] # EpL to EpM (Elbow width) : 85

excel_ulna_y = df.iloc[15, 3] # Mid RS-US to Mid EpL-EpM (ulnar length) : 240

###############################################################################
# Lengths form the generic opensim model
###############################################################################

# load the osim file
#fileDasName = 'model-original-das3_withMarkersScaling.osim'
fileDasName = 'model-das3-MOD2_withMarkersScaling.osim'
fileDas = CURR_DIR+'/'+fileDasName
tree_osim = osim = das = xml.ElementTree(file=fileDas)

# function to read the osim xml file and extract the local coordinates of markers 
def get_marker(targetMarker, verbose=False):

    # Loop through the markers in the opensim file
    for this_marker in tree_osim.findall("./Model/MarkerSet/objects/"):
        name = this_marker.attrib['name']
        
        # If the name matches our target
        if name == targetMarker:

            # Extract the 'location' tag, convert to a numpy array, return it
            locationInBody_str = this_marker.find("location").text
            locationInBody     = np.fromstring(locationInBody_str, dtype=float, sep=' ')
            if verbose: print("Name {} \nLocation {}".format(name,locationInBody))
    
            return locationInBody

# Function to get the length of a vector in 3d space, defined by 2 points
def length(a, b, millimeters=True):
    this_length = np.linalg.norm(a - b)
    if millimeters:
        this_length = this_length * 1000 # from meters to millimeters
    return this_length

# Required markers
ij                = get_marker('ij')
ij_clavicle       = get_marker('ij_clavicle')
c7                = get_marker('c7')
px                = get_marker('px')
ac                = get_marker('ac')
aa                = get_marker('aa')
ac                = get_marker('ac')
ac_thorax         = get_marker('ac_thorax')
ac_scapula        = get_marker('ac_scapula')
ai                = get_marker('ai')
ts                = get_marker('ts')
gu                = get_marker('gu')
EpM               = get_marker('EpM')
EpL               = get_marker('EpL')
centelbow_humerus = get_marker('centelbow_humerus')
mid_RS_US_humerus = get_marker('mid_RS_US_humerus')

osim_thorax_x = length(ij, ac_thorax)   # ideally would be ij_proj_c7 -> ij
osim_thorax_y = length(ij, px)          # ideally would be ij_proj_px -> ij
osim_thorax_z = length(ij, c7)          # ideally would be ij_proj_ac -> ij

osim_scapula_x = length(ts, aa)         # ideally would be ai_proj_aa_z -> ai
osim_scapula_y = length(ai, ts)         # ideally would be ai_proj_ts   -> ai
osim_scapula_z = length(aa, ac_scapula) # ideally would be ai_proj_aa_x -> ai

osim_clavicle_z = length(ij_clavicle, ac)

osim_humerus_y = length(gu,  centelbow_humerus)
osim_humerus_z = length(EpL, EpM)

osim_ulna_y = length(centelbow_humerus, mid_RS_US_humerus)

###############################################################################
# Scaling factors from osim to excel - difference in the two measurements
###############################################################################

# the missing ones are taken from the available factors on other axes

scale_thorax_x = excel_thorax_x / osim_thorax_x
scale_thorax_y = excel_thorax_y / osim_thorax_y
scale_thorax_z = excel_thorax_z / osim_thorax_z

scale_scapula_x = excel_scapula_x / osim_scapula_x
scale_scapula_y = excel_scapula_y / osim_scapula_y
#scale_scapula_z = excel_scapula_z / osim_scapula_z
scale_scapula_z = (scale_scapula_x + scale_scapula_y) / 2 # average of y and z

scale_clavicle_z = excel_clavicle_z / osim_clavicle_z
scale_clavicle_x = scale_clavicle_z # from clavicle z
scale_clavicle_y = scale_clavicle_z # from clavicle z

scale_humerus_y = excel_humerus_y / osim_humerus_y
scale_humerus_z = excel_humerus_z / osim_humerus_z
scale_humerus_x = (scale_humerus_y + scale_humerus_z) / 2 # average of y and z

scale_ulna_y = excel_ulna_y / osim_ulna_y
scale_ulna_x = scale_ulna_y   # from ulna y
scale_ulna_z = scale_ulna_y   # from ulna y

scale_radius_x = scale_ulna_y # from ulna y
scale_radius_y = scale_ulna_y # from ulna y
scale_radius_z = scale_ulna_y # from ulna y

###############################################################################
# Table summary
###############################################################################

print("\n##########################################################")
print("    Excel   |   Osim")
print("##########################################################")
print("Thorax")
print("x = excel {}     |   osim {}   =   scale factor {}".format(excel_thorax_x, osim_thorax_x, scale_thorax_x))
print("y = excel {}     |   osim {}   =   scale factor {}".format(excel_thorax_y, osim_thorax_y, scale_thorax_y))
print("z = excel {}     |   osim {}   =   scale factor {}".format(excel_thorax_z, osim_thorax_z, scale_thorax_z))
print("Scapula")
print("x = excel {}     |   osim {}   =   scale factor {}".format(excel_scapula_x, osim_scapula_x, scale_scapula_x))
print("y = excel {}     |   osim {}   =   scale factor {}".format(excel_scapula_y, osim_scapula_y, scale_scapula_y))
print("z = excel {}     |   osim {}   =   scale factor {}".format(excel_scapula_z, osim_scapula_z, scale_scapula_z))
print("Clavicle")
print("x = excel {}     |   osim {}   =   scale factor {}".format(excel_clavicle_z, osim_clavicle_z, scale_clavicle_z))
print("Humerus")
print("y = excel {}     |   osim {}   =   scale factor {}".format(excel_humerus_y, osim_humerus_y, scale_humerus_y))
print("z = excel {}     |   osim {}   =   scale factor {}".format(excel_humerus_z, osim_humerus_z, scale_humerus_z))
print("Ulna")
print("y = excel {}     |   osim {}   =   scale factor {}".format(excel_ulna_y, osim_ulna_y, scale_ulna_y))
print("##########################################################\n")

###############################################################################
# Generate config file
###############################################################################

# https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scale+Setup+File

# Load generic config file
scale_configFile   = 'config-scale-generic-20211103.xml'
scaleXML = CURR_DIR+"/"+scale_configFile
tree1    = double = xml.ElementTree(file=scaleXML)
root     = tree1.getroot()

### Change the tags

height = 1730 # in mm
weight = 63.5 # in kg

## Model file
scale_modelFile = fileDasName
tree1.find(".ScaleTool/GenericModelMaker/model_file").text = scale_modelFile

## ModelScaler configuration

# apply : true
tree1.find(".ScaleTool/ModelScaler/apply").text = 'true'

# scaling_order = manualScale
tree1.find(".ScaleTool/ModelScaler/scaling_order").text = 'manualScale'

# preserve_mass_distribution = true
tree1.find(".ScaleTool/ModelScaler/preserve_mass_distribution").text = 'true'

# test1 : output_model_file = 1.model_scaled_bones.osim
tree1.find(".ScaleTool/ModelScaler/output_model_file").text = '1.model_scaled_bones.osim'
# test2 : output_scale_file = two.xml
tree1.find(".ScaleTool/ModelScaler/output_scale_file").text = ''

# scaleset : all of our scaling factors

obj = tree1.find(".ScaleTool/ModelScaler/ScaleSet/objects")

for x in obj:

    name = x.find("segment").text 
    
    if name == 'thorax':
        x.find("scales").text = str(scale_thorax_x) + ' ' + str(scale_thorax_y) + ' ' + str(scale_thorax_z)
        #x.find("scales").text = "1 1 4"

    if name == 'scapula_r':
        x.find("scales").text = str(scale_scapula_x) + ' ' + str(scale_scapula_y) + ' ' + str(scale_scapula_z)

    if name == 'clavicle_r':
        x.find("scales").text = str(scale_clavicle_x) + ' ' + str(scale_clavicle_y) + ' ' + str(scale_clavicle_z)

    if name == 'humerus_r':
        x.find("scales").text = str(scale_humerus_x) + ' ' + str(scale_humerus_y) + ' ' + str(scale_humerus_z)

    if name == 'ulna_r':
        x.find("scales").text = str(scale_ulna_x) + ' ' + str(scale_ulna_y) + ' ' + str(scale_ulna_z)

    if name == 'radius_r':
        x.find("scales").text = str(scale_radius_x) + ' ' + str(scale_radius_y) + ' ' + str(scale_radius_z)

## Save the modified config file

# nameNoExt  = configFile.rpartition(".")[0] # https://stackoverflow.com/questions/7351744/split-string-in-to-2-based-on-last-occurrence-of-a-separator#7351789
# temp_config_name = nameNoExt+"_tmp.xml" # ik_outputFolder #tmp1 = file.split("/", 1) ; #tmp2       = tmp1[1][0:-3]
temp_config_name = 'tmp_config_file.xml'
tree1.write(CURR_DIR+"/"+temp_config_name,xml_declaration=True, encoding='utf-8') # save with the header

###############################################################################
# Run opensim's scaling tool and generate a scaled model + visualization API
###############################################################################

import opensim

# Load the model
configFolder = CURR_DIR+"/"
model = opensim.Model(configFolder+scale_modelFile)
print("loaded model at {}".format(configFolder+scale_modelFile))

# Run the scaling tool
print("now will load the config file {}".format(configFolder+temp_config_name))
scale_tool = opensim.ScaleTool(configFolder+temp_config_name)
scale_tool.setDebugLevel(3)
scale_tool.run()

######
# Visualization of scaled model
######

model_scaled = opensim.Model(CURR_DIR+"/1.model_scaled_bones.osim")
for body in model_scaled.getBodyList():
    print(body.getName())
model_scaled.setUseVisualizer(True)
state = model_scaled.initSystem()
viz = model_scaled.updVisualizer().updSimbodyVisualizer()
#viz.setBackgroundColor(osim.Vec3(0)) # white
viz.setGroundHeight(-2)
import time
while 1:
    time.sleep(1)
    model_scaled.getVisualizer().show(state)
























