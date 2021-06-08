#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Edits the muscles parameters in the DAS4.4 OSIM model.
Input is a spreadsheet file containing the oxford testing scores.

Dependencies:
pandas, xmltodict, xldr and openpyxl (for excel xls and xlsx files)

"""
import pandas as pd
import os
import xmltodict
import math
CURR_DIR = os.path.dirname(os.path.realpath(__file__))

######################
# Parameters
######################

input_oxford_file = 'oxford_original.xls'
input_osim_file   = 'das4_v4_original.osim'
output_osim_file  = 'generated.osim'

######################
# Read Oxford test into a dict
######################

oxford_data = []
# uses engine specific to spreadsheet file version
filename, file_extension = os.path.splitext(input_oxford_file)
if file_extension.lower() == ".xls":
    oxford_data = pd.read_excel(CURR_DIR+'/'+input_oxford_file) # uses the default xldr engine for excel 97-2003 files
elif file_extension.lower() == ".xlsx":
    oxford_data = pd.read_excel(CURR_DIR+'/'+input_oxford_file, engine='openpyxl') # uses the openpyxl engine for excel 2007-365 files
elif file_extension.lower() == ".ods":
    oxford_data = pd.read_excel(CURR_DIR+'/'+input_oxford_file, engine='odf') # uses the pdf engine for libreoffice/openoffice files

# create a dict from https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html
oxford_values_tmp = oxford_data.to_dict('index')
oxford_dict = {}
for x in oxford_values_tmp:
    name  = oxford_values_tmp[x]['Testing']
    score = oxford_values_tmp[x]['Oxford score (/5)']
    if math.isnan(score) :
        pass #print('nan detected, skipping')
    else: 
        oxford_dict[name] = score # adding score to dict

# Correspondence table (movement->muscles)
table = {
    'Shoulder Vertical Flexion'     : ['pect_maj_c','pect_maj_t','delt_clav','bic_l','coracobr'],
    'Shoulder Vertical Extension'   : ['lat_dorsi','ter_maj','ter_min','pect_maj_t','delt_scap','tric_long','infra'],
    'Shoulder ABduction'            : ['supra','delt_clav','delt_scap','trap_clav','trap_scap','lev_scap'],
    'Shoulder Horizontal ABduction' : ['delt_scap'],
    'Shoulder ADduction'            : ['coracobr','pect_maj_t','pect_maj_c','lat_dorsi','ter_maj','rhomboid','serr_ant'],
    'Shoulder Horizontal ADduction' : ['pect_maj_c','pect_min','delt_clav'],
    'Shoulder Medial Rotation'      : ['subscap','ter_maj','lat_dorsi','pect_maj_t','pect_maj_c','delt_clav'],
    'Shoulder Lateral Rotation'     : ['ter_min','infra','delt_scap'],
    'Elbow Flexion'                 : ['bic_l','bic_b','brachialis','brachiorad','pron_teres_1','pron_teres_2'],
    'Elbow Extension'               : ['tric_long','tric_med','tric_lat','anconeus'],
    'Wrist Flexion'                 : ['FCU','FCR','PL'],
    'Wrist Extension'               : ['ECU','ECRL','ECRB'],
    'Wrist Radial deviation'        : ['FCR'],
    'Wrist Ulnar deviation'         : ['FCU','ECU'],
    'Wrist Pronation'               : ['pron_teres_1','pron_teres_2','pron_quad','brachiorad'],
    'Wrist Supination'              : ['supinator','bic_l','bic_b','brachiorad']
    }

######################
# Load the original osim file
######################

myxml = open(CURR_DIR+'/'+input_osim_file).read()
osim_file_loaded = xmltodict.parse(myxml)
osim_muscle_data = osim_file_loaded['OpenSimDocument']['Model']['ForceSet']['objects']['Schutte1993Muscle_Deprecated'] # List of dicts (1dict / muscle)

######################
# Read through table, get the score of each muscle from oxford_data, then modify in osim file
######################

# Read through groups of muscles (from table)
for muscle_group in table:
    print("\n1) Doing group of muscles: "+muscle_group)

    # What is the oxford score for this muscle group?
    oxford_score = int(oxford_dict[muscle_group])
    oxford_ratio = 0
    if oxford_score == 0:
        oxford_ratio = 0.00 # 0/5 No contraction
    if oxford_score == 1:
        oxford_ratio = 0.01 # 1/5 Visible/palpable muscle contraction but no movement
    if oxford_score == 2:
        oxford_ratio = 0.20 # 2/5 Movement with gravity eliminated
    if oxford_score == 3:
        oxford_ratio = 0.50 # 3/5 Movement against gravity only
    if oxford_score == 4:
        oxford_ratio  = 0.75 # 4/5 Movement against gravity with some resistance
    if oxford_score == 5:
        oxford_ratio = 1.00 # 5/5 Movement against gravity with full resistance

    # Read through muscles of this group of muscles (from table)
    for muscle in table[muscle_group]:
        print("\n2) Processing muscle {}".format(muscle))

        # Find this muscle in the osim file
        print("looking for it in osim file...")
        for x in osim_muscle_data:

            name = x['@name']

            # When found, edit the osim file wth the score (% of current value)
            # There can be several muscle fibres to represent the same muscle
            if muscle in name:

                print('Muscle ({}) found in osim file (under name of: {})'.format(muscle,name))              

                # Get the original values
                fmax_original = float(x['max_isometric_force'])

                # Now, decide the new values based on the oxford score
                fmax_new = fmax_original * oxford_ratio
                # We also add "min_control" and "max_control" tags (from 0 to 1)
                # min_control is always 0, max control is proportional to the oxford score
                min_control_new = 0
                max_control_new = 1 * oxford_ratio

                # Save it to the list of dicts
                x['max_isometric_force'] = fmax_new
                x['min_control'] = min_control_new
                x['max_control'] = max_control_new

                print("Max iso force goes from {} to {} (ratio {} from oxford score {}). Min_control is {} and max_control is {}".format(fmax_original,fmax_new,oxford_ratio,oxford_score, min_control_new, max_control_new))

######################
# Saving as new file
######################
#print(xmltodict.unparse(test, pretty = True))

try:
    f = open(CURR_DIR+'/'+output_osim_file, "w")
    print(xmltodict.unparse(osim_file_loaded, pretty = True), file=f)
    f.close()
    print("File saved as {}".format(CURR_DIR+'/'+output_osim_file))
except:
    print("Error saving the file ({})".format(CURR_DIR+'/'+output_osim_file))