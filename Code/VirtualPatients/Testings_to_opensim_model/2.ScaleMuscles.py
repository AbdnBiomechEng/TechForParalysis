#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Info:
Edits the muscles parameters in the DAS4.4 OSIM model.
Input is a spreadsheet file containing the oxford testing scores.
Output is a DAS4.4 osim file with scaled muscle parameters of activation (max_activation) and maximum force (fmax ou max_isometric_force).

Dependencies:
pandas, xmltodict, sympy, xldr and openpyxl (for excel xls and xlsx files respectively)

#############
Walkthrough:
#############

The input spreadsheet contains two tables.
The first table has the names of group of muscles performing motions around the shoulder, elbow and wrist (e.g. "shoulder abduction"), as well as and Oxford score performed without stimulation (0/5).
The second table has the names of individual muscles that are tested with electrical stimulation (0/5).

This code associates [groups of muscles] to [lists of individual muscles],
and works out each muscle's maximum force and maximum activation based on the oxford score(s) available.

We have the following relations:
muscle_force = theoretical_max_force_for_this_muscle * muscle_weakness_ratio
with muscle_weakness = oxford_ratio_stimul / fes_activation

max_control  = innervation factor (if it's a stimulated muscle)
with innervation_factor = (oxford_ratio_vol * fes_activation) / oxford_ratio_stimul

###########################
Example 1: with FES testing
For example, when processing the muscle delt_scap_1,
if the participant has a weak shoulder Abduction Oxford score without stimulation  (oxford_score_volunt = 1, therefore oxford_ratio_volunt = 0.1)
and it has a stronger shoulder Abduction Oxford score when stimulating this muscle (oxford_score_stimul = 3, therefore oxford_ratio_stimul = 0.5)

Then the maximum muscle force of this muscle is:
muscle_force = theoretical_max_force_for_this_muscle * muscle_weakness_ratio
with muscle_weakness = oxford_ratio_stimul / fes_activation
with fes_activation being a constant of how much excitation the FES can provide (simulate) to the muscle (around 0.7 or 0.8)
// muscle_force = theoretical_max_force_for_this_muscle * (oxford_ratio_stimul / fes_activation)
// muscle_force = 289                                   * (0.5                 / 0.8           )
// muscle_force = 180.625

And the maximum ability of the nerve of this muscle to send excitation is:
max_control  = innervation factor (if it's a stimulated muscle)
with innervation_factor = (oxford_ratio_vol * fes_activation) / oxford_ratio_stimul
// max_control = (oxford_ratio_vol * fes_activation) / oxford_ratio_stimul)
// max_control = (0.1              * 0.8           ) / 0.5                )
// max_control = 0.16

###########################
Example 2: without FES testing
For example. when processing the muscle supra_1
the participant has a weak shoulder Abduction Oxford score without stimulation (oxford_score_volunt = 1, therefore oxford_ratio_volunt = 0.1)

We have to either assume that the muscle weakness is equals to the nerve limitation,
!2DO! or have some sort of reference based on the level of injury in order to distribute unevenly the muscle_weakness_factor and the innervation_factor

The maximum muscle force of this muscle is:
muscle_force = theoretical_max_force_for_this_muscle * muscle_weakness_ratio
with muscle_weakness = oxford_ratio_volunt
// muscle_force = theoretical_max_force_for_this_muscle * oxford_ratio_volunt
// muscle_force = 120                                   * 0.1
// muscle_force = 12

And the maximum ability of the nerve of this muscle to send excitation is:
max_control  = innervation factor
with innervation_factor = oxford_ratio_volunt
// max_control = 0.1

#######
2do
#######

[o] Choose how to handle muscles that have no FES testing : either a straight scaling factor based on the Oxford score, or some sort of decision-tree based on the innervation root

"""
import pandas as pd
import os
import xmltodict
import math
import sympy
CURR_DIR = os.path.dirname(os.path.realpath(__file__))

######################
# Parameters
######################

input_oxford_file = 'oxford_original.xls'
#input_osim_file   = 'das4_v4_original.osim'
input_osim_file = '1.model_scaled_bones.osim'
output_osim_file  = '2.model_scaled_muscles.osim'
fes_activation    = 0.8 # constant defining the activation level simulated under FES stimulation

######################
# Read Oxford test into dicts
######################

# Read spreadsheet with engine specific to file version
engine = ''
filename, file_extension = os.path.splitext(input_oxford_file)
if file_extension.lower() == ".xls":
    engine = None # uses the default xldr engine for excel 97-2003 files
elif file_extension.lower() == ".xlsx":
    engine='openpyxl' # uses the openpyxl engine for excel 2007-365 files
elif file_extension.lower() == ".ods":
    engine='odf' # uses the pdf engine for libreoffice/openoffice files
oxford_data = pd.read_excel(CURR_DIR+'/'+input_oxford_file, engine=engine) # uses the odf engine for libreoffice/openoffice files

# create a dict of the oxford score performed without FES stimulation (from https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_dict.html)
oxford_manual_values_tmp = oxford_data.to_dict('index')
oxford_manual_dict = {}
for x in oxford_manual_values_tmp:
    name  = oxford_manual_values_tmp[x]['Testing without stimulation']
    score = oxford_manual_values_tmp[x]['Oxford score (/5)']
    if math.isnan(score) :
        pass #print('nan detected, skipping')
    else: 
        oxford_manual_dict[name] = score # adding score to dict

# create a dict for the oxford score performed with FES stimulation
oxford_stimul_values_tmp = oxford_data.to_dict('index')
oxford_stimul_dict = {}
for x in oxford_stimul_values_tmp:
    muscle_name  = oxford_stimul_values_tmp[x]['Testing with stimulation']
    score = oxford_stimul_values_tmp[x]['Stim score (oxford) (/5)']
    if math.isnan(score) :
        pass #print('nan detected, skipping')
    else: 
        oxford_stimul_dict[muscle_name] = score # adding score to dict

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
osim_muscle_data = osim_file_loaded['OpenSimDocument']['Model']['ForceSet']['objects']['Schutte1993Muscle_Deprecated'] # List of dicts (one dict per muscle)

######################
# Read through table (correspondence movement/muscles), get the score of each muscle from oxford_data (the loaded spreadsheet), then modify muscle parameters in osim file
######################

# Define an Oxford Ratio based on the Oxford Score
def get_oxford_ratio(oxford_score):
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
    return oxford_ratio

# Read through groups of muscles (from table)
for muscle_group in table:
    print("\n1) Doing group of muscles: "+muscle_group)

    # What are the manual oxford score and ratio for this muscle group?
    oxford_score_manual = int(oxford_manual_dict[muscle_group])
    oxford_ratio_manual = get_oxford_ratio(oxford_score_manual)

    # Read through muscles of this group of muscles (from table)

    for muscle in table[muscle_group]:
        if muscle == "delt_scap": print("\n\n\n\n\n##############################################\n\n\n\n\n\n\n")

        print("\n2) Processing muscle {}".format(muscle))

        # Find this muscle in the osim file
        print("looking for it in osim file...")
        for x in osim_muscle_data:

            name = x['@name']

            # When found, edit the muscle parameters of the osim file
            # There can be several muscle fibres to represent the same muscle so we do a partial match on the name
            if muscle in name:

                if muscle == "delt_scap":
                    print("\n")
                    print(x)

                print('Muscle ({}) found in osim file (under name of: {})'.format(muscle,name))    

                # Get the original max force value (fmax)
                print("max iso force is: {}".format(x['max_isometric_force']))
                fmax_original = float(x['max_isometric_force'])

                # Initialise the new fmax and min/max_control
                fmax_new = min_control_new = max_control_new = ''
                muscle_weakness_ratio = ''
                innervation_factor = ''
                has_fes_testing = False

                # We need to know if this specific muscle has been tested with FES.
                muscle_has_fes_testing = False                
                for name_tmp in oxford_stimul_dict:
                    if muscle in name_tmp:
                        muscle_has_fes_testing = True


                # If it's the case, we use the Oxford Testing with FES to determine the new fmax and max_control
                if muscle_has_fes_testing:
                    print("This muscle has a testing with FES. Using it to scale muscle parameters.")
                    oxford_score_stimul = oxford_stimul_dict[name_tmp]
                    print("Its oxford_score_stimul is {}".format(oxford_score_stimul))

                    #the maximum muscle force of this muscle is:
                    #muscle_force = theoretical_max_force_for_this_muscle * muscle_weakness_ratio
                    #with muscle_weakness = oxford_ratio_stimul / fes_activation
                    #with fes_activation being a constant of how much excitation the FES can provide (simulate) to the muscle (around 0.7 or 0.8)

                    # The estimated force ability fmax_new is based on the oxford ratio under stimulation
                    oxford_ratio_stimul = get_oxford_ratio(oxford_score_stimul)
                    muscle_weakness_ratio = oxford_ratio_stimul / fes_activation
                    fmax_new = fmax_original * muscle_weakness_ratio

                    #And the maximum ability of the nerve of this muscle to send excitation is:
                    #max_control  = innervation factor (if it's a stimulated muscle)
                    #with innervation_factor = (oxford_ratio_vol * fes_activation) / oxford_ratio_stimul

                    # The estimated innervation ability is:

                    innervation_factor = (oxford_ratio_manual * fes_activation) / oxford_ratio_stimul
                    max_control_new = innervation_factor
                    min_control_new = 0

                    print("\nMax iso force goes from {} to {} (ratio {} from oxford score {}). Min_control is {} and max_control is {}\n".format(fmax_original,fmax_new,oxford_ratio_stimul,oxford_score_stimul, min_control_new, max_control_new))

                    #fmax_new = (fmax_original * oxford_ratio_stimul) / fes_activation
                    # The estimated innervation ability is:
                    # We know oxford_score_manual and oxford_score_stimul
                    # oxford_score_manual = fmax_new * max_control_new
                    # oxford_score_stimul = fmax_new * fes_activation
                    # Solve for this equation:

                    # From https://stackoverflow.com/questions/6789927/is-there-a-python-module-to-solve-linear-equations
                    """
                    print("We have two unknown : fmax_new and max_control_new")
                    
                    print("oxford_score_manual = fmax_new * max_control_new")
                    print("{} = fmax_new * max_control_new".format(oxford_score_manual))
                    print("oxford_score_stimul = fmax_new * fes_activation")
                    print("{} = fmax_new * {}".format(oxford_score_stimul,fes_activation))

                    print("Rearanged as :")
                    print("0 = fmax_new * max_control_new - oxford_score_manual")
                    print("0 = fmax_new * max_control_new - {}".format(oxford_score_manual))
                    print("and")
                    print("0 = fmax_new * fes_activation - oxford_score_stimul")
                    print("0 = fmax_new * {} - {}".format(fes_activation,oxford_score_stimul))
                    """
                    """                    
                    # To solve manually we would do, for fes_activation = 0.8, oxford_score_stimul = 4, oxford_score_manual = 2:
                    # oxford_score_manual = fmax_new * max_control_new
                    # oxford_score_stimul = fmax_new * fes_activation
                    # soit
                    # equation A : 2 = fmax_new * max_control_new
                    # equation B : 4 = fmax_new * 0.8              --> fmax_new = 4 / 0.8 = 5
                    # solving equation B allows us to can complete equation A:
                    # equation A : 2 = 5 * max_control_new         --> max_control_new = 2 / 5 = 0.4
                    """

                    # Solve for this equation:
                    """
                    x0, x1 = sympy.symbols(['x0', 'x1']) # x0 is fmax_new , x1 is max_control_new and the only one we consider here
                    result = sympy.solve([x0 * x1 - oxford_score_manual, x0 * fes_activation - oxford_score_stimul], [x0, x1])
                    max_control_new = round(float(result[0][1]),2) # result rounded to two decimals
                    min_control_new = 0
                    
                    print("\nMax iso force goes from {} to {} (ratio {} from oxford score {}). Min_control is {} and max_control is {}\n".format(fmax_original,fmax_new,oxford_ratio_stimul,oxford_score_stimul, min_control_new, max_control_new))
                    """
                # If not, we fall back to using the Oxford Testing without FES
                else:
                    print("This muscle does not have a testing with FES, using non-FES manual testing to scale muscle parameters.")
                    
                    #We have to either assume that the muscle weakness is equals to the nerve limitation,
                    #!2DO! or have some sort of reference based on the level of injury  to distribute unevenly 

                    #The maximum muscle force of this muscle is:
                    #muscle_force = theoretical_max_force_for_this_muscle * muscle_weakness_ratio
                    #with muscle_weakness = oxford_ratio_volunt

                    muscle_weakness_ratio = oxford_ratio_manual
                    fmax_new = fmax_original * muscle_weakness_ratio

                    #And the maximum ability of the nerve of this muscle to send excitation is:
                    #max_control  = innervation factor
                    #with innervation_factor = oxford_ratio_volunt

                    innervation_factor = oxford_ratio_manual
                    max_control_new = innervation_factor
                    min_control_new = 0

                    print("\nMax iso force goes from {} to {} (ratio {} from oxford score {}). Min_control is {} and max_control is {}\n".format(fmax_original,fmax_new,oxford_ratio_manual,oxford_score_manual, min_control_new, max_control_new))

                # Save the calculated values to the list of dicts
                x['max_isometric_force_original'] = fmax_original
                x['max_isometric_force'] = fmax_new
                x['min_control'] = min_control_new
                x['max_control'] = max_control_new
                x['muscle_weakness_ratio'] = muscle_weakness_ratio
                x['innervation_factor'] = innervation_factor
                x['muscle_has_fes_testing'] = muscle_has_fes_testing

######################
# Saving as new file
######################

try:
    f = open(CURR_DIR+'/'+output_osim_file, "w")
    print(xmltodict.unparse(osim_file_loaded, pretty = True), file=f)
    f.close()
    print("File saved as {}".format(CURR_DIR+'/'+output_osim_file))
except:
    print("Error saving the file ({})".format(CURR_DIR+'/'+output_osim_file))