#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Optimisation to compute muscle weakness (atrophy and reduction in nerve input) based on volutional (and stimulated) testings

# 1 : load das3 and target values (from experiment)

# 2 : perform optimisation for selected muscles, with flag for STIM or not (2 different cost functions)

# 3 : save as new osim file with muscles scaled

2do:
    [o] Lock ROM : fix the xml parser

"""

# Spyder only: Rest the vars at run script
from IPython import get_ipython
try:
    __IPYTHON__
except NameError:
    #print("Not in IPython")
    pass
else:
    print("ipython console: resetting")
    def __reset__():
        get_ipython().magic("reset -sf")
    __reset__()

import os
from time import sleep
import numpy as np
np.seterr(divide='ignore', invalid='ignore') # do not display error messages likely to happen during optimisation
import pandas as pd
from collections import OrderedDict
from scipy.optimize import minimize, Bounds, LinearConstraint
import re # regex

import opensim
import opensim as osim
import xmltodict

import plotly.graph_objects as go
import plotly.express as px
import plotly
from termcolor import colored

CURR_DIR = os.path.dirname(os.path.realpath(__file__))
os.chdir(os.path.abspath(os.path.dirname(__file__)))   # change current working directory

def show(data):
    print(np.around(data, 4))

def title(text, color="blue"):
    n = len(text) + 4
    horizontal = '#' * n
    print("\n##############################\n# {}\n##############################\n".format(colored(text, color, attrs=['bold']),))
    #print("\n{}\n# {} #\n{}\n".format(horizontal,text,horizontal))

#####################################
#
# Step 1 : compute the moment arm and max force of the target muscles from the opensim file, with the elbow in 90° flexion and in neutral prono-supination
# Each muscle is defined as several muscles in the opensim model
#
#####################################

title("Step 1 - load model, parameters, and compute muscles moment arms")

#####################################
# ------ Model parameters
#####################################

# Scale factor for Fmax : an manual average computed from bone scaling
# There are several options to do this, we set it at 1.4 for now.
scaling_fmax_main = 1.4

# Model
#input_osim_file  = '../../Data/Models/das3.osim'
input_osim_file  = '1.model_scaled_bones.osim'   # the initial file to load : scaled skeleton
input_osim_file_scaled = 'das3_scaled_temp.osim' # the temporary file that scaled all muscle forces by same factor
output_osim_file = '2.model_scaled_bonesAndMuscles.osim'             # the optimised file, that has modified each muscle's max force and activation

#####################################
# ------ Scale the muscle forces (based on geometry and one scaling factor)
#####################################

scaleInitialMuscleForce = True

if scaleInitialMuscleForce:

    # Load the original osim file in external xml parser
    
    myxml            = open(CURR_DIR + "/" + input_osim_file).read()
    osim_file_loaded = xmltodict.parse(myxml)
    osim_muscle_data = osim_file_loaded["OpenSimDocument"]["Model"]["ForceSet"]["objects"]["Schutte1993Muscle_Deprecated"] # List of dicts (one dict per muscle)
    
    for x in osim_muscle_data:
    
        # get the original fmax, and set the new scaled force
        name                     = x["@name"]
        fmax_original            = float(x["max_isometric_force"])
        x["max_isometric_force"] = fmax_original * scaling_fmax_main
        x["max_isometric_force_original"] = fmax_original
    
    # Save as a new file

    try:
        f = open(CURR_DIR+'/'+input_osim_file_scaled, "w")
        print(xmltodict.unparse(osim_file_loaded, pretty = True), file=f)
        f.close()
        print("File saved as {}".format(CURR_DIR+'/'+input_osim_file_scaled))
    except:
        print("Error saving the file ({})".format(CURR_DIR+'/'+input_osim_file_scaled))

#####################################
# ------ Load the (scaled) model, and initialise states
#####################################

# Load in external parser
myxml            = open(CURR_DIR + "/" + input_osim_file_scaled).read()
osim_file_loaded = xmltodict.parse(myxml)
osim_muscle_data = osim_file_loaded["OpenSimDocument"]["Model"]["ForceSet"]["objects"]["Schutte1993Muscle_Deprecated"] # List of dicts (one dict per muscle)

# Load in osim
model            = opensim.Model(input_osim_file_scaled)
visualizeFlag    = False
model.setUseVisualizer(visualizeFlag)

# Initialize the model and create the states
state = model.initSystem()

# Set the elbow flexion to 90° and supination to maximum

elbow = model.getJointSet().get("hu") # Get the elbow joint
angle = 90
elbow.getCoordinate().setValue(state, np.deg2rad(angle)) # 0.5 * osim.SimTK_PI if using radians and pi

ulna_radius = model.getJointSet().get("ur")    # Get the ulna/radius joint # coordinates : PS_y
ulna_radius.getCoordinate().setValue(state, 0) # Set the coordinates to max supination (values go from 0 to 3 so set to 0)

# Configure the visualizer for graphical confirmation of position

if visualizeFlag: 
    viz = model.updVisualizer().updSimbodyVisualizer()
    viz.setBackgroundColor(osim.Vec3(0)) # white
    viz.setGroundHeight(-2)
    model.getVisualizer().show(state)

# Main function calculating the moment arm and fmax at this angle, then performing the optimisation

def performOptimisation(jointCoordinatesName, muscles_all, muscles_vol, muscles_stim, jointMoment_vol, jointMoment_stim, muscles_vol_osimNames, muscles_stim_osimNames, optimise=True, withFES=True):
    #####################################
    # ------ Get moment arm of all the muscles involved in this motion (jointCoordinates)
    #####################################
    
    jointCoordinates = model.updCoordinateSet().get(jointCoordinatesName) # define coordinates (ex: 'EL_x' for elbow flexion)
    
    model.equilibrateMuscles(state) # make sure states are in equilibrium
    
    def get_momentArm_forMuscles(muscles_list, data, verbose=False):
        # locate the muscles in the xml document
        for muscle in muscles_list:
        
            if verbose: print("\nLooking for muscle {}".format(muscle))
        
            for x in osim_muscle_data:
        
                name = x["@name"]
    
                # There can be several muscle 'fibres' to represent the same muscle so we do a partial match on the name
                string_match = re.search("^"+muscle, name) # regex to check if our target muscle is at beginning of muscle name
    
                if string_match:
        
                    # When found a matching muscle, compute moment arm and save in the momentarm_dict
                    if verbose: print("found matching name: {}".format(name))
        
                    # get the muscle
                    this_muscle = model.getMuscles().get(name)
                    
                    # compute moment arm
                    this_moment_arm = abs(this_muscle.computeMomentArm(state, jointCoordinates))
                    #print('moment arm is {}'.format(this_moment_arm))
                    # keep it positive to have consistent graphs (otherwise elbow extensors' new Fmax will appear negative)
                    # this_moment_arm = abs(this_moment_arm)
                    if verbose: print('moment arm is {}'.format(this_moment_arm))
        
                    # get the original Fmax
                    fmax_original = float(x["max_isometric_force"])
        
                    # compute fmax at this angle            
                    this_muscle.setActivation(state, 1) # Set the activation to 1 (madesx)
                    model.equilibrateMuscles(state)     # make sure states are in equilibrium
                    active  = this_muscle.getActiveForceLengthMultiplier(state)
                    passive = this_muscle.getPassiveForceMultiplier(state)
                    this_maxisoforce = this_muscle.getMaxIsometricForce() # get max isometric force
                    this_pennationAngle = this_muscle.getCosPennationAngle(state)
                    fmax_thisAngle = this_maxisoforce * (active + passive) * this_pennationAngle
                    if verbose: print("Fmax {} at angle {} is {}".format(this_maxisoforce, angle, fmax_thisAngle))
            
                    # save to dict
                    data[name] = {'moment_arm' : this_moment_arm, 'fmax_thisAngle' : fmax_thisAngle, 'fmax_original' : fmax_original}
        
        '''
        # Result example: data =
        
        {'bic_l'       : {'moment_arm': 0.05143613931184783,  'fmax_thisAngle': 332.4, 'fmax_original' : 347.0},
         'bic_b_1'     : {...},
         'bic_b_2'     : {...},
         'brachialis_1': {...},
        '''
    
    # Get moment of the muscles, and add to data_all
    get_momentArm_forMuscles(muscles_all, data_all)
    
    #####################################
    #
    # Step 2 : Optimise
    # from the data (joint moment under volutional and stimulated control, and fmax of each muscle at this angle),
    # run optimisation to compute the weakness factor (reduction in nerve input as well as muscle atrophy)
    #
    #####################################
    
    title("Step 2 - Optimise")
    
    #####################################
    # ------ First step is to compute the weakness factor for the muscle(s) that have been tested under stimulation
    #####################################
    
    # We need two lists representing parameters of the stimulated muscles:
    Fmax_stim = [] # the max force at this angle
    r_stim = []    # the moment arm at this angle
    
    for muscle in muscles_stim: # for each muscle that is stimulated
        
        for key, value in data_all.items(): # look into our data
            name = key
            string_match = re.search("^"+muscle, name) # regex to check if our target muscle is at beginning of muscle name
            
            if string_match: # if the name matches:
                
                # add the parameters to our lists
                Fmax_stim.append(value['fmax_thisAngle'])
                r_stim.append(value['moment_arm'])
                
                # also save the name of each muscle as present in the osim file
                muscles_stim_osimNames.append(name)
    
    # Given that:
    # jointMoment_stim = sum( activation * (scale * Fmax_stim) * r_stim )
    # with activation = 0.5 under FES stimulation
    # we can work out the variable "scale" =  jointMoment_stim / (0.5 * Fmax_stim * r_stim)
    
    Sstim = jointMoment_stim / ( sum( 0.5 * np.array(Fmax_stim) * np.array(r_stim) ) )
    #print("Sstim = {}".format(Sstim))
    
    # Considering this scaling factor, the maximum moment these muscles can generate in a voluntary situation is:
    
    max_force_patho_stim = moment_muscles_stim = sum( Sstim * np.array(Fmax_stim) * Sstim * np.array(r_stim) ) # = 6.659993919241959 Nm. So there's only ~2.1 Nm to be shared with the other muscles!
    #print(max_force_patho_stim)
    # This is the contribution of the stimulated muscles to the jointMoment during a VOLUNTARY maximum contraction (without stimulation).
    # The remaining of the jointMoment is to be shared between the other muscles.
    # The next step is therefore to optimise the values that scale both the activation and the fmax of the remainng muscles, in order to match the total voluntary jointmoment that has been recorded.
    
    #####################################
    # ------ We need two lists representing parameters of the voluntary muscles:
    #####################################
    
    Fmax_vol = [] # the max force at this angle
    r_vol = []    # the moment arm at this angle
    
    for muscle in muscles_vol: # for each muscle that is stimulated
        
        for key, value in data_all.items(): # look into our data
            name = key
            string_match = re.search("^"+muscle, name) # regex to check if our target muscle is at beginning of muscle name
            
            if string_match: # if the name matches:
                
                # add the parameters to our lists
                Fmax_vol.append(value['fmax_thisAngle'])
                r_vol.append(value['moment_arm'])
    
                # also save the name of each muscle as present in the osim file
                muscles_vol_osimNames.append(name)
    
    number_muscles_vol = len(r_vol)
    
    #####################################
    # ------ Run the optimisation
    #####################################
    
    # in some cases (e.g. participant 1 elbow extension), we do not want 

    if optimise:
        # Objective function to minimize / cost function
        def fun(x):
            result = sum(1/x**2) # maximise the square of the weakness factor (squared as appears twice : once for the activation, once for the fmax)
            return result
        
        # Equality constraints
        def equality_constraint(x):
        
            # The equality constraint is that the sum of all the muscles moments has to match the joint moment.
            # We already know the contribution of some muscles (the ones that have been tested under stimulation) : max_force_patho_stim   
            # Therefore, the equality constraint is : max_force_patho_stim + moment_all_other_muscles = jointMoment
            # or max_force_patho_stim + moment_all_other_muscles - jointMoment = 0
        
            result = moment_muscles_stim + sum( x**2 * Fmax_vol * r_vol ) - jointMoment_vol
        
            return result
        
        constraint1 = ({'type':'eq','fun': equality_constraint})
        
        def equality_constraint_noFES(x):
        
            # this one has no FES included : all the muscles that have to match the jointMoment are unknown

            result = sum( x**2 * Fmax_vol * r_vol ) - jointMoment_vol
        
            return result
        
        constraint2 = ({'type':'eq','fun': equality_constraint_noFES})

        # now depending on the type of analysis, we use different constraints
        constraints = []
        if withFES:
            constraints = [constraint1]
        else:
            constraints = [constraint2]
        
        # Bounds for the output values
        bounds_lower   = np.zeros([number_muscles_vol,]) # zeros
        bounds_upper   = np.ones([number_muscles_vol,])  # ones
        bounds         = Bounds(bounds_lower, bounds_upper)
                
        def run():
            
            # Initial guess : as many values as we have muscles with unknown scaling factors
            # they are in the function run() so that we can generate new starting values in case it fails
            x0 = np.random.uniform(low=0, high=1, size=number_muscles_vol) # between 0 and 1, as many values as muscles : so we have voluntary and known_stimulated_muscle
            #x0 = np.ones([number_muscles_vol,])

            global res
            res = minimize(
                fun,
                x0=x0,
                constraints=constraints,
                bounds=bounds,
                options={'disp':True}
            )

        # make a loop to re-run the optimisation if it has failed :
        optSuccess = False
        nb_tries=0
        while not optSuccess:
            nb_tries+=1
            run()
            optSuccess = res.success
            if not optSuccess:
                title("OPT FAILED", 'red')
                print("try", nb_tries)
        
        title("OPT SUCCESS", 'green')
        print("try", nb_tries)

        #####################################
        # ------ Results
        #####################################
        
        title("Results")
        
        ######
        # Sanity check 
        ######
        
        # Using these scaling factors on the unknown muscles, and adding the known muscles contribution, should match the voluntary joint moment (i.e. follow the equality constraint!)
        result_cumulated_unknownMuscleMoments = + sum( res.x**2 * Fmax_vol * r_vol )
        isClose = np.isclose(moment_muscles_stim + result_cumulated_unknownMuscleMoments, jointMoment_vol) # moment_muscles_stim + result_cumulated_unknownMuscleMoments should equal (or close enough) jointMoment_vol
        color = 'red'
        if isClose: color='green'
        print("Is the result ({} + {}) close to our target jointMoment ({})? {}\n".format(
            colored('~'+str(np.around(moment_muscles_stim, 3)), color, attrs=['bold']),
            colored('~'+str(np.around(result_cumulated_unknownMuscleMoments, 3)), color, attrs=['bold']),
            colored('~'+str(np.around(jointMoment_vol, 6)), 'blue', attrs=['bold']),
            colored(isClose, color, attrs=['bold'])
            ))
        
        if not isClose:
            print("ERROR!")
            input()
        return res, Sstim, muscles_stim_osimNames, muscles_vol_osimNames

    # Else, we just return the Sstim scaling factor that will be true for all muscles (fmax) but all muscles will have excitation 0!
    else:
        return Sstim, muscles_stim_osimNames, muscles_vol_osimNames

def add_to_data(res, Sstim, muscles_stim_osimNames, muscles_vol_osimNames):
    
    ######
    # Tidy up : add the result into our initial dict to display with the names
    ######
    
    # Add the result of the known muscles (the ones tested with FES)
    for muscle in muscles_stim_osimNames:
        print(muscle)
        data_all[muscle]['weakness']        = Sstim
        data_all[muscle]['fmax_new']        = data_all[muscle]['fmax_original'] * Sstim
        data_all[muscle]['mc_new'] = Sstim
        data_all[muscle]['fes']             = True

    # Add the results of the unknown muscles
    i=0
    for muscle in muscles_vol_osimNames:
        print(muscle)
        data_all[muscle]['weakness'] = res.x[i]
        data_all[muscle]['fmax_new'] = data_all[muscle]['fmax_original'] * res.x[i]
        data_all[muscle]['mc_new'] = res.x[i]
        data_all[muscle]['fes'] = False
        i+=1

#####################################
# Main parameters, and batch running optimisations
#####################################

data_all = OrderedDict() # Initialize our dictionary that will receive the moment arms of the matching mucles

#####################################
# Ebow flexors

# The joint of interest as defined in opensim ('coordinates')
jointCoordinatesName = 'EL_x'

# The muscles and moments parameters as defined in opensim
muscles_all      = ["bic", "brachialis", "brachiorad"] # Each of these muscles are represented as several muscles in opensim. biceps is "bic" to get a partial match on "bic_l" and "bic_b_1/2"
muscles_vol      = ["brachialis", "brachiorad"]
muscles_stim     = ["bic"]
jointMoment_vol  = 8.75 # the joint moment under voluntary action only (all muscles)
jointMoment_stim = 8.75 # the joint moment under stimulated action (stimulated muscles only)
muscles_vol_osimNames  = [] # will be populated later, only used in final summary
muscles_stim_osimNames = []

# Compute moment arm and do the optimisation
title("Optimize elbow flexors", 'magenta')
res, Sstim, muscles_stim_osimNames, muscles_vol_osimNames = performOptimisation(jointCoordinatesName, muscles_all, muscles_vol, muscles_stim, jointMoment_vol, jointMoment_stim, muscles_vol_osimNames, muscles_stim_osimNames)
# Add the result data to the main dict
add_to_data(res, Sstim, muscles_stim_osimNames, muscles_vol_osimNames)

title("Finished optimisation for elbow flexors", 'magenta')

#####################################
# Ebow extensors

# special case for participant 1, the elbow extensors have a known activation of 0 because the jointMoment_vol is 0!
# so we only need to compute fmax based on fes, and distribute this to all the other muscles : no optimisation is required...

# get new list of muscles and moment, send to function again
muscles_all  = ['tric_long', 'tric_med', 'tric_lat', 'anconeus']
muscles_vol  = ['anconeus']
muscles_stim = ['tric_long', 'tric_med', 'tric_lat']
jointMoment_vol  = 0        # the joint moment under voluntary action only (all muscles)
jointMoment_stim = 3.75     # the joint moment under stimulated action (stimulated muscles only)
muscles_vol_osimNames  = [] # will be populated later, only used in final summary
muscles_stim_osimNames = []

# Compute moment arm and do NOT perform the optimisation : only use Sstim to compute the atrophy, and force activation ability to 0 as no voluntary movement
# Therefore assume that the atrophy factor is the same for all extensors, and activation is 0
title("Optimize elbow extensors", 'magenta')
Sstim, muscles_stim_osimNames, muscles_vol_osimNames = performOptimisation(jointCoordinatesName, muscles_all, muscles_vol, muscles_stim, jointMoment_vol, jointMoment_stim, muscles_vol_osimNames, muscles_stim_osimNames, optimise=False)
# Add the result data to the main dict : force activation to be 0

# Add the result of the known muscles (the ones tested with FES)
for muscle in muscles_stim_osimNames:
    print(muscle)
    data_all[muscle]['weakness'] = Sstim
    data_all[muscle]['fmax_new'] = data_all[muscle]['fmax_original'] * Sstim
    data_all[muscle]['mc_new']   = 0
    data_all[muscle]['fes']      = True

# Add the results of the unknown muscles
i=0
for muscle in muscles_vol_osimNames:
    print(muscle)
    data_all[muscle]['weakness'] = Sstim
    data_all[muscle]['fmax_new'] = data_all[muscle]['fmax_original'] * Sstim
    data_all[muscle]['mc_new']   = 0
    data_all[muscle]['fes']      = False
    i+=1

#####################################
# Wrist flexors and extensors

# Skipping for now - need to integrate wendy's model first

#####################################
# other muscles of the shoulder

# special case also : this time, we did not test the stimulation for these muscles.
### we can't distribute properly : if we separate the muscles into one specific motion, then we require values that are obviously too high to match the jointmoment.
### if we use the same muscle twice, then there is no good way to calculate the weakness factor as no jointMoment will be matched either!
### so restarting the shoulder muscles : just taking the average of the weakness factors and applying them to the shoulder muscles!
# NOT basing on innervation level as elbow flexion/extension is roughly C6 (with secondary c5,c6,c7,c8) which covers all of the shoulder...

def get_momentArm_forMuscles(muscles_list, data, verbose=False, returnOsimNames=False):
    # locate the muscles in the xml document
    for muscle in muscles_list:
    
        if verbose: print("\nLooking for muscle {}".format(muscle))
    
        for x in osim_muscle_data:
    
            name = x["@name"]

            # There can be several muscle 'fibres' to represent the same muscle so we do a partial match on the name
            string_match = re.search("^"+muscle, name) # regex to check if our target muscle is at beginning of muscle name

            if string_match:
    
                # When found a matching muscle, compute moment arm and save in the momentarm_dict
                if verbose: print("found matching name: {}".format(name))
    
                # get the muscle
                this_muscle = model.getMuscles().get(name)
                
                # compute moment arm
                this_moment_arm = abs(this_muscle.computeMomentArm(state, jointCoordinates))
                #print('moment arm is {}'.format(this_moment_arm))
                # keep it positive to have consistent graphs (otherwise elbow extensors' new Fmax will appear negative)
                # this_moment_arm = abs(this_moment_arm)
                if verbose: print('moment arm is {}'.format(this_moment_arm))
    
                # get the original Fmax
                fmax_original = float(x["max_isometric_force"])
    
                # compute fmax at this angle            
                this_muscle.setActivation(state, 1) # Set the activation to 1 (madesx)
                model.equilibrateMuscles(state)     # make sure states are in equilibrium
                active  = this_muscle.getActiveForceLengthMultiplier(state)
                passive = this_muscle.getPassiveForceMultiplier(state)
                this_maxisoforce = this_muscle.getMaxIsometricForce() # get max isometric force
                this_pennationAngle = this_muscle.getCosPennationAngle(state)
                fmax_thisAngle = this_maxisoforce * (active + passive) * this_pennationAngle
                if verbose: print("Fmax {} at angle {} is {}".format(this_maxisoforce, angle, fmax_thisAngle))
        
                # save to dict
                data[name] = {'moment_arm' : this_moment_arm, 'fmax_thisAngle' : fmax_thisAngle, 'fmax_original' : fmax_original}

# Display directly using pandas 
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 200)
df = pd.DataFrame(data_all).T
# What is the average weakness factor?
weakness_average_elbow = df['weakness'].mean()

# Scale all the shoulder muscles by this, then add to osim file and save it!
muscles_all  = ['delt_clav', 'coracobr', 'lat_dorsi', 'supra', 'pect_maj_t', 'pect_maj_c', 'subscap', 'teres_maj', 'teres_min', 'infra', 'delt_scap', 'trap_clav', 'lev_scap', 'pect_min', 'trap_scap', 'rhomboid', 'serr_ant']

# get moment arm to add to data_all:
# The joint of interest as defined in opensim ('coordinates')
jointCoordinatesName = 'GH_z'

jointCoordinates = model.updCoordinateSet().get(jointCoordinatesName) # define coordinates (ex: 'EL_x' for elbow flexion, 'GH_z' for shoulder elevation)
model.equilibrateMuscles(state) # make sure states are in equilibrium

get_momentArm_forMuscles(muscles_all, data_all, returnOsimNames=True)

# get osim names
muscles_all_osimNames = []
for muscle in muscles_all: # for each muscle that is stimulated
    
    for key, value in data_all.items(): # look into our data
        name = key
        string_match = re.search("^"+muscle, name) # regex to check if our target muscle is at beginning of muscle name
        if string_match: # if the name matches:
            # also save the name of each muscle as present in the osim file
            muscles_all_osimNames.append(name)


for muscle in muscles_all_osimNames:
    # set weakness, fmax_new, mc_new, fes=False
    #print(muscle)
    data_all[muscle]['weakness'] = weakness_average_elbow
    data_all[muscle]['fmax_new'] = data_all[muscle]['fmax_original'] * weakness_average_elbow
    data_all[muscle]['mc_new']   = weakness_average_elbow
    data_all[muscle]['fes']      = False

df2 = pd.DataFrame(data_all).T
print(df2)

######
# Add to osim file xml
######

for muscle in data_all:
    
    for x in osim_muscle_data:
        name = x["@name"]
        string_match = re.search("^"+muscle, name) # regex to check if our target muscle is at beginning of muscle name
        if string_match:
            # what is the original (scaled geometry) fmax in the current osim document?
            fmax_original = float(x["max_isometric_force"])

            # update values and create tags
            x["max_isometric_force"] = data_all[muscle]['fmax_new']
            x['min_control'] = 0
            x['max_control'] = data_all[muscle]['mc_new']

            x['weakness'] = data_all[muscle]['weakness']
            x['max_isometric_force_original_scaled'] = fmax_original
            x['muscle_has_fes_testing'] = data_all[muscle]['fes']

######
# LOCK the ROM :
######

# 2do : fix the xml parser

# flexion = 110  = GH_z +
# abduction = 90 = GH_z +
# extension = 45 = GH_z -
# adduction = 30 = GH_z -
# med rot   = 90 = GH_yy -
# lat rot   = 15 = GH_yy +

# --> using:
# GH_z  = -45 -> 110 (in rad = -0.7853981633974483 -> 1.9198621771937625)
# GH_y  = untouched
# GH_yy = -90 --> 15 (in rad = -1.5707963267948966 -> 0.2617993877991494)

# set clamped true
# save as new osim file

#osim_file_loaded["OpenSimDocument"]["Model"]["BodySet"]["objects"]['Body']['gh3']['CoordinateSet']['objects']['range']

osim_body_data = osim_file_loaded["OpenSimDocument"]["Model"]["BodySet"]["objects"]['Body'] # List of dicts (one dict per muscle)
print("\n@@@@@@@@@@@@\n")
for x in osim_body_data:
    print(x['@name'])

try:
    f = open(CURR_DIR+'/'+output_osim_file, "w")
    print(xmltodict.unparse(osim_file_loaded, pretty = True), file=f)
    f.close()
    #print("File saved as {}".format(CURR_DIR+'/'+output_osim_file))
    print("File saved successfully")
except:
    print("Error saving the file ({})".format(CURR_DIR+'/'+output_osim_file))



















