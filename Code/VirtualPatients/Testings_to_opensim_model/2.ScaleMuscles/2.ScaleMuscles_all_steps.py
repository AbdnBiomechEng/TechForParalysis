#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script:
    - lists main elbow flexor muscles from osim file,
    - sets the model in position (arm flexed 90°)
    - calculates the moment arms and max forces in this position
    - calls matlab's optimizer fmincon to compute newFmax and newMaxControl
    - edits the osim file with the new values and save as new model file

2do:
    [o] scale parameter for final fmax should be based on geometry (opensim style)
    [o] compute max force of a muscle at specific angle using the force length curve directly, because opensim has variable support depending on the muscle model used (shuttle 1993 is deprecated)
"""

#####################################
# Imports and setting up
#####################################

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
import numpy as np
from collections import OrderedDict

CURR_DIR = os.path.dirname(os.path.realpath(__file__))
os.chdir(os.path.abspath(os.path.dirname(__file__)))   # change current working directory

#####################################
#
# Step 1 : compute the moment arm and max force of the target muscles from the opensim file, with the elbow in 90° flexion and in neutral prono-supination
# Each muscle is defined as several muscles in the opensim model
#
#####################################

import opensim
import opensim as osim
import xmltodict

print("Step 1")

#####################################
# ------ Load the model
#####################################

input_osim_file  = 'Data/Models/das4_v4_original.osim'
#input_osim_file  = 'Data/Models/1.model_scaled_bones.osim'
#input_osim_file  = 'Data/Models/model-das3-MOD2-noOffsetBestPractice-noRotation-newVisualization.osim'
model = opensim.Model(input_osim_file)
visualizeFlag = True
model.setUseVisualizer(visualizeFlag)

#####################################
# ------ Initialize the model and create the states
#####################################

state = model.initSystem()

#####################################
# ------ Set the elbow flexion to 90° and supination to maximum
#####################################

elbow = model.getJointSet().get("hu") # Get the elbow joint
angle = 90
elbow.getCoordinate().setValue(state, np.deg2rad(angle)) # 0.5 * osim.SimTK_PI)

ulna_radius = model.getJointSet().get("ur") # Get the ulna/radius joint # coordinates : PS_y
ulna_radius.getCoordinate().setValue(state, 0) # Set the coordinates to max supination (values go from 0 to 3 so set to 0)

#####################################
# ------ Configure the visualizer for visual confirmation of position
#####################################

if visualizeFlag: 
    viz = model.updVisualizer().updSimbodyVisualizer()
    viz.setBackgroundColor(osim.Vec3(0)) # white
    viz.setGroundHeight(-2)
    model.getVisualizer().show(state)

#####################################
# ------ Get moment arm of all the muscles required for the elbow flexion
#####################################

elbowFlexCoord = model.updCoordinateSet().get('EL_x') # define coordinates

model.equilibrateMuscles(state) # make sure states are in equilibrium

muscles_target = ["bic", "brachialis", "brachiorad"] # Each of these muscles are represented as several muscles in opensim. biceps is "bic" to get a partial match on "bic_l" and "bic_b_1/2"

data = OrderedDict() # Initialize our dictionary that will receive the moment arms of the matching mucles

#####################################
# Load the original osim file in external xml parser
#####################################

myxml            = open(CURR_DIR + "/" + input_osim_file).read()
osim_file_loaded = xmltodict.parse(myxml)
osim_muscle_data = osim_file_loaded["OpenSimDocument"]["Model"]["ForceSet"]["objects"]["Schutte1993Muscle_Deprecated"] # List of dicts (one dict per muscle)

# locate the muscles in the xml document
for muscle in muscles_target:

    print("\nLooking for muscle {}".format(muscle))

    for x in osim_muscle_data:

        name = x["@name"]

        # There can be several muscle 'fibres' to represent the same muscle so we do a partial match on the name
        if muscle in name:

            # When found a matching muscle, compute moment arm and save in the momentarm_dict
            print("found matching name: {}".format(name))

            # get the muscle
            this_muscle = model.getMuscles().get(name)
            
            # compute moment arm
            this_moment_arm = this_muscle.computeMomentArm(state, elbowFlexCoord)
            print('moment arm is {}'.format(this_moment_arm))

            # get the original Fmax
            fmax_original = float(x["max_isometric_force"])

            # compute fmax at this angle            
            this_muscle.setActivation(state, 1) # Set the activation to 1 (madesx)
            model.equilibrateMuscles(state) # make sure states are in equilibrium
            active  = this_muscle.getActiveForceLengthMultiplier(state)
            passive = this_muscle.getPassiveForceMultiplier(state)
            this_maxisoforce = this_muscle.getMaxIsometricForce() # get max isometric force
            this_pennationAngle = this_muscle.getCosPennationAngle(state)
            fmax_thisAngle = this_maxisoforce * (active + passive) * this_pennationAngle
            print("Fmax {} at angle {} is {}".format(this_maxisoforce, angle, fmax_thisAngle))
    
            # save to dict
            data[name] = {'moment_arm' : this_moment_arm, 'fmax_thisAngle' : fmax_thisAngle, 'fmax_original' : fmax_original}

'''
# Result example: data =

{'bic_l'       : {'moment_arm': 0.05143613931184783,  'fmax_thisAngle': 332.4, 'fmax_original' : 347.0},
 'bic_b_1'     : {...},
 'bic_b_2'     : {...},
 'brachialis_1': {...},
 'brachialis_2': {...},
 'brachialis_3': {...},
 'brachialis_4': {...},
 'brachialis_5': {...},
 'brachialis_6': {...},
 'brachialis_7': {...},
 'brachiorad_1': {...},
 'brachiorad_2': {...},
 'brachiorad_3': {...},
'''

# print("\nEnd of Step 1 : press key")
# Pause
# input()

#####################################
#
# Step 2 : compute newFmax and control_max of the elbow flexor muscles, using optimizer (matlab fmincon)
#
# From : https://github.com/AbdnBiomechEng/TechForParalysis/blob/muscle_params_optimisation/Code/Model/CustomiseRoutines/MuscleParameterOptimisation.m
#
#####################################

import matlab.engine

print("\nStep 2")

'''
    Commments from the matlab file:
        The elbow is actively flexed by the (1) biceps, (2) brachialis and (3) brachioradialis muscles.
        The resultant moment is measured and we want to scale the muscle forces to the measured data.
        Moment arms of the muscles are below
'''

#####################################
# -------- Data
#####################################

stimulated_muscle_numbers = [0,1,2] # the muscles tested are biceps (1,2,3) = 0,1,2 in index python
M_vol                     = 8.75 # the joint moment under voluntary action only
M_stim                    = 8.75 # the joint moment under stimulated action

# Make lists of the moment_arm and Fmax values
r_list              = []
Fmax_thisAngle_list = []
Fmax_original_list  = []
for x in data:
    r_list.append(data[x]['moment_arm'])
    Fmax_thisAngle_list.append(data[x]['fmax_thisAngle'])
    Fmax_original_list.append(data[x]['fmax_original'])

# list to numpy array
r_mus_python  = np.array(r_list)
Fmax_thisAngle_python   = np.array(Fmax_thisAngle_list)

# Convert data from numpy to matlab format so can send to optimizer
r_mus_matlab  = matlab.double(np.array(r_mus_python).tolist())
Fmax_thisAngle_matlab   = matlab.double(np.array(Fmax_thisAngle_python).tolist())

#####################################
# -------- Parameters
#####################################

# How many values (muscles) are we dealing with?
number_muscles = len(r_mus_python)

# Initial guess, or starting point
x0_python = np.zeros([number_muscles,]) # zeros
# Linear equality constraints
lec1_python = np.array([])
lec2_python = np.array([])
# A and B in A[x] = B, the linear equality contraint
A_python = np.array([])
B_python = np.array([])
# Lower bound
lb_python = np.zeros([number_muscles,]) # zeros
# Upper bound
ub_python = np.ones([number_muscles,]) # ones

# Convert parameters to matlab
x0_matlab   = matlab.double(np.array(x0_python).tolist())
lec1_matlab = matlab.double(np.array(lec1_python).tolist())
lec2_matlab = matlab.double(np.array(lec2_python).tolist())
A_matlab    = matlab.double(np.array(A_python).tolist())
B_matlab    = matlab.double(np.array(B_python).tolist())
lb_matlab   = matlab.double(np.array(lb_python).tolist())
ub_matlab   = matlab.double(np.array(ub_python).tolist())

#####################################
# -------- Run matlab's fmincon
#####################################

print('--- starting matlab engine ---')

eng = matlab.engine.start_matlab()

print('--- calling fmincon ---')

# ---- set options ----
'''
options = {'Algorithm': 'active-set', 'AlwaysHonorConstraints': 'bounds',
    'display': 'iter-detailed', 'MaxIter': 1000, 'MaxFunEvals': 10000,
    'TolCon': 1e-6, 'TolFun': 1e-6, 'Diagnostics': 'on'}
'''

mc_matlab, J, exitflag, output = eng.optimize(M_vol, M_stim, r_mus_matlab, Fmax_thisAngle_matlab, x0_matlab, lb_matlab, ub_matlab, nargout=4)

# close matlab engine
print('--- closing matlab engine ---')
eng.quit() 

#####################################
# -------- Get scale, max_control and NewFmax
#####################################

# ---- max control ----
# already have "mc" from the optimizer
# convert back to python format
mc_python = np.array(mc_matlab._data.tolist())
mc_python = mc_python.reshape(mc_matlab.size).transpose()

# ---- Fmax ----
# reshape from opensim
Fmax_thisAngle_python = Fmax_thisAngle_python.reshape(number_muscles,1)

# ---- scale ----

sum_Moments_stimulatedMuscles = 0 # the moments of the biceps that was stimulated under FES

for i in stimulated_muscle_numbers:
    sum_Moments_stimulatedMuscles = sum_Moments_stimulatedMuscles + (r_mus_python[i] * mc_python[i] * Fmax_thisAngle_python[i])

scale = 2*M_stim / (sum_Moments_stimulatedMuscles)

# ---- newFmax ----
newFmax = scale * mc_python * Fmax_thisAngle_python

# Make a new dict from the output of the optimizer to clean things up
optimizer_output = OrderedDict()
i = 0
for name in data:
    optimizer_output[name] = { 'newFmax' :  float(newFmax[i]) , 'newMax_control' : float(mc_python[i]) }    
    i+=1

# Plot this https://plotly.com/python/horizontal-bar-charts/#colored-horizontal-bar-chart
import plotly.graph_objects as go
import plotly.express as px
import plotly

fig = go.Figure()

# Muscle names and values to plot
plot_list_names    = []
plot_original_fmax = []
plot_newFmax       = []
for x in data:
    plot_list_names.append(x)
    plot_original_fmax.append(data[x]['fmax_original'])
    plot_newFmax.append(optimizer_output[x]['newFmax'])

# Add the original Fmax
fig.add_trace(go.Bar(
    y=plot_list_names,
    x=plot_original_fmax,
    name='Original Fmax',
    orientation='h',
    marker=dict(
        color='rgba(246, 78, 139, 0.6)',
        line=dict(color='rgba(246, 78, 139, 1.0)', width=3)
    )))

# Add the newFmax
fig.add_trace(go.Bar(
    y=plot_list_names,
    x=plot_newFmax,
    name='new Fmax',
    orientation='h',
    marker=dict(
        color='rgba(58, 71, 80, 0.6)',
        line=dict(color='rgba(58, 71, 80, 1.0)', width=3)
    )))

#######
# Customise the plot
fig.update_layout(
                  showlegend=True,
                  title="Muscle Force (N)",
                  barmode='overlay', # stack https://plotly.com/python-api-reference/generated/plotly.express.bar
                  font=dict(
                      #family="Courier New, monospace",
                      size=32,
                      # color="RebeccaPurple"
                      ),
                  )
########
# Show fig
config_plot = {'displayModeBar': False}
fig.show(renderer="png",config=config_plot)
plotly.offline.plot(fig, config=config_plot, filename="/tmp/temp_plot.html")

#####################################
#
# Step 3 : save in the opensim file
#
#####################################

'''
This section:
    - opens the opensim model file of our scaled (skeleton) subject
    - scales all the muscle force by the computed scale factor
    - updates values of the computed Fmax and max_control tags for the elbow flexor muscles
'''

#####################################
# -------- Load the original osim file
#####################################

#input_osim_file = "Data/Models/1.model_scaled_bones.osim"
output_osim_file = "Data/Models/2.model_scaled_muscles.osim"

myxml = open(CURR_DIR + "/" + input_osim_file).read()
osim_file_loaded = xmltodict.parse(myxml)
osim_muscle_data = osim_file_loaded["OpenSimDocument"]["Model"]["ForceSet"]["objects"]["Schutte1993Muscle_Deprecated"]  # List of dicts (one dict per muscle)

#####################################
# -------- Scale all muscle forces by 'scale' factor
#####################################

for x in osim_muscle_data:
    name = x["@name"]
    # print("Scaling muscle {}".format(name))

    # Scale Fmax
    fmax_original = float(x["max_isometric_force"])
    x["max_isometric_force_original_notScaled"] = fmax_original

    # Add some tags to ease reading later on
    x["max_isometric_force"]    = fmax_original * scale
    x["min_control"]            = 0
    x["max_control"]            = 1
    x["muscle_has_fes_testing"] = 0

#####################################
# -------- Edit our target muscles
#####################################

# locate them in the xml document
for muscle in optimizer_output:

    print("doing muscle {}".format(muscle))

    for x in osim_muscle_data:

        name = x["@name"]

        # When found, edit the muscle parameters of the osim file
        # There can be several muscle fibres to represent the same muscle so we do a partial match on the name
        if muscle in name:

            print("found matching name: {}".format(name))
            
            # edit fmax
            print("setting old value {} to new value {}".format(x['max_isometric_force'][0] , optimizer_output[muscle]['newFmax']))
            x['max_isometric_force'] = optimizer_output[muscle]['newFmax']
            # edit max_control
            x["max_control"] = optimizer_output[muscle]['newMax_control']

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
    

