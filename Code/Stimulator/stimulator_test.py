# -*- coding: utf-8 -*-
"""
The python wrapper is coming from: https://github.com/humancomputerintegration/rehamove-integration-lib
    # please read the license and limitations

Notes:
    - had to disable USB sleep in power settings on a windows laptop
    - one of the usb port is not working reliably (need to use the left port on ABD's laptop)
    - the pre-compiled python library for windows that is available on github requires python 3.7 (otherwise recompile)
    - keep using the low-level mode, to send individual pulses

2do:
    - make anaconda env including serial-tool
    - include test to load the lib. windows/linux
    - add a loop to wait for the hardware to be connected

"""

# Debug: clear variables of IDE
from IPython import get_ipython
try: __IPYTHON__
except: pass
else: get_ipython().run_line_magic('reset','-sf')

# imports (load pre-compiled windows library in python3.7)
import sys, os
dirpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dirpath+'\\windows\\libs') # windows
#sys.path.append(dirpath+'/libs/linux') # linux (need to give permission in /dev/ttyUSB*)
import rehamove

# Detect which COM port is in use
comPort = ''
from serial.tools import list_ports
ports = list_ports.comports()
for port, desc, hwid in sorted(ports):
        if 'PID=0403:6014' in hwid: # this is the ID of our rehamove
            # if it's a match, we save the COM port number
            comPort = port

# create instance
r = rehamove.Rehamove(comPort)

# get info
r.info()
r.battery()
r.version()

# set to low-level mode to control each individual pulse
r.change_mode(0)

# send pulse
r.pulse("blue", 30, 50) # cable_color, milliAmps, duration in microseconds



