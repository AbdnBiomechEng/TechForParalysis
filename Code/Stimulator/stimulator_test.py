# -*- coding: utf-8 -*-
"""
https://github.com/humancomputerintegration/rehamove-integration-lib

Notes:
    - does not run on university computer that restricts access to USB/COM ports from scripts
    - the pre-compiled python library for windows that is available on github requires python 3.7 (otherwise recompile)
    - keep using the low-level mode, to send individual pulses
"""
# Debug: clear variables of IDE
from IPython import get_ipython
try: __IPYTHON__
except: pass
else: get_ipython().run_line_magic('reset','-sf')

# imports (load pre-compiled windows library in python3.7)
import sys, os
dirpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dirpath+'\\libs')
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


