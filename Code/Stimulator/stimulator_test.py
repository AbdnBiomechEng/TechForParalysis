# -*- coding: utf-8 -*-
"""
https://github.com/humancomputerintegration/rehamove-integration-lib

Notes:
    - does not run on university computer that restricts access to USB/COM ports from scripts
    - the pre-compiled python library for windows that is available on github requires python 3.7 (otherwise recompile)
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

# create instance
r = rehamove.Rehamove("COM5")

# get info
r.info()

# send pulse
r.pulse("blue", 30, 50)

