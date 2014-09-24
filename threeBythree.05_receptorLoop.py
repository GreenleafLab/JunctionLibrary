#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

# This function calls dependent functions to make the subset of
# final library that is all three by three junctions and less


##### IMPORT #####
import numpy as np
import os
import sys

# load custom libraries
import create_library
import globalvars
parameters = globalvars.Parameters()
from hjh.helix import Helix
from hjh.junction import Junction
import create_library


##### MODULE #####
"""
SETUP
set working directory, create filename to save all subsequent
sequences to.
Initialize count.
"""
#wd = os.path.join(os.getcwd(), 'libraries') # working directory
wd = parameters.wd

# check if working directory exists and if not, creates it
if not os.path.exists(wd):
    os.mkdir(wd)
    
# initialize file
filename = os.path.join(wd, 'all3x3junctions.05_receptor_loops.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 1


"""
ITERATING THROUGH KISSING LOOP/RECEPTORS
have many possible kissing loops, iterated through all of them
and change helix length as well for each
"""
junctionMotif = ('',)
receptorNames = parameters.receptorDict.keys()
loopNames     = parameters.loopDict.keys()
helixName   = 'rigid'
junction = Junction(junctionMotif)
totalLengths = [8, 9, 10, 11, 12]
helices = Helix(parameters.helixDict[helixName], junction.length).formatHelices(totalLengths, [5]*len(totalLengths))

for receptorName in receptorNames:
    for loopName in loopNames:
        count = create_library.saveSet(junction, helices, helixName, receptorName, loopName, f, logfile, count)

# close
f.close()
logfile.close()
