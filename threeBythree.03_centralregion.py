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
wd = os.path.join(os.getcwd(), 'libraries') # working directory

# check if working directory exists and if not, creates it
if not os.path.exists(wd):
    os.mkdir(wd)
    
# initialize file
filename = os.path.join(wd, 'all3x3junctions.03_centralregion.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 1

"""
JUNCTIONS IN CENTRAL REGION
Save subset of junctions located in 'central' set
"""
# save one helix context, many junctions, in many different locations
print 'Doing 13 different positions of subset of junctions'
junctionMotifs = parameters.sameLoopJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixNames    = parameters.standardHelixNames
cutOffNumber = 64
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    
    # take subset of junctions spanning 64 junctions. 
    if junction.howManyPossibilities() > cutOffNumber:
        subsetIndex = np.around(np.linspace(0, junction.howManyPossibilities()-1, cutOffNumber)).astype(int)
        junction.sequences = junction.sequences[subsetIndex]

    # now save
    for helixName in helixNames:
        helices = Helix(parameters.helixDict[helixName], junction.length).centralRegion()
        count = create_library.saveSet(junction, helices, helixName, receptorName, loopName, f, logfile, count)
    
"""
JUNCTIONS IN CENTRAL REGION WITH DIFFERENT LOOP
Save subset of jucntions located in 'central' set with different loop
"""
junctionMotifs = parameters.differentLoopJunctions
receptorName = 'R1'
loopNames     = ['goodLoop', 'badLoop']
helixNames    = parameters.standardHelixNames
cutOffNumber = 64
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    
    # take subset of junctions spanning 64 junctions. 
    if junction.howManyPossibilities() > cutOffNumber:
        subsetIndex = np.around(np.linspace(0, junction.howManyPossibilities()-1, cutOffNumber)).astype(int)
        junction.sequences = junction.sequences[subsetIndex]
    
    # do two different helices
    for helixName in helixNames:
    
        # do two different loops
        for loopName in loopNames:
            # now save
            helices = Helix(parameters.helixDict[helixName], junction.length).centralRegion()
            count = create_library.saveSet(junction, helices, helixName, receptorName, loopName, f, logfile, count)

# close
f.close()
logfile.close()
