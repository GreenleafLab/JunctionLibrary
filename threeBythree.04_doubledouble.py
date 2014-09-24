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
filename = os.path.join(wd, 'all3x3junctions.04_double_double.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

# initialize counts
count = 1

receptorName = 'R1'
loopName     = 'goodLoop'
helixName    = 'wc'

# junctions of length 0
numConfigurations = 4
anybases = ['A', 'G', 'T', 'C']
sequences = np.empty(numConfigurations*len(anybases), dtype ={'names':('j1_side1', 'j2_side1', 'j2_side2', 'j1_side2'),
    'formats':('S100', 'S100', 'S100', 'S100')})
junction_length = 0
indx = 0
for i in range(numConfigurations):
    for anybase in anybases:
        if i==0:
            sequences[indx] = ('A', '', anybase, '')
        elif i==1:
            sequences[indx] = ('A', anybase, '', '')
        elif i==2:
            sequences[indx] = (anybase, 'A', '', '')
        elif i==3:
            sequences[indx] = (anybase, '', 'A', '')
        print '%s\t%s\t%s\t%s'%(sequences[indx]['j1_side1'], sequences[indx]['j2_side1'], sequences[indx]['j2_side2'], sequences[indx]['j1_side2'])
        indx += 1
helices = Helix(parameters.helixDict[helixName], junction_length).doubleDouble()
count = create_library.saveSequenceDoubleDouble([helixName, receptorName, loopName],
                             sequences,
                             helices,
                             parameters.receptorDict[receptorName],
                             parameters.loopDict[loopName],
                             parameters.tectoBase,
                             parameters.sequencingAdapters,
                             f)

# now do single mismaches, keeping two of the bases G's to reduce computational space
indx = 0
anybases = ['A', 'G', 'T']
sequences = np.empty(numConfigurations*len(anybases), dtype ={'names':('j1_side1', 'j2_side1', 'j2_side2', 'j1_side2'),
    'formats':('S100', 'S100', 'S100', 'S100')})
for i in range(numConfigurations):
    for anybase in anybases:
        if i==0:
            sequences[indx] = ('G', anybase, 'G', anybase)
        elif i==1:
            sequences[indx] = (anybase, anybase, 'G', 'G')
        elif i==2:
            sequences[indx] = ('G', 'G', anybase, anybase)
        elif i==3:
            sequences[indx] = (anybase, 'G', anybase, 'G')
        print '%s\t%s\t%s\t%s'%(sequences[indx]['j1_side1'], sequences[indx]['j2_side1'], sequences[indx]['j2_side2'], sequences[indx]['j1_side2'])
        indx += 1
helices = Helix(parameters.helixDict[helixName], junction_length).doubleDouble()
count = create_library.saveSequenceDoubleDouble([helixName, receptorName, loopName],
                             sequences,
                             helices,
                             parameters.receptorDict[receptorName],
                             parameters.loopDict[loopName],
                             parameters.tectoBase,
                             parameters.sequencingAdapters,
                             f)

# close
f.close()
