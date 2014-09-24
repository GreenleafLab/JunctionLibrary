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
filename = os.path.join(wd, 'all3x3junctions.03_centralregion.fasta')
print 'saving to %s'%filename
f = open(filename, 'w')

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

junctionMotif = ('B1', 'M')

junction = Junction(junctionMotif)
    
# take subset of junctions spanning 64 junctions. 
if junction.howManyPossibilities() > cutOffNumber:
    subsetIndex = np.around(np.linspace(0, junction.howManyPossibilities()-1, cutOffNumber)).astype(int)
    junction.sequences = junction.sequences[subsetIndex]

junctionSequence = junction.sequences[0]
helixName = 'rigid'
# now save

helices = Helix(parameters.helixDict[helixName], junction.length).centralRegion()
for helixSequence in helices:
    names = [helixName, receptorName, loopName]
    name = '%s.%s.%d_%d'%('.'.join(names),
              '_'.join(junctionSequence),
                        len(helixSequence['h1_side1']),
                        len(helixSequence['h2_side1']))
    f.write('>%s\n'%name)
    sequenceList = create_library.twoWayJunctionOrder(helixSequence,
                                                      junctionSequence,
                                                      parameters.receptorDict[receptorName],
                                                      parameters.tectoBase,
                                                      parameters.sequencingAdapters)
    sequence = parameters.loopDict[loopName]
    for outside in sequenceList[:-1]:
        sequence = create_library.threadTogether(sequence, outside)
    f.write('%s\n'%sequence)
 

# close
f.close()

