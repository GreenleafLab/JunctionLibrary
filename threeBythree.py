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

##### FUNCTIONS #####
def getFilename(helixName, description, loopName, receptorName, junctionMotifs):
    fileName = [helixName,
                description,
                loopName,
                receptorName,
                '_'.join([''.join(result) for result in junctionMotifs])]
    return '.'.join(fileName)+'.txt'

def saveSet(junction, helices, helixName, receptorName, loopName, f, countAll):
    """
    this is just a space saver. Calls the saveSequences command for a given
    set of junctions, helices,  receptor, and loop)
    """
    f, count = create_library.saveSequences([helixName, receptorName, loopName],
                         junction.sequences,
                         helices,
                         parameters.receptorDict[receptorName],
                         parameters.loopDict[loopName],
                         parameters.tectoBase,
                         parameters.sequencingAdapters,
                         f)
    print 'Saved %d sequences...'%(count+countAll)
    logfile.write('%s\t%s\t%s\t%s\t%d\t%d\t%d\t\t%d\n'%(receptorName, loopName, helixName,
                                                      ','.join(junction.motif),
                                          junction.howManyPossibilities(),
                                          len(junction.sequences),
                                          len(helices),
                                          count
                                          ))
    return f, count

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
filename = os.path.join(wd, 'all3x3junctions.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 1


"""
STANDARD: One position, Rigid and Watson Crick
Save all junctions in one position in two different helix contextx. 
"""
junctionMotifs = parameters.allJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixNames   = parameters.standardHelixNames
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    for helixName in helixNames:
        # helices in default location
        helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
        f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)

"""
DIFFERENT HELIX CONTEXT: One position, subset of junctions, ten other helix contexts
Save subset of junctions in one position in ten different helix contextx. 
"""
junctionMotifs = parameters.allJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixNames   = parameters.otherHelixNames
cutOffNumber = 12
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    
    # take a subset of junctions of each junction topology
    if junction.howManyPossibilities() > cutOffNumber:
        subsetIndex = np.around(np.linspace(0, junction.howManyPossibilities()-1, cutOffNumber)).astype(int)
        junction.sequences = junction.sequences[subsetIndex]
    
    # for each helix name
    for helixName in helixNames:
        # helices in default location
        helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
        f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)

"""
JUNCTIONS IN ONE POSITION DIFFERENT LOOP
Save subset of junctions with different loop.
"""
junctionMotifs = parameters.allJunctions
receptorName = 'R1'
loopName     = 'badLoop'
helixName    = parameters.standardHelixNames
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    # for each helix name
    for helixName in helixNames:
        helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
        f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)

"""
JUNCTIONS IN ALL 20 POSITIONS
Save subset of junctions located in 'along' set
"""
# save one helix context, many junctions, in many different locations
print 'Doing 20 different positions of subset of junctions'
junctionMotifs = parameters.alongJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixName    = 'rigid'
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    # if junction is 'W', do it in all different helix contexts. Everythin else, do in two
    if junctionMotif = 'W':
        helixNames = parameters.allHelixNames
        for helixName in helixNames:

    else:
        helixNames = parameters.standardHelixNames
    helices = Helix(parameters.helixDict[helixName], junction.length).alongHelix()
    f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)
"""
JUNCTIONS IN CENTRAL REGION
Save subset of junctions located in 'central' set
"""
# save one helix context, many junctions, in many different locations
print 'Doing 13 different positions of subset of junctions'
junctionMotifs = parameters.centralJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixName    = 'rigid'
cutOffNumber = 576
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    
    # take subset of junctions spanning 576 junctions. Should only affect the 3x3 loop
    if junction.howManyPossibilities() > cutOffNumber:
        subsetIndex = np.around(np.linspace(0, junction.howManyPossibilities()-1, cutOffNumber)).astype(int)
        junction.sequences = junction.sequences[subsetIndex]
        
    # now save
    helices = Helix(parameters.helixDict[helixName], junction.length).centralRegion()
    f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)

# close
f.close()
logfile.close()
