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
import create_library

# load stored sequences from globalvars file
import globalvars
parameters = globalvars.Parameters()

##### FUNCTIONS #####
def getFilename(helixName, description, loopName, receptorName, junctionMotifs):
    fileName = [helixName,
                description,
                loopName,
                receptorName,
                '_'.join([''.join(result) for result in junctionMotifs])]
    return '.'.join(fileName)+'.txt'

def saveSet(junction, helices, receptorName, loopName, f, countAll):
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
    return f, count

def saveSetNoJunction(junction, helixName, receptorName, loopName, f, countAll):
    """
    given a class of junctions, helixName, receptorName, and loopName,
    find helices associated with NO junction (1 helix).
    save all possible to open file f.
    """
    helices = Helix(parameters.helixDict[helixName], junction.length).noJunctionLocation()
    f, count = saveSet(junction, helices, receptorName, loopName, f, countAll)
    return f, count+countAll

def saveSetAlong(junction, helixName, receptorName, loopName, f, countAll):
    """
    given a class of junctions, helixName, receptorName, and loopName,
    find helices associated with junctions positioned in 20 DIFFERENT locations (20 helices).
    save all possible to open file f.
    """
    helices = Helix(parameters.helixDict[helixName], junction.length).alongHelix()
    f, count = saveSet(junction, helices, receptorName, loopName, f, countAll)
    return f, count+countAll

def saveSetCentral(junction, helixName, receptorName, loopName, f, countAll):
    """
    given a class of junctions, helixName, receptorName, and loopName,
    find helices associated with junctions positioned in 
    """
    helices = Helix(parameters.helixDict[helixName], junction.length).centralRegion()
    f, count = saveSet(junction, helices, receptorName, loopName, f, countAll)
    return f, count+countAll


##### MODULE #####

"""
First step: do all junctions that are going to be in 'along' configuration of helix
"""
# open filename
wd = os.path.join(os.getcwd(), 'libraries')
if not os.path.exists(wd): os.mkdir(wd)

filename = os.path.join(wd, 'all3x3junctions.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

count = 1


"""
order of operations for no junctions
"""
print 'Doing no Junction set'
junctionMotif = ('',)
junction = Junction(junctionMotif)

# saving all helixes, no junction, correct receptors..
print '\tsaving all helixes, no junction, correct receptors..'
receptorName = 'R1'
loopName     = 'GGAA'
helixNames   = parameters.helixDict.keys()  # all helices
for helixName in helixNames:
    f, count = saveSetNoJunction(junction, helixName, receptorName, loopName, f, count)

# saving only one helix, no junction, different loop   
print '\tsaving only one helix, no junction, different loop ..'
receptorName = 'R1'
loopName     = 'GAAA'
helixName    = 'rigid'
f, count = saveSetNoJunction(junction, helixName, receptorName, loopName, f, count)

# saving only one helix, no junction, 6 different receptors
print '\tsaving only one helix, no junction, 6 different receptors'
receptorNames = ['KL1', 'KL2']
loopName      = 'GGAA'
helixName     = 'rigid'
for receptorName in receptorNames:
    f, count = saveSetNoJunction(junction, helixName, receptorName, loopName, f, count)

"""
order of operations for subset of junctions located in 'along' set
"""
# save one helix context, many junctions, in many different locations
print 'Doing 20 different positions of subset of junctions'
junctionMotifs = [('M',), ('M', 'M'), ('B1',), ('B2',), ('W',)]
receptorName = 'R1'
loopName     = 'GGAA'
helixName    = 'rigid'
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    f, count = saveSetAlong(junction, helixName, receptorName, loopName, f, count)

# save one helix context, many junctions, in many different locations
print 'Doing 13 different positions of subset of junctions'
junctionMotifs = [('M',), ('M', 'M'), ('B1',), ('B2',), ('W',)]
receptorName = 'R1'
loopName     = 'GGAA'
helixName    = 'rigid'

print 'Doing 1 different positions of subset of junctions in one location, all different helices'
