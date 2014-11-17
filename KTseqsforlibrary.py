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
import globalvars_KT
parameters = globalvars_KT.Parameters()
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
    logfile.write('%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n'%(helixName, receptorName, loopName,
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
filename = os.path.join(wd, 'KinkTurnSequences_9_29_14.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 1


#"""
#JUNCTIONS IN ONE POSITION MANY HELICES
#Save all junctions with 10 different helix contexts. 
#"""
#junctionMotifs = parameters.allJunctions
#receptorName = 'R1'
#loopName     = 'goodLoop'
#helixNames   = parameters.allHelixNames # all helices
#for junctionMotif in junctionMotifs:
#    junction = Junction(junctionMotif)
#    for helixName in helixNames:
#        # helices in default location
#        helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
#        #f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, countAll)
#        f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)
#
#"""
#JUNCTIONS IN ONE POSITION DIFFERENT TERTIARY CONTACT
#Save subset of junctIons with different loop.
#"""
#junctionMotifs = [('',)]
#receptorName = 'R1'
#loopName     = 'badLoop'
#helixName    = 'rigid'
#for junctionMotif in junctionMotifs:
#    helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
#    f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)
#
#"""
#JUNCTIONS IN ONE POSITION DIFFERENT TERTIARY CONTACT
#Save subset of junctions with different receptors
#"""
#junctionMotifs = [('',)]
#receptorNames = ['KL1', 'KL2']
#loopName      = 'goodLoop'
#helixName     = 'rigid'
#for junctionMotif in junctionMotifs:
#    for receptorName in receptorNames:
#        helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
#        f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)
"""
6 DIFFERENT KINK TURNS  IN 5 positions 
Save subset of junctions located in 'along' set
"""
# save one helix context, many junctions, in many different locations
count = 1
print 'Doing 5 different positions for 6 different kinkturns'
KTMotifs = parameters.KinkTurnsdict
KTnames = parameters.allKT_Names
receptorName = 'R1'
loopName     = 'goodLoop'
helixName    = 'KThelix1'
junctionMotifs = [('',)]


junction = Junction(junctionMotifs[0])

for KTname in KTnames:
    countAll = count + 1
    KTseq = KTMotifs[KTname]
    junction.sequences[0] = KTseq
    junction.length = 6
    #helices = Helix(parameters.KThelixDict[helixName],0).alongHelix()
    helices = Helix(parameters.helixDict['rigid'], junction.length).alongHelix()
    f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)
   
"""
JUNCTIONS IN CENTRAL REGION
Save subset of junctions located in 'central' set
"""
# save one helix context, one KT sequence, with different helix lengths -2,-1,0,+1,+2
count = 1
print 'Doing  different helix lengths for KT_1'
KTMotifs = parameters.KinkTurnsdict
KTnames = parameters.allKT_Names
receptorName = 'R1'
loopName     = 'goodLoop'
helixName    = 'KThelix1'
KTmotif_base = parameters.allKT_Names[1]
junctionMotif = [('',)]

junction = Junction(junctionMotifs[0])
    
KTseq = KTMotifs[KTmotif_base]  
junction.sequences[0] = KTseq
junction.length = 6
        
# now save
#helices = Helix(parameters.KThelixDict[helixName], 0).centralRegion()
helices = Helix(parameters.helixDict['rigid'], junction.length).centralRegion()
f, count = saveSet(junction, helices, helixName, receptorName, loopName, f, count)

"""
Iterate through all posible mismatches, WC bases, and bulges for the three regions of the kink turn
"""
junctionMotifs = parameters.KTjunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helix_KTs = parameters.KThelixDict;
for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)    
    if junction.motif == ('M', 'M'):
        helices = Helix(parameters.helixDict['rigid'], junction.length).KThelixMM()
        helixname = 'rigid'
        f, count = saveSet(junction, helices, helixname, receptorName, loopName, f, count)
        
    if junction.motif == ('W', 'W'):
        helices = Helix(parameters.helixDict['rigid'], junction.length).KThelixWC()
        helixname = 'rigid'
        f, count = saveSet(junction, helices, helixname, receptorName, loopName, f, count)
        
    if junction.motif == ('B2', 'B2', 'B2'):
        helices = Helix(parameters.helixDict['rigid'], junction.length).KThelixBulge()
        helixname = 'rigid'
        f, count = saveSet(junction, helices, helixname, receptorName, loopName, f, count)




# close
f.close()
logfile.close()
