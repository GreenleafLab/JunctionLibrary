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
import pandas as pd

# load custom libraries
import create_library_dataframe
import globalvars_KT
parameters = globalvars_KT.Parameters()
from hjh.helix import Helix
from hjh.junction import Junction


##### FUNCTIONS #####
def getFilename(helixName, description, loopName, receptorName, junctionMotifs):
    fileName = [helixName,
                description,
                loopName,
                receptorName,
                '_'.join([''.join(result) for result in junctionMotifs])]
    return '.'.join(fileName)+'.txt'

def saveSet(junction, helices, helixName, receptorName, loopName, f,f2, countAll):
    """
    this is just a space saver. Calls the saveSequences command for a given
    set of junctions, helices,  receptor, and loop)
    """
    count = create_library_dataframe.saveSequences([helixName, receptorName, loopName],
                         junction.sequences,
                         helices,
                         parameters.receptorDict[receptorName],
                         parameters.loopDict[loopName],
                         parameters.tectoBase,
                         parameters.sequencingAdapters, junction,
                         f,f2)
                 
  
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
wd = parameters.wd


# check if working directory exists and if not, creates it
if not os.path.exists(wd):
    os.mkdir(wd)
    
# initialize file
filename = os.path.join(wd, 'allKTjunctionsNB.06.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

filename2 = os.path.join(wd, 'allKTjunctionsNB.06_characterization.txt')
print 'saving to %s'%filename
f2 = open(filename2, 'w')

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 1
df = pd.DataFrame(columns = ['Sequence', 'Junction Topology','Loop', 'Receptor', 'HelixContext', 'Junction Sequence',
            'HelixSequence', 'HelixLength1', 'Helixlength2','Junctionlength', 'Totallength']) 
#add header
df.to_csv(f2, sep='\t')


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
junction.motif = [('KT')]
for KTname in KTnames:
    countAll = count + 1
    KTseq = KTMotifs[KTname]
    junction.sequences[0] = KTseq
    print(junction.sequences)    
    helices = Helix(parameters.KThelixDict[helixName],0 ).alongHelix()
    f, count = saveSet(junction, helices, helixName, receptorName, loopName, f,f2, count)
    #what is this count variable?
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

junction = Junction(junctionMotif[0])
    
KTseq = KTMotifs[KTmotif_base]  
junction.sequences[0] = KTseq
junction.motif =[('KT')]
       
# now save
helices = Helix(parameters.KThelixDict[helixName], 0).centralRegion()
f, count = saveSet(junction, helices, helixName, receptorName, loopName, f,f2, count)





# close
f.close()
f2.close()
logfile.close()
