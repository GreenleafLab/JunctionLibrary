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

import globalvars
parameters = globalvars.Parameters()
from hjh.helix import Helix
from hjh.junction import Junction
import create_library_dataframe

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
filename = os.path.join(wd, 'all3x3junctions_NB.01_helixcontext.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

filename2 = os.path.join(wd, 'all3x3junctions_NB.01_helixcontext_characterization.txt')
print 'saving to %s'%filename
f2 = open(filename2, 'w')

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 1

"""
DIFFERENT HELIX CONTEXT: One position, subset of junctions, ten other helix contexts
Save subset of junctions in one position in ten different helix contextx. 
"""
junctionMotifs = parameters.differentHelixJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixNames   = parameters.otherHelixNames
cutOffNumber = 12

df = pd.DataFrame(columns = ['Sequence', 'Junction Topology','Loop', 'Receptor', 'HelixContext', 'Junction Sequence',
            'HelixSequence', 'HelixLength1', 'Helixlength2','Junctionlength', 'Totallength']) 
#add header
df.to_csv(f2, sep='\t')

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
        count = create_library_dataframe.saveSet(junction, helices, helixName, receptorName, loopName, f,f2, logfile, count)

# close
f.close()
f2.close()
logfile.close()
