#!/usr/bin/env python

# Author: Sarah Denny, Namita Bisaria Stanford University 

# Takes helix DNA library and outputs text file with sequence, followed by 
# features of each variant



# This function calls dependent functions to make the subset of
# final library that is all three by three junctions and less

##### IMPORT #####
import numpy as np
import pandas as pd
import os
import sys



# load custom libraries
import create_library_dataframe
import globalvars
parameters = globalvars.Parameters()
from hjh.helix import Helix
from hjh.junction import Junction



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
    
# initialize file and data file
filename = os.path.join(wd, 'all3x3junctions_NB.00_standard.txt')
filename_data = os.path.join(wd, 'all3x3junctions_NB_characterization.txt')
print 'saving to %s'%filename
f = open(filename, 'w')

f2 = open(filename_data, 'w') 

# initalize log file
logfile = open(os.path.join(wd, '%s.log'%os.path.splitext(filename)[0]), 'w')

# initialize counts
count = 0
   

"""
STANDARD: One position, Rigid and Watson Crick
Save all junctions in one position in two different helix contextx. 
"""
junctionMotifs = parameters.differentHelixJunctions
receptorName = 'R1'
loopName     = 'goodLoop'
helixNames   = parameters.standardHelixNames

df = pd.DataFrame(columns = ['Sequence', 'Junction Topology','Loop', 'Receptor', 'HelixContext', 'Junction Sequence',
            'HelixSequence', 'HelixLength1', 'Helixlength2','Junctionlength', 'Totallength']) 
#add header
df.to_csv(f2, sep='\t')

for junctionMotif in junctionMotifs:
    junction = Junction(junctionMotif)
    for helixName in helixNames:
        # helices in default location
        helices = Helix(parameters.helixDict[helixName], junction.length).centerLocation()
           

        count = create_library_dataframe.saveSet(junction, helices, helixName, receptorName, loopName, f, f2, logfile, count)
        
# close
f2.close()
f.close()
logfile.close()
