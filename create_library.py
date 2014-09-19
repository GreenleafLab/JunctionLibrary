#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

##### IMPORT #####
import numpy as np
from hjh.helix import Helix
from hjh.junction import Junction
from hjh.seqfun import reverseComplement
import subprocess

def threadTogether(inside, outside):
    """
    Given a sequence 'inside' and a tuple 'outside' ('side1', 'side2'),
    add  the sequences from 'outside' to either side of 'inside'.
    
    outside(side 1) + inside + outside(side 2)
    """
    
    seq = (outside[0] + inside + outside[1])
    return seq

def countSequences(junctionSequences, helixSequences, receptor,  loop, base):
    # count the number of possible sequences
    count = 0
    for junctionSequence in junctionSequences:
        for helixSequence in helixSequences:
            sequence = threadTogether(base,
                                      receptorDict[receptorName],
                                      helixSequence,
                                      junctionSequence,
                                      loopDict[loopName])
            count += 1           
    return count

def twoWayJunctionOrder(helixSequence, junctionSequence, receptor, base, sequencingAdapters):
    """
    for a single two-way junction (one entry each side) and single helix breakpoint
    (two entries each side), this is the way to thread together, adjoining each entry to either
    side of central point (i.e. loop)
    """
    print helixSequence
    print junctionSequence
    sequenceList = [(helixSequence['h2_side1'], helixSequence['h2_side2']),
                    (junctionSequence['side1'], junctionSequence['side2']),
                    (helixSequence['h1_side1'], helixSequence['h1_side2']),
                    (receptor[0], receptor[1]),
                    (base[0], base[1]),
                    (sequencingAdapters[0], sequencingAdapters[1])]
    
    return sequenceList
    
def saveSequences(names, junctionSequences, helixSequences, receptor,  loop, base, sequencingAdapters, f):
    # f = open fileID. Save the possible sequences
    count = 0

    # for every junction sequence included by 'junctionMotif'
    for junctionSequence in junctionSequences:
        
        # for every helix orientation/sequence given by 'helices'
        for helixSequence in helixSequences:
            
            # give a name to each junction/helix orientation
            name = '%s.%s.%d_%d'%('.'.join(names),
                                  '_'.join(junctionSequence),
                                            len(helixSequence['h1_side1']),
                                            len(helixSequence['h2_side1']))
            
            # find the list of tuples that are added to either side of a central region to
            # obtain the final sequence.
            sequenceList = twoWayJunctionOrder(helixSequence, junctionSequence, receptor, base, sequencingAdapters)
            
            # starting with the loop, progressively add each tuple in sequence list to either side
            # only add up to tecto Base before assessing structure
            sequence = loop
            for outside in sequenceList[:-1]:
                sequence = threadTogether(sequence, outside)
                    
            # call RNAFold (vienna package) to get dot bracket notation
            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
            dotbracket = structure[0]
            energy = float(structure[1].strip('()'))
            
            # now add sequencing adapters
            sequence = threadTogether(sequence, sequenceList[-1])
            
            # save to open file
            f.write('%s\t%s\t%d\t%4.2f\t%s\n'%(name, sequence, len(sequence), energy, dotbracket))
            count += 1
            
    return f, count

if __name__ == '__main__':
    
    """
    Example run shown below. Actual functionality should come from different file, i.e.
    'threeBythree.py' that call global parameters for the helix and loop sequences. 
    
    """
    helixDict      = {'rigid':('AAGATCCTGG', 'CTGGGATCTT')}
    base           =  ('CTAGGA', 'TCCTAG')
    loopDict       = {'GGAA':'GGAA'}
    receptorDict   = {'R1':('TATGG', 'CCTAAG')}

    
    # 2) Along Helix
    
    # define parameters
    junctionMotifs = [('M',),
        ('B1',),
        ('B2',),
        ('W',)]
    receptorName = 'R1'
    loopName     = 'GGAA'
    helixName    = 'rigid'
    
    # open file to save
    fileName = [helixName,
                'along',
                loopName,
                receptorName,
                '_'.join([''.join(result) for result in junctionMotifs])]
    f = open('.'.join(fileName)+'.txt', 'w')
    
    # start count
    countAll = []

    # loop through junctions and possible locations for junctions.
    # save all to open file id 'f'
    for junctionMotif in junctionMotifs:
        junction = Junction(junctionMotif)
        helices = Helix(helixDict[helixName], junction.length).alongHelix()
        f, count = saveSequences(junction.sequences,
                                helices,
                                receptorDict[receptorName],
                                loopDict[loopName],
                                base,
                                f)
        countAll.append(count)
    f.close()
   
    