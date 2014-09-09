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
    
def saveSequences(junctionSequences, helixSequences, receptor,  loop, base, f):
    # f = open fileID. Save the possible sequences
    count = 0

    # for every junction sequence included by 'junctionMotif'
    for junctionSequence in junction.sequences:
        
        # for every helix orientation/sequence given by 'helices'
        for helixSequence in helices:
            
            # give a name to each junction/helix orientation
            name = 'junction:%s rigid:%d_%d'%('_'.join(junctionSequence),
                                            len(helixSequence['h1_side1']),
                                            len(helixSequence['h2_side1']))
            
            # find the list of tuples that are added to either side of a central region to
            # obtain the final sequence.
            sequenceList = [(helixSequence['h2_side1'], helixSequence['h2_side2']),
                            (junctionSequence['side1'], junctionSequence['side2']),
                            (helixSequence['h1_side1'], helixSequence['h1_side2']),
                            (receptor[0], receptor[1]),
                            (base[0], base[1])]
            
            # starting with the loop, progressively add each tuple in sequence list to either side
            sequence = loop
            for outside in sequenceList:
                sequence = threadTogether(sequence, outside)
                    
            # call RNAFold (vienna package) to get dot bracket notation
            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
            dotbracket = structure[0]
            energy = float(structure[1].strip('()'))
            
            # save to open file
            f.write('%s\t%s\t%4.2f\t%s\n'%(name, sequence, energy, dotbracket))
            count += 1
            
    return f, count

if __name__ == '__main__':
    
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
        
    # 3) Alternate Helix sequence
    # AU rich vs GC rich
    # same AU-GC, 5 others
    # change everything but GU + tandem mismatches there.
    # 
    
    # next step: change receptors
    # program in specific sequences of junctions (i.e. kink turns)
    # allow changing helix content
    
#def saveSequences(junctionMotifs, helixContext, helixDict, receptorDict,  loopDict, base):
#    
#    for loopName in loopDict.keys():
#        for receptorName in receptorDict.keys():
#            for helixName in helixDict.keys():
#                # initialize file saving
#                fileName = [helixName,
#                            helixContext,
#                            '_'.join([''.join(result) for result in junctionMotifs]),
#                            loopName,
#                            receptorName]
#                f = open('.'.join(fileName)+'.txt', 'w')
#                
#                # also Count the number of motifs
#                count = 0 
#                
#                # cycle through junction Motifs
#                for motif in junctionMotifs:
#                    
#                    # for each of the junction motifs, find all possible junction sequences
#                    junction = Junction(motif)
#                    
#                    # given the topology of the junction (i.e. length), what is the
#                    # surrounding sequence?
#                    if helixContext == 'central':
#                        helices = Helix(helixDict[helixName], junction.length).centralRegion()
#                    elif helixContext == 'default':
#                        helices = Helix(helixDict[helixName], junction.length).defaultLocation()
#                    elif helixContext == 'along':
#                        helices = Helix(helixDict[helixName], junction.length).alongHelix()
#                    elif helixContext == 'doubledouble':
#                        helices = Helix(helixDict[helixName], junction.length).doubleDouble()
#                    else:
#                        print "Error: check preference for helix context.\nEither 'central', 'default', 'along', or 'doubledouble'"
#            
#                    # for each junction sequence and helix sequence, print line
#                    for junctionSequence in junction.sequences:
#                        for helix in helices:
#                            name = 'junction:%s rigid:%d_%d'%('_'.join(junctionSequence),
#                                                            len(helix[0]),
#                                                            len(helix[1]))
#                            sequence = threadTogether(base,
#                                                      receptorDict[receptorName],
#                                                      helix,
#                                                      junctionSequence,
#                                                      loopDict[loopName])
#                            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
#                            dotbracket = structure[0]
#                            energy = float(structure[1].strip('()'))
#                            f.write('%s\t%s\t%4.2f\t%s\n'%(name, sequence, energy, dotbracket))
#                            count += 1
#                f.close()            
#    return count   
                
        
        
    