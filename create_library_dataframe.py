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
import globalvars
import pandas as pd
parameters = globalvars.Parameters()

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

def doubleDoubleJunctionOrder(helixSequence, junctionSequence, receptor, base, sequencingAdapters):
    """
    In the case where you have two junctions per helix and three subhelices,
    i.e. double double configuration,
    Threading order is more complicated
    """
    sequenceList = [(helixSequence['h3_side1'], helixSequence['h3_side2']),
                    (junctionSequence['j2_side1'], junctionSequence['j2_side2']),
                    (helixSequence['h2_side1'], helixSequence['h2_side2']),
                    (junctionSequence['j1_side1'], junctionSequence['j1_side2']),
                    (helixSequence['h1_side1'], helixSequence['h1_side2']),
                    (receptor[0], receptor[1]),
                    (base[0], base[1]),
                    (sequencingAdapters[0], sequencingAdapters[1])]
    return sequenceList

def twoWayJunctionOrder(helixSequence, junctionSequence, receptor, base, sequencingAdapters):
    """
    for a single two-way junction (one entry each side) and single helix breakpoint
    (two entries each side), this is the way to thread together, adjoining each entry to either
    side of central point (i.e. loop)
    """
    sequenceList = [(helixSequence['h2_side1'], helixSequence['h2_side2']),
                    (junctionSequence['side1'], junctionSequence['side2']),
                    (helixSequence['h1_side1'], helixSequence['h1_side2']),
                    (receptor[0], receptor[1]),
                    (base[0], base[1]),
                    (sequencingAdapters[0], sequencingAdapters[1])]
    return sequenceList

def saveSequenceDoubleDouble(names, junctionSequences, helixSequences, receptor,  loop, base, sequencingAdapters, f,f2):
    # f = open fileID. Save the possible sequences
    df = pd.DataFrame(columns = ['Sequence', 'Junction Topology','Loop', 'Receptor', 'HelixContext', 'Junction Sequence',
            'HelixSequence', 'HelixLength1', 'Helixlength2','Junctionlength', 'Totallength']) 
    count = 0      
    data = [];      

    # for every junction sequence included by 'junctionMotif'
    for junctionSequence in junctionSequences:
        
        # for every helix orientation/sequence given by 'helices'
        for helixSequence in helixSequences:
            
            # give a name to each junction/helix orientation
            name = '%s.%s.%d_%d'%('.'.join(names),
                                  '_'.join(junctionSequence),
                                            len(helixSequence['h2_side1']),
                                            len(''.join(helixSequence))/2)
            
            # find the list of tuples that are added to either side of a central region to
            # obtain the final sequence.
            sequenceList = doubleDoubleJunctionOrder(helixSequence, junctionSequence, receptor, base, sequencingAdapters)
            
            # starting with the loop, progressively add each tuple in sequence list to either side
            # only add up to tecto Base before assessing structure
            sequence = loop
            for outside in sequenceList[:-1]:
                sequence = threadTogether(sequence, outside)
            sequence_i = sequence
                
            # call RNAFold (vienna package) to get dot bracket notation
            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
            dotbracket = structure[0]
            energy = float(structure[1].strip('()'))
            
            # now add sequencing adapters
            sequence = threadTogether(sequence, sequenceList[-1])
            
            # save to open file
            f.write('%s\t%s\t%d\t%4.2f\t%s\n'%(name, sequence, len(sequence), energy, dotbracket))
            
            
            #this puts an "N' for the 
            newtup = tuple('N' if x == '' else x for x in junctionSequence)

            junctionsequence_parse = '_'.join(newtup)
            helixsequence_parse = '_'.join(helixSequence)
             
            
            #should have done this with an array, can fix?
            lengths = {'len_h1s1': len(helixSequence['h1_side1']),'len_h2s1': len(helixSequence['h2_side1']), 'len_h2s2': len(helixSequence['h2_side2']),
            'len_h1s2': len(helixSequence['h1_side2']), 'len_js1': len(junctionSequence[0]) + len(junctionSequence[1]), 'len_js2' : len(junctionSequence[2]) + len(junctionSequence[3])}
            total_length = len(helixSequence['h1_side1']) + len(helixSequence['h1_side2']) + 2
            junctionlength = 2
            #junction motif is D_D for double-double
            data = (sequence_i, 'D_D', names[2], names[1], names[0], junctionsequence_parse, helixsequence_parse, lengths['len_h1s1'],
            lengths['len_h2s1'], junctionlength,total_length )
            
            df.loc[count] = data
            count += 1
    df.to_csv(f2, header = False,sep='\t')        
    return count
  
def saveSequences(names, junctionSequences, helixSequences, receptor,  loop, base, sequencingAdapters, junction, f,f2):
    # f = open fileID. Save the possible sequences
    df = pd.DataFrame(columns = ['Sequence', 'Junction Topology','Loop', 'Receptor', 'HelixContext', 'Junction Sequence',
            'HelixSequence', 'HelixLength1', 'Helixlength2','Junctionlength', 'Totallength']) 
    count = 0      
    data = [];                           
    # for every junction sequence included by 'junctionMotif'
    
    for junctionSequence in junctionSequences:
        
        # for every helix orientation/sequence given by 'helices'
        for helixSequence in helixSequences:
           # import pdb; pdb.set_trace()
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
            sequence_i = sequence       
            # call RNAFold (vienna package) to get dot bracket notation
            
            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
            dotbracket = structure[0]
            energy = float(structure[1].strip('()'))
            
            # now add sequencing adapters
            
            sequence = threadTogether(sequence, sequenceList[-1])
            
            # save to open file
            f.write('%s\t%s\t%d\t%4.2f\t%s\n'%(name, sequence, len(sequence), energy, dotbracket))
            
            junctionsequence_parse = '_'.join(junctionSequence)
            helixsequence_parse = '_'.join(helixSequence)
            
            
            #should have done this with an array, can fix?
            lengths = {'len_h1s1': len(helixSequence['h1_side1']),'len_h2s1': len(helixSequence['h2_side1']), 'len_h2s2': len(helixSequence['h2_side2']),
            'len_h1s2': len(helixSequence['h1_side2']), 'len_js1': len(junctionSequence[0]), 'len_js2' : len(junctionSequence[1])}
            total_length = len(helixSequence['h1_side1']) + len(helixSequence['h1_side2']) + junction.length
            
            
            data = (sequence_i, '_'.join(junction.motif), names[2], names[1], names[0], junctionsequence_parse, helixsequence_parse, lengths['len_h1s1'],
            lengths['len_h2s1'], junction.length,total_length )
            
            df.loc[count] = data
            
            
            count += 1
            #df = pd.DataFrame( columns = ['Sequence', 'Junction Topology','Loop', 'Receptor', 'Base', 'HelixContext', 'Junction Sequence',
            # 'HelixSequence', 'HelixLength1', 'Helixlength2',Junctionlength',' 'Totallength'])    
            
    df.to_csv(f2, header = False,sep='\t')
    return count

def saveSet(junction, helices, helixName, receptorName, loopName, f, f2, logfile, countAll, savetype=None):
    """
    this is just a space saver. Calls the saveSequences command for a given
    set of junctions, helices,  receptor, and loop)
    """
    
    if savetype is None:
        count = saveSequences([helixName, receptorName, loopName],
                             junction.sequences,
                             helices,
                             parameters.receptorDict[receptorName],
                             parameters.loopDict[loopName],
                             parameters.tectoBase,
                             parameters.sequencingAdapters, junction,
                             f,f2)
    elif savetype == 'double':
        count = saveSequenceDoubleDouble([helixName, receptorName, loopName],
                             junction.sequences,
                             helices,
                             parameters.receptorDict[receptorName],
                             parameters.loopDict[loopName],
                             parameters.tectoBase,
                             parameters.sequencingAdapters, junction,
                             f,f2)
        
    print 'Saved %d sequences...'%(count+countAll)
    logfile.write('%s\t%s\t%s\t%s\t%d\t%d\t%d\t\t%d\n'%(receptorName, loopName, helixName,
                                                      ','.join(junction.motif),
                                          junction.howManyPossibilities(),
                                          len(junction.sequences),
                                          len(helices),
                                          count
                                          ))
    return count

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
   
    