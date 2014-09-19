#!/usr/bin/env python

# Author: Sarah Denny, Namita BisariaStanford University 

# How do you debug? How do you run just a section of code?
#Why are these tuples?
# Why return f?
# Creates library of construct to test NN rules

##### IMPORT #####
import numpy as np
from hjh.helix import Helix
from hjh.junction import Junction
from Bio.Seq import Seq
import subprocess

#def threadTogether(base, receptor, helix, junction, loop):
#    """
#    Given a receptor, h1, junction, h2, and loop, this will thread them
#    together into a single string.
#    
#    Each is a tuple, first entry is the side1 sequence, second is the side 2 sequence
#    
#    receptor(side 1) + h1(side 1) + junction(side 1) +
#        h2(side 1) + loop + h2(side2) + junction(side2) + h1(side 2) + receptor(side 2)
#    """
#    
#    seq = (base + receptor[0] + helix[0] + junction[0] + helix[1] + loop +
#                                helix[2] + junction[1] + helix[3] + receptor[1] +str(Seq(base).reverse_complement()))
#    return seq
#
#def countSequences(junctionSequences, helixSequences, receptor,  loop, base):
#    # count the number of possible sequences
#    count = 0
#    for junctionSequence in junctionSequences:
#        for helixSequence in helixSequences:
#            sequence = threadTogether(base,
#                                      receptorDict[receptorName],
#                                      helixSequence,
#                                      junctionSequence,
#                                      loopDict[loopName])
#            count += 1           
#    return count  
#    
#def saveSequences(junctionSequences, helixSequences, receptor,  loop, base, f):
#    # f = open fileID. Save the possible sequences
#    count = 0
#    for junctionSequence in junctionSequences:
#        for helixSequence in helixSequences:
#            name = 'junction:%s rigid:%d_%d'%('_'.join(junctionSequence),
#                                                    len(helixSequence[0]),
#                                                    len(helixSequence[1]))
#            sequence = threadTogether(base,
#                                      receptor,
#                                      helixSequence,
#                                      junctionSequence,
#                                      loop)
#            
#            # call RNAFold (vienna package) to get dot bracket notation
#            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
#            dotbracket = structure[0]
#            energy = float(structure[1].strip('()'))
#            f.write('%s\t%s\t%4.2f\t%s\n'%(name, sequence, energy, dotbracket))
#            count += 1           
#    return f, count



def satmutagenesis(base, pos):
    
    #Where do i put new subfunctions?
    # f = open fileID. Save the possible sequences
   def makeotherbases(base,pos,others)
   seqs = [];

   for num in others
    seq1 = list(base);
    seq1[pos] = others[num]; 
   
   return seqs
   
############################
   count = 0
   sequence = [];
    for pos in positions:
            orig = base(pos_s)
        if orig == 'A':
            others = 'TCG';
            seqs = makeotherbases(base,pos,others)
        elif orig == 'T'
            others = 'ACG';
        elif orig == 'C'
            others = 'ATG'
        elif orig == 'G'
            others = 'ATC'
        sequences.append
        
        

    for junctionSequence in junctionSequences:
        for helixSequence in helixSequences:
            name = 'junction:%s rigid:%d_%d'%('_'.join(junctionSequence),
                                                    len(helixSequence[0]),
                                                    len(helixSequence[1]))
            sequence = threadTogether(base,
                                      receptor,
                                      helixSequence,
                                      junctionSequence,
                                      loop)
            
            # call RNAFold (vienna package) to get dot bracket notation
            structure = subprocess.check_output('echo '+sequence+' | RNAFold --noPS', shell=True).split('\n')[1].split(' (')
            dotbracket = structure[0]
            energy = float(structure[1].strip('()'))
            f.write('%s\t%s\t%4.2f\t%s\n'%(name, sequence, energy, dotbracket))
            count += 1           
    return f, count
    
    
if __name__ == '__main__':
    #WHY TUPLE FORMAT?
    hairpin1 = {'h1': 'GCACTAGGCTTCTAGTGC'}
    buf1 = {'b1':'TCCA'}
    buf2 = {'b2':'TAG'}
    seqbase = {'let7':'ACTCCATCA'}
    hairpin2 = {'h2':'CGATGTTGCTTAACATCG'}
    
  seqsbase =    
       
    
    helixDict      = {'rigid':('AAGATCCTGG', 'CTGGGATCTT')}
    base           =  'CTAGGA'
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
                             loopDict[loopName], base, f)
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
                
        
        
    