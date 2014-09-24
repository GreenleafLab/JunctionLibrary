#!/usr/bin/env python

# Author: Sarah Denny, Stanford University

# How to find new helix sequences?
# Want helices that cover a range of free energies according to Fang's simu_tecto code

from hjh.seqfun import reverseComplement

def getHelixGC(percentGC, helixLength):
    bits = np.random.binomial(1, percentGC, helixLength)
    helixSequence = ['']*helixLength
    for i, base in enumerate(bits):
        coinflip = np.random.binomial(1, 0.5)
        if coinflip == 0:
            if base == 0 :
                helixSequence[i] = 'A'
            else:
                helixSequence[i] = 'G'
        else:
            if base == 0 :
                helixSequence[i] = 'T'
            else:
                helixSequence[i] = 'C'
        
    return  ''.join(helixSequence)

def findGCcontent(helixSequence):
    count = 0
    for base in helixSequence:
        if base == 'G' or base == 'C':
            count +=1
    return count/float(len(helixSequence))

gc_contents = np.linspace(0.25, 0.75, 3)
numTrials = 10
helixLength = 10
for gc_content in gc_contents:
    for trial in range(numTrials):
        helixSequence = getHelixGC(gc_content, helixLength)
        print '%4.2f\t%s\t%s\t%4.2f'%(gc_content, helixSequence, reverseComplement(helixSequence), findGCcontent(helixSequence))
        

