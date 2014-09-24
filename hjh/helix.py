#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

##### IMPORT #####
import numpy as np
import globalvars
parameters = globalvars.Parameters()

class Helix(object):
    def __init__(self, helixSequence, junctionLength):
        
        self.sequence = helixSequence
        if self.testHelix():
            print 'helix is good'
        self.effectiveLength = len(helixSequence[0])-junctionLength
        
    
    def testHelix(self):
        """
        Make sure the Helix is properly formatted
        """
        sequence = self.sequence
        isHelixGood = True
        
        if not isinstance(sequence, tuple):
            print 'Check input: helix is not in tuple format.\n\tex: (\'AAAA\', \'TTTT\')'
            isHelixGood = False
        if len(sequence[0]) != len(sequence[1]):
            print 'Check input: helix sides are not of equal length'
            isHelixGood = False
        return isHelixGood
        

    def formatHelices(self, totalLengths, lengthsOneHelix):
        """
        given the sequence of the whole helix, the total
        length of the desired helix (i.e. can be 7-12nt),
        and the desired length on the receptor-side of the junction,
        (i.e. helix One), give a new sequence of the form
        (helix1-side1, helix1-side2, helix2-side2, helix1-side2)
        """
        sequence = self.sequence
        side1 = sequence[0]
        side2 = sequence[1]
        
        if len(totalLengths) != len(lengthsOneHelix):
            print 'Error: total and side one should be same length'
        
        # initializ new helix giving all the possible lengths of
        # helix around junction
        numCombos = len(totalLengths)
        helixAll = np.array(np.empty(numCombos),
                            dtype={'names':('h1_side1', 'h2_side1', 'h2_side2', 'h1_side2'),
                                'formats':('S1000', 'S1000', 'S1000', 'S1000')})
        
        for i in range(numCombos):
            # being careful to get everything in the proper orientation
            # in linear space the order is: helix1 side1, helix2 side1, helix2 side2, helix1 side1

            helixAll['h1_side1'][i] = side1[:lengthsOneHelix[i]]
            helixAll['h1_side2'][i] = side2[::-1][:lengthsOneHelix[i]][::-1]
            
            helixAll['h2_side1'][i] = side1[::-1][:(totalLengths[i]-lengthsOneHelix[i])][::-1]
            helixAll['h2_side2'][i] = side2[:totalLengths[i]-lengthsOneHelix[i]]

        self.printHelixOneAndTwo(helixAll)
        return helixAll
        
        
    def printHelixOneAndTwo(self, helixAll):
        """
        print in the 5'-3' direction of the two helices
        """
        numCombos = len(helixAll)
        for i in range(numCombos):
            #print '%s %s\t%s %s'%(helixOne['side1'][i], helixTwo['side1'][i],
            #                      helixTwo['side2'][i], helixOne['side2'][i])
            
            print '%s %s\t%s %s'%(helixAll['h1_side1'][i],
                                  helixAll['h2_side1'][i],
                                  helixAll['h2_side2'][i],
                                  helixAll['h1_side2'][i])
            
        return
        
    def defaultLocation(self):
        """
        Return helix 1 and 2 when junciton is placed in the middle.
        If helix length is odd, returns two helices with different offsets
        """
        sequence = self.sequence
        helixLength = self.effectiveLength
    
        helixOneLengths = []
        helixTwoLengths = []
        if helixLength%2 == 0: # is even
            helixOneLengths = np.array([helixLength/2])
        else: # total length is odd
            helixOneLengths = np.array([np.floor(helixLength/2.0), np.ceil(helixLength/2.0)]).astype(int)
        
        totalLengths = np.array([helixLength]*len(helixOneLengths))

        return self.formatHelices(totalLengths, helixOneLengths)
    
    def centerLocation(self):
        """
        Return helix 1 and 2 when junction is placed in the middle.
        Returns only one helix, with helix one longer by one bp if necessary
        """
        sequence = self.sequence
        helixLength = self.effectiveLength

        return self.formatHelices([helixLength], [int(np.ceil(helixLength*0.5))])
    
    def centralRegion(self):
        """
        Return all possible helix 1 and 2 given the junction length
        by perturbing +/-2 all combos
        """
        sequence = self.sequence
        helixLength = self.effectiveLength
        
        possibleLengths = np.arange(helixLength-1, helixLength+3)
        helixOneLengths = []
        helixTwoLengths = []
        for i, possLength in enumerate(possibleLengths):
            if possLength%2 == 0: # is even
                oneSide = np.array([possLength/2-1, possLength/2, possLength/2+1])
   
            else: # total length is odd
                oneSide = np.array([np.floor(possLength/2.0), np.ceil(possLength/2.0)])
            helixOneLengths = np.append(helixOneLengths, oneSide).astype(int)
            helixTwoLengths = np.append(helixTwoLengths, oneSide[::-1]).astype(int)
            
        totalLengths = helixOneLengths + helixTwoLengths

        return self.formatHelices(totalLengths, helixOneLengths)
    
    def alongHelix(self):
        """
        Return all possible helix 1 and 2 given the junction length
        by moving junction positions to every point in helix
        """
        sequence = self.sequence
        helixLength = self.effectiveLength
        
        # for standard helix length
        helixOneLength = range(helixLength+1)
        totalLengths = [helixLength]*len(helixOneLength)
        
        # for increase helix length by one
        for i in range(1, helixLength+1):
            helixOneLength.append(i)
            totalLengths.append(helixLength+1)
        
        # for increase helix length by two
        for i in range(1, helixLength+2):
            helixOneLength.append(i)
            totalLengths.append(helixLength+2)
        
        return self.formatHelices(totalLengths, helixOneLength)
    
    def doubleDouble(self):
        """
        Return helix split into 3 parts now
        """
        sequence = self.sequence
        helixLength = self.effectiveLength
        
        side1 = sequence[0]
        side2 = sequence[1]
        middleHelixSequence = parameters.middleHelix
        side1_middle = middleHelixSequence[0]
        side2_middle = middleHelixSequence[1]
        
        # for standard helix length
        separationLength = [3,4,5,6,7]
        totalLengths = [helixLength]*len(separationLength)
        
        # for increase helix length by two
        for i in range(3, 9):
            separationLength.append(i)
            totalLengths.append(helixLength+2)
            
        # now format helix
        numCombos = len(totalLengths)
        helixAll = np.array(np.empty(numCombos),
                            dtype={'names':('h1_side1', 'h2_side1', 'h3_side1', 'h3_side2', 'h2_side2', 'h1_side2'),
                                'formats':('S1000', 'S1000', 'S1000', 'S1000', 'S1000', 'S1000')})   
        for i in range(numCombos):
            # being careful to get everything in the proper orientation
            # in linear space the order is: helix1 side1, helix2 side1, helix2 side2, helix1 side1
            break1 = int(np.floor((totalLengths[i] - separationLength[i])/2.0))
            break2 = break1 + separationLength[i]
            middleHelixInd = np.arange((len(side1)-separationLength[i])/2, (len(side1)+separationLength[i])/2)
            helixAll['h1_side1'][i] = side1[:break1]
            helixAll['h2_side1'][i] = side1_middle[(len(side1)-separationLength[i])/2:(len(side1)+separationLength[i])/2]
            helixAll['h3_side1'][i] = side1[-(totalLengths[i]-break2):]
            
            helixAll['h3_side2'][i] = side2[::-1][-(totalLengths[i]-break2):][::-1]
            helixAll['h2_side2'][i] = side2_middle[::-1][(len(side1)-separationLength[i])/2:(len(side1)+separationLength[i])/2][::-1]
            helixAll['h1_side2'][i] = side2[::-1][:break1][::-1]
  
        return helixAll