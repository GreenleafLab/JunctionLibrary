#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

##### IMPORT #####
import numpy as np
import itertools
from Bio.Seq import Seq

##### Define CLASSES #####


class Junction(object):
    def __init__(self, junctionMotif):
        
        self.motif = junctionMotif
        if self.testJunction(): # is junction properly formatted?
            self.sequences = self.mapJunctionToSequence()
            self.length = self.findEffectiveJunctionLength()
        
    def testJunction(self):
        """
        Is Junction properly formatted? i.e. are all submotifs allowable?
        """
        isJunctionGood = True
        motif = self.motif
        # first test if it is a tuple:
        if not isinstance(motif, tuple):
            print 'Check input: junction is not in tuple format.\n\tex: (\'M\',) or (\'B1\', \'B1\')\n'
            isJunctionGood = False
        for submotif in motif:
            if submotif == 'M' or submotif == 'B1' or submotif == 'B2' or submotif == 'W':
                pass
            else:
                print('%s\n%s\n%s\n%s')%('Improperly formatted junction subMotif. Proper notation is:',
                           '\t\'M\' for Mismatch',
                           '\t\'B1\' for Bulge side 1, \'B2\' for Bulge Side 2',
                           '\t\'W\' for Watson Crick bp\n')
                isJunctionGood = False
                break
        return isJunctionGood
        
    def mapSubmotifToNucleotide(self, submotif):
        """
        A 'submotif' is a mismatch ('M'), a Watson-Crick base pair ('W'),
        or a bulge ('B'). Return possible bases on each side of junction.
        """
        if submotif == 'M':
            # if mismatch, list all possible bp, then get rid of the ones that are watson crick
            numberOfPossibilities = 12
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                      dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A']*3 + ['T']*3 + ['G']*3 + ['C']*3
            possibleBases['side2'] = ['A', 'G', 'C', 'T', 'G', 'C', 'A', 'T', 'G', 'A', 'T', 'C']

        elif submotif == 'W':
            # if watson crick, just list Watson crick bp
            numberOfPossibilities = 4
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                      dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A', 'T', 'G', 'C']
            possibleBases['side2'] = ['T', 'A', 'C', 'G']
        
        elif submotif == 'B1':
            # if bulge on side 1, side 1 is every base, opposite a blank
            numberOfPossibilities = 4
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A', 'T', 'G', 'C']
            possibleBases['side2'] = ['']*4
            
        elif submotif == 'B2':
            # if bulge on side 2, side 2 is every base, opposite a blank 
            numberOfPossibilities = 4
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['']*4
            possibleBases['side2'] = ['A', 'T', 'G', 'C']

        return possibleBases
    
    def howManyPossibilities(self):
        """
        Given a junction motif in the form ('M','M','M'),
        return how many possible combos there are
        """
        motif = self.motif
        numberOfCombos = 1
        for submotif in motif:
            if submotif == 'M':
                numberOfCombos *= 12
            elif submotif == 'B1' or submotif == 'B2' or submotif =='W':
                numberOfCombos *= 4
                
        return numberOfCombos
    
    def mapJunctionToSequence(self):
        """
        Given a junction motif in the form ('M','M','M'),
        return all possible side 1 and 2. Both are read 5'-3'
        """
        motif = self.motif
        side1 = ['']
        side2 = ['']
        for submotif in motif:
            possibleBases = self.mapSubmotifToNucleotide(submotif)
            side1 = [''.join(result) for result in itertools.product(possibleBases['side1'], side1)]
            side2 = [''.join(result[::-1]) for result in itertools.product(possibleBases['side2'], side2)]
        
        numberOfCombos = self.howManyPossibilities()       
        junctionSequences = np.array(np.empty(numberOfCombos),
                                     dtype={'names':('side1', 'side2'), 'formats':('S1000','S1000')})
        junctionSequences['side1'] = side1
        junctionSequences['side2'] = side2
        
        # print
        print 'Number of possible combinations: %d'%(numberOfCombos)
        print 'First 16 combos: '
        for i in range(min(16, numberOfCombos)):
            print '%s\t%s'%(junctionSequences['side1'][i], junctionSequences['side2'][i])
        
        return junctionSequences
    
    def findEffectiveJunctionLength(self):
        """
        Given a junction motif, how long does it count for?
        Mismatch is 1, BP is 1, bulge is 0
        """
        motif = self.motif
        effectiveLength = 0
        for submotif in motif:
            if submotif == 'M' or submotif == 'W':
                effectiveLength += 1
        
        return effectiveLength
        
        

