#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

##### IMPORT #####
import numpy as np
import itertools
import pandas as pd

##### Define CLASSES #####


class Junction(object):
    def __init__(self, junctionMotif=None, sequences=None):
        
        if junctionMotif is not None:
            self.motif = junctionMotif
            if self.testJunction(): # is junction properly formatted?
                self.sequences = self.mapJunctionToSequence()
                self.length = self.findEffectiveJunctionLength(motif=junctionMotif)
        if sequences is not None:
            self.sequences = sequences
        
        
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
            if submotif == 'M' or submotif == 'B1' or submotif == 'B2' or submotif == 'W' or submotif == '' or submotif == 'CG' or submotif == 'GC' or submotif == 'AU' or submotif == 'UA' or submotif == 'N'   :
                pass
            else:
                print('%s\n%s\n%s\n%\n%s\n%s\n')%('Improperly formatted junction subMotif. Proper notation is:',
                           '\t\'M\' for Mismatch',
                           '\t\'B1\' for Bulge side 1, \'B2\' for Bulge Side 2',
                           '\t\'W\' for Watson Crick bp',
                           "\t'N' for any mismatch OR WC base pair",
                           "\t'CG', 'GC', 'AU', 'UA' for any particular WC base pair")
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
            possibleBases['side1'] = ['A']*3 + ['U']*3 + ['G']*3 + ['C']*3
            possibleBases['side2'] = ['A', 'G', 'C', 'U', 'G', 'C', 'A', 'U', 'G', 'A', 'U', 'C']

        elif submotif == 'W':
            # if watson crick, just list Watson crick bp
            numberOfPossibilities = 4
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                      dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A', 'U', 'G', 'C']
            possibleBases['side2'] = ['U', 'A', 'C', 'G']
        
        elif submotif == 'B1':
            # if bulge on side 1, side 1 is every base, opposite a blank
            numberOfPossibilities = 4
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A', 'U', 'G', 'C']
            possibleBases['side2'] = ['']*numberOfPossibilities
            
        elif submotif == 'B2':
            # if bulge on side 2, side 2 is every base, opposite a blank 
            numberOfPossibilities = 4
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['']*numberOfPossibilities
            possibleBases['side2'] = ['A', 'U', 'G', 'C']
        
        elif submotif == 'GC':
            # if bulge on side 2, side 2 is every base, opposite a blank 
            numberOfPossibilities = 1
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['G']
            possibleBases['side2'] = ['C']
             
        elif submotif == 'CG':
            # if bulge on side 2, side 2 is every base, opposite a blank 
            numberOfPossibilities = 1
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['C']
            possibleBases['side2'] = ['G']
            
        elif submotif == 'AU':
            # if bulge on side 2, side 2 is every base, opposite a blank 
            numberOfPossibilities = 1
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A']
            possibleBases['side2'] = ['U']

        elif submotif == 'UA':
            # if bulge on side 2, side 2 is every base, opposite a blank 
            numberOfPossibilities = 1
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['U']
            possibleBases['side2'] = ['A']
        
        elif submotif == '':
            # if no junction, i.e. for straight 'rigid'
            numberOfPossibilities = 1
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['']*numberOfPossibilities
            possibleBases['side2'] = ['']*numberOfPossibilities

        elif submotif == 'N':
            # if no junction, i.e. for straight 'rigid'
            numberOfPossibilities = 16
            possibleBases = np.array(np.empty(numberOfPossibilities),
                                    dtype={'names':('side1', 'side2'), 'formats':('S1','S1')})
            possibleBases['side1'] = ['A']*4 + ['U']*4 + ['G']*4 + ['C']*4
            possibleBases['side2'] = ['A', 'G', 'C', 'U',  'A', 'G', 'C', 'U', 'A', 'G', 'C', 'U','A', 'G', 'C', 'U',]
            
        return pd.DataFrame(possibleBases)
    
    def howManyPossibilities(self):
        """
        Given a junction motif in the form ('M','M','M'),
        return how many possible combos there are
        """
        motif = self.motif
        numberOfCombos = 1
        for submotif in motif:
            if submotif == 'N':
                numberOfCombos *= 16
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
            """
            given possibile nucleotides for each junction side, what are all possible permutations of nucleotides?
            """
            possibleBases = self.mapSubmotifToNucleotide(submotif)
            side1 = pd.Series([''.join(result) for result in itertools.product(possibleBases['side1'], side1)])
            side2 = pd.Series([''.join(result[::-1]) for result in itertools.product(possibleBases['side2'], side2)])
        
        numberOfCombos = self.howManyPossibilities()
        junctionSequences = pd.concat([side1, side2], axis=1, keys=['side1', 'side2']).astype(str)
        #junctionSequences = np.array(np.empty(numberOfCombos),
        #                             dtype={'names':('side1', 'side2'), 'formats':('S1000','S1000')})
        #junctionSequences['side1'] = side1
        #junctionSequences['side2'] = side2
        
        # print
        print 'Number of possible combinations: %d'%(numberOfCombos)
        print 'First 16 combos: '
        for i in range(min(16, numberOfCombos)):
            print '%s\t%s'%(junctionSequences['side1'][i], junctionSequences['side2'][i])
        
        return junctionSequences
    
    def findEffectiveJunctionLength(self, motif=None, sequence=None):
        """
        Given a junction motif, how long does it count for?
        Mismatch is 1, BP is 1, bulge is 0
        """
        effectiveLength = 0
        if motif is not None:
            for submotif in motif:
                if submotif == 'M' or submotif == 'W':
                    effectiveLength += 1
        if sequence is not None:
            effectiveLength = np.min([len(sequence.loc[side]) for side in ['side1', 'side2']])
            
        
        return effectiveLength
        
        

