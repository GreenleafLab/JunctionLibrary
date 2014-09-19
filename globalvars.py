#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

# This contains all of the hard-coded variables called during the library generation

##### IMPORT #####
import numpy as np

class Parameters():
    """
    stores helix sequences, loop sequences, base sequences, and library construct sequences
    """
    def __init__(self):
        # tecto RNA parameters
        """
        helices include rigid, watson crick, and then 10 other essentially random helices.
        helices were generated to have particular GC content (0.25, 0.5, 0.75), and then
        run through Fang's code ('simu_tecto.py') to find energetics of tecto construct.
        Final energy was sorted from lowest to highest, and these are a chosen subset of
        the original 30 helices. See HJH_project/data/simulations/9_15_14_tectoRNA_sims.xlxs
        """
        self.helixDict = {'rigid':  ('AAGATCCTGG', 'CTGGGATCTT'),
                          'wc':     ('AAGATCCTCG', 'CGAGGATCTT'),
                            'h02': ('CATATGTACT', 'AGTACATATG'),
                            'h06': ('CCTGATTTGT', 'ACAAATCAGG'),
                            'h08': ('TAAATCAGCC', 'GGCTGATTTA'),
                            'h10': ('CAGCACTGCC', 'GGCAGTGCTG'),
                            'h12': ('GATTGCCCTT', 'AAGGGCAATC'),
                            'h14': ('AAAGCCTTGA', 'TCAAGGCTTT'),
                            'h16': ('AGCCGGCAAG', 'CTTGCCGGCT'),
                            'h21': ('AGGCCGGCTG', 'CAGCCGGCCT'),
                            'h25': ('GCGGCTGTAA', 'TTACAGCCGC'),
                            'h28': ('AAGAAACGAC', 'GTCGTTTCTT'),
        }
        self.standardHelixNames = np.array(['rigid', 'wc'])
        
        self.allHelixNames = np.sort(self.helixDict.keys())
        self.otherHelixNames = self.allHelixNames[np.logical_not(np.in1d(self.allHelixNames, self.standardHelixNames))]

        self.loopDict = {'goodLoop': 'GGAA',
                         'badLoop': 'GAAA',
                         'KL1': 'AAGATTCCA',
                         'KL2': 'AAGATTCA',
                         'KL3': 'AAGATTA',
                         'KL4': 'ACACAA',
                         'KL5': 'AAGATA',
                         'KL6': 'ACACAGGA'}
        """
        loops: KL1-4 are designed to hybridize to R2-receptor that binds GGAA tetraloop.
        With A's added to either side.
        """
        
        
        self.receptorDict   = ({'R1':('TATGG', 'CCTAAG'),
            'R1_2':('', 'CCTAAG'),
                                'KL1': ('AGGTCGGA', 'A'),
                                'KL2': ('AU', 'AGGTCGGA'),
                                'KL3': ('AU', 'AGTCGGA'),
                                'KL4': ('AU', 'AGGTCGA'),
                                'KL5': ('AU', 'AGTCGA'),
                                'KL6': ('A', 'AGGTCGGA'),
                                'KL7': ('A', 'AGTCGGA'),
                                'KL8': ('A', 'AGGTCGA'),
                                'KL9    ': ('A', 'AGTCGA'),
                                    })
        """
        receptors KL1 through KL9 are different binding parters of flow piece 3.
        R1_2 is different binding parnter to flow piece 2. Needs different tectoBase.
        
        R1_2_base = ('CTAG', 'CTAG')
        R1_2 = ('ATGG', 'CCTAAGTC')
        
        """

        self.tectoBase = ('CTAGGA', 'TCCTAG')
        
        self.allJunctions = [('',),('M',), ('M', 'M'), ('B1',), ('B2',), ('W',),
                            ('M','M','M'),
                            ('B1', 'B1'), ('B1', 'B1', 'B1'), ('B2', 'B2'), ('B2', 'B2', 'B2'),
                            ('M',  'B1', 'B1'), ('M',  'B1'), ('M',  'M', 'B1'),
                            ('B2', 'M',  'M' ), ('B2', 'M' ), ('B2', 'B2', 'M')]
        
        self.alongJunctions = [('M',), ('M', 'M'), ('B1',), ('B2',), ('W',)]
        self.centralJunctions =    [('M','M','M'),
                                    ('B1', 'B1'), ('B1', 'B1', 'B1'), ('B2', 'B2'), ('B2', 'B2', 'B2'),
                                    ('M',  'B1', 'B1'), ('M',  'B1'), ('M',  'M', 'B1'),
                                    ('B2', 'M',  'M' ), ('B2', 'M' ), ('B2', 'B2', 'M')]
        
        # sequencing library parameters
        self.sequencingAdapters = ('TTGTATGGAAGACGTTCCTGGATCC', 'AGATCGGAAGAGCGGTTCAGC')