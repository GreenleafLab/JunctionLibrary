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
        self.helixDict = {'rigid':  ('AAGATCCTGG', 'CTGGGATCTT'),
                          'wc':     ('AAGATCCTCG', 'CGAGGATCTT'),
                          'h01':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h02':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h03':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h04':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h05':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h06':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h07':    ('AAAAAAAAAA', 'TTTTTTTTTT'),
                          'h08':    ('AAAAAAAAAA', 'TTTTTTTTTT'),

        }
        self.allHelixNames = np.sort(self.helixDict.keys())
        
        
        self.loopDict = {'goodLoop': 'GGAA',
                         'badLoop': 'GAAA'}
        self.receptorDict   = ({'R1':('TATGG', 'CCTAAG'),
                                'KL1': ('', ''),
                                'KL2': ('', ''),
                                'KL3': ('', ''),
                                'KL4': ('', ''),
                                'KL5': ('', ''),
                                'KL6': ('', '')
                                    })

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