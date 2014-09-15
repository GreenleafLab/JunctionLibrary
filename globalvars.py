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
                          'h3':     ('', ''),
                          'h4':     ('', ''),
                          'h5':     ('', ''),
                          'h6':     ('', ''),
                          'h7':     ('', ''),
                          'h8':     ('', ''),
                          'h9':     ('', ''),
                          'h10':    ('', ''),
        }
        
        
        self.loopDict = {'GGAA': 'GGAA',
                         'GAAA': 'GAAA'}
        self.receptorDict   = ({'R1':('TATGG', 'CCTAAG'),
                                'KL1': ('', ''),
                                'KL2': ('', ''),
                                'KL3': ('', ''),
                                'KL4': ('', ''),
                                'KL5': ('', ''),
                                'KL6': ('', '')
                                    })

        self.tectoBase = ('CTAGGA', 'TCCTAG')
        
        # sequencing library parameters
        self.sequencingAdapters = ('TTGTATGGAAGACGTTCCTGGATCC', 'AGATCGGAAGAGCGGTTCAGC')