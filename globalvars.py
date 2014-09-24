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
        """
        where to save files
        """
        self.wd = '~/Dropbox/HJH_project/HJH_project/libraries/v1'
    
        
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
        
        # for if there are two junctions
        self.middleHelix = ('CAGCACTGCC', 'GGCAGTGCTG')
        
        # list all helix names
        self.allHelixNames = np.sort(self.helixDict.keys())
        
        # the helices that aren't 'wc' and 'rigid'
        self.otherHelixNames = self.allHelixNames[np.logical_not(np.in1d(self.allHelixNames, self.standardHelixNames))]
        
        """
        LOOPS:
        - goodLoop is the GGAA tetraloop, and badLoop is a nonbinding control
        - R2 * loops are designed to hybridize to R2-receptor that binds GGAA tetraloop.
        - can be complementary to either side of the R2 receptor (side1 and side2)
        - A's appended to either side.
        """
        
        self.loopDict = {   'goodLoop': 'GGAA',
                            'badLoop': 'GAAA',
                            'R2_side1_0_6': 'GGAA',
                            'R2_side1_0_5': 'ACCAGAA',
                            'R2_side1_1_6': 'ACAGATA',
                            'R2_side1_2_7': 'AAGATTA',
                            'R2_side1_0_4': 'ACCAGA',
                            'R2_side1_1_5': 'ACAGAA',
                            'R2_side1_2_6': 'AAGATA',
                            'R2_side1_3_6': 'AGATA',
                            'R2_side2_1_7': 'ACACAGGA',
                            'R2_side2_1_6': 'ACACAGA',
                            'R2_side2_2_7': 'AACAGGA',
                            'R2_side2_0_5': 'AACACAA',
                            'R2_side2_1_5': 'ACACAA',
                            'R2_side2_2_6': 'AACAGA',
                            'R2_side2_3_7': 'ACAGGA',
                            'R2_side2_1_4': 'ACACA',
                        }

        """
        RECEPTORS
        -receptors HIV *  are different binding parters of flow piece 3 (HIV kissing loop).
        -receptors FF * are different binding partners of the 'For Free' flow piece
            that should bind to the 11nt receptor. These pieces should make it even easier to
            bind there, and have different lengths, etc.
        
        """
        
        self.receptorDict   = { 'R1': ('TATGG', 'CCTAAG'),
                                'HIV_side1_A_0_6': ('AGGTCGGAGG', 'CCA'),
                                'HIV_side1_AA_0_6': ('AGGTCGGAGG', 'CCAA'),
                                'HIV_side2_A_0_6': ('AGG', 'CCAGGTCGGA'),
                                'HIV_side2_A_0_5': ('AGG', 'CCAGGTCGA'),
                                'HIV_side2_A_1_6': ('AGG', 'CCAGTCGGA'),
                                'HIV_side2_A_0_4': ('AGG', 'CCAGGTCA'),
                                'HIV_side2_A_1_5': ('AGG', 'CCAGTCGA'),
                                'HIV_side2_A_2_6': ('AGG', 'CCATCGGA'),
                                'HIV_side2_AA_0_6': ('AAGG', 'CCAGGTCGGA'),
                                'HIV_side2_AA_0_5': ('AAGG', 'CCAGGTCGA'),
                                'HIV_side2_AA_1_6': ('AAGG', 'CCAGTCGGA'),
                                'HIV_side2_AA_0_4': ('AAGG', 'CCAGGTCA'),
                                'HIV_side2_AA_1_5': ('AAGG', 'CCAGTCGA'),
                                'HIV_side2_AA_2_6': ('AAGG', 'CCATCGGA'),
                                'HIV_side2_AA_0_3': ('AAGG', 'CCAGGTA'),
                                'HIV_side2_AA_3_6': ('AAGG', 'CCACGGA'),
                                'FF_side1_A_0_6': ('ATAAGTCAGG', 'CCA'),
                                'FF_side1_AA_0_6': ('ATAAGTCAGG', 'CCAA'),
                                'FF_side2_A_0_6': ('AGG', 'CCATAAGTCA'),
                                'FF_side2_A_0_5': ('AGG', 'CCATAAGTA'),
                                'FF_side2_A_1_6': ('AGG', 'CCAAAGTCA'),
                                'FF_side2_A_0_4': ('AGG', 'CCATAAGA'),
                                'FF_side2_A_1_5': ('AGG', 'CCAAAGTA'),
                                'FF_side2_A_2_6': ('AGG', 'CCAAGTCA'),
                                'FF_side2_AA_0_6': ('AAGG', 'CCATAAGTCA'),
                                'FF_side2_AA_0_5': ('AAGG', 'CCATAAGTA'),
                                'FF_side2_AA_1_6': ('AAGG', 'CCAAAGTCA'),
                                'FF_side2_AA_0_4': ('AAGG', 'CCATAAGA'),
                                'FF_side2_AA_1_5': ('AAGG', 'CCAAAGTA'),
                                'FF_side2_AA_2_6': ('AAGG', 'CCAAGTCA'),
                                'FF_side2_AA_0_3': ('AAGG', 'CCATAAA'),
                                'FF_side2_AA_3_6': ('AAGG', 'CCAGTCA'),
                                'BR_AA_AA': ('AAGG', 'CCAA')}
        

        """
        The base of the tectoRNA forms a hairpin below the receptor.
        """
        self.tectoBase = ('CTAGGA', 'TCCTAG')
        
        
        
        self.allJunctions = [('',),('M',), ('M', 'M'), ('B1',), ('B2',), ('W',),
                            ('M','M','M'),
                            ('B1', 'B1'), ('B1', 'B1', 'B1'), ('B2', 'B2'), ('B2', 'B2', 'B2'),
                            ('M',  'B1', 'B1'), ('M',  'B1'), ('M',  'M', 'B1'),
                            ('B2', 'M',  'M' ), ('B2', 'M' ), ('B2', 'B2', 'M')]
        self.differentHelixJunctions = [('',),('M',), ('M', 'M'), ('B1',), ('B2',), 
                            ('M','M','M'),
                            ('B1', 'B1'), ('B1', 'B1', 'B1'), ('B2', 'B2'), ('B2', 'B2', 'B2'),
                            ('M',  'B1', 'B1'), ('M',  'B1'), ('M',  'M', 'B1'),
                            ('B2', 'M',  'M' ), ('B2', 'M' ), ('B2', 'B2', 'M')]
        
        self.alongJunctions = [('M',), ('M', 'M'), ('B1',), ('B2',), ('W',),
            ('M', 'M', 'B1'), ('B2', 'M', 'M')]
        
        self.centralJunctions =    [('M','M','M'),
                                    ('B1', 'B1'), ('B1', 'B1', 'B1'), ('B2', 'B2'), ('B2', 'B2', 'B2'),
                                    ('M',  'B1', 'B1'), ('M',  'B1'), ('M',  'M', 'B1'),
                                    ('B2', 'M',  'M' ), ('B2', 'M' ), ('B2', 'B2', 'M')]
        
        # these are the junctions that are both central and different loop mutants
        self.differentLoopJunctions = [('',),('M',), ('M', 'M'), ('W',),
                            ('M','M','M'),
                            ('B1', 'B1'), ('B1', 'B1', 'B1'), ('B2', 'B2'), ('B2', 'B2', 'B2'),
                            ('M',  'B1', 'B1'), ('M',  'M', 'B1'),
                            ('B2', 'M',  'M' ), ('B2', 'B2', 'M')]
        self.sameLoopJunctions = [('B1', ), ('B2',), ('M', 'B1'), ('M', 'B2')]
        
        # sequencing library parameters
        self.sequencingAdapters = ('TTGTATGGAAGACGTTCCTGGATCC', 'AGATCGGAAGAGCGGTTCAGC')