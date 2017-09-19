#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

# This function calls dependent functions to make the subset of
# final library that is all three by three junctions and less


##### IMPORT #####
import numpy as np
import pandas as pd
import os
import sys
import argparse
import itertools
from operator import mul
import subprocess
import cPickle as pickle
from matplotlib import gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import functools
import multiprocessing
import ipdb

# load custom libraries
from fittinglibs import seqfun
from hjh import processing
import hjh.tecto_assemble
from hjh.helix import Helix
from hjh.junction import Junction
from hjh.tecto_assemble import TectoSeq


#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-p', '--predefined_helix_seqs', help='filename of helix contexts with the '
                    'junction location already defined (h1_side1, h2_side1, h2_side2, h2_side1).'
                    ' If given, all other options (i.e. length, offsets) are ignored.')

parser.add_argument('-r', '--junction_seqs', help='side1 and side2 of the junction sequence', required=True)
parser.add_argument('-out','--out_file', help='file to save output', )

def find_length(seq):
    """find length of shortest side of junction seq"""
    return min([len(s) for s in seq.split('_')])




if __name__ == '__main__':
    args = parser.parse_args()
    helix_seqs = processing.load_file(args.predefined_helix_seqs)
    junction_seqs = processing.load_file(args.junction_seqs)
        
    all_seqs = []
    for (idx1, helix_row), (idx2, junction_row) in itertools.product(helix_seqs.iterrows(), junction_seqs.iterrows()):

        seq = pd.Series({'side1':helix_row.h1_side1 + junction_row.side1 + helix_row.h2_side1,
                         'side2':helix_row.h2_side2 + junction_row.side2 + helix_row.h1_side2})
        seq_data = pd.concat([helix_row.drop(['h1_side1', 'h1_side2', 'h2_side1', 'h2_side2']),
                              junction_row.drop(['side1', 'side2']),
                              pd.Series({'helix_seq':helix_row.h1_side1 + '_' + helix_row.h2_side1 + '&' + helix_row.h2_side2 + '_' + helix_row.h1_side2,
                                         'junction_seq':junction_row.side1 + '_' + junction_row.side2})])
        
        all_seqs.append(pd.concat([seq, seq_data]))
            
            
    all_seqs = pd.concat(all_seqs, axis=1).transpose()
    all_seqs.to_csv(args.out_file, index=False, sep='\t' )