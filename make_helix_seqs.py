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
parser.add_argument('-c', '--helix_seqs', help='filename of helix contexts (side1 and side2)')
parser.add_argument('-m', '--mode',help='mode in which to insert junctions. Can be '
                    '"along" = all possible locations in 8, 9, 10, 11 bp helices, or "standard" = one length in 8, 11, and three lengths in 9 and 10 bp',
                    default='standard')
parser.add_argument('-l', '--lengths', nargs='+', type=int)
parser.add_argument('-o', '--offsets', nargs='+', type=int)
parser.add_argument('-r', '--junction_seqs', help='side1 and side2 of the junction sequence', required=True)
parser.add_argument('-f', '--flank_to_add', help='bases to add above and below junction seq. default="U_U"', default='U_U')
parser.add_argument('-out','--out_file', help='file to save output', )

def find_length(seq):
    """find length of shortest side of junction seq"""
    return min([len(s) for s in seq.split('_')])

def convert_offset_to_h1length(effective_length, offset):
        return int(np.floor(effective_length/2.) + offset)

def find_h1_lengths_all(total_length, j_len):
    """For a given junction length and total length, find the offsets for 'along'"""
    effective_length = total_length - j_len
    helix_one_lengths = range(effective_length)
    return helix_one_lengths

def find_lengths_h1lengths_along(j_len):
    """Return all lengths corresponding to 'along' mode"""
    h1_lengths = []
    lengths = []
    for total_length in [8, 9, 10, 11]:
        helix_one_lengths = find_h1_lengths_all(total_length, j_len)
        h1_lengths += helix_one_lengths
        lengths += [total_length]*len(helix_one_lengths)
    return lengths, h1_lengths

def find_lengths_h1lengths_standard(j_len):
    """Return all lengths corresponding to 'standard' mode"""
    lengths = [8,9,9,9, 10,10,10,11]
    offsets = [0,-1, 0, 1, -1, 0, 1, 0]
    h1_lengths = [convert_offset_to_h1length(length-j_len, i) for i, length in zip(offsets, lengths)]
    return lengths, h1_lengths


if __name__ == '__main__':
    args = parser.parse_args()
    helix_seqs = processing.load_file(args.helix_context)
    junction_seqs = processing.load_file(args.junction_seqs)
    
    if args.flank_to_add:
        base_before, base_after = args.flank_to_add.split('_')
        junction_seqs.loc[:, 'no_flank'] = junction_seqs.side1 + '_' + junction_seqs.side2
        junction_seqs.loc[:, 'side1'] = base_before + junction_seqs.side1 + base_after
        junction_seqs.loc[:, 'side2'] = seqfun.rc(base_after, rna=True) + junction_seqs.side1 + seqfun.rc(base_before, rna=True)
        junction_seqs.loc[:, 'flank'] = args.flank_to_add
        

    all_seqs = []
    for (idx1, helix_row), (idx2, junction_row) in itertools.product(helix_seqs.iterrows(), junction_seqs.iterrows()):
        j_len = junction_row.str.len().min()
       
        # check user-suppied args
        if args.lengths and args.offsets:
            lengths = args.lengths
            h1_lengths = [convert_offset_to_h1length(length-j_len, i) for i, length in zip(args.offsets, lengths)]
        else:
            # if not supplied, assume one of mode
            if args.mode=='standard':
                lengths, h1_lengths = find_lengths_h1lengths_standard(j_len)
            elif args.mode=='along':
                lengths, h1_lengths = find_lengths_h1lengths_along(j_len)
            else:
                print 'ERROR mode "%s" not recognized'%args.mode
                sys.exit()
            
        for h1_length, length in zip(h1_lengths, lengths):
            helix = Helix(helix_row, j_len)
            helix_df = helix.formatHelix(helix_row, length-j_len, h1_length)
            
            seq = pd.Series({'side1':helix_df.before.side1 + junction_row.side1 + helix_df.after.side1,
                             'side2':helix_df.before.side2 + junction_row.side2 + helix_df.after.side2})
            seq_data = pd.concat([helix_row.drop(['side1', 'side2']),
                                  junction_row.drop(['side1', 'side2']),
                                  pd.Series({'length':length,
                                             'helix_one_length':h1_length,
                                             'junction_length':j_len,
                                             'junction_seq':junction_row.side1 + '_' + junction_row.side2})])
            
            
            all_seqs.append(pd.concat([seq, seq_data]))
            
            
    all_seqs = pd.concat(all_seqs, axis=1).transpose()