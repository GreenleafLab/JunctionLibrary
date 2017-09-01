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
import hjh.tecto_assemble
from hjh.helix import Helix
from hjh.junction import Junction
from hjh.tecto_assemble import TectoSeq


#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-l', '--loops', help='loop sequences',
                    default='seq_params/loop_1.dat')
parser.add_argument('-a', '--adapter_seqs', help='side1 and side2 of the adapter sequences',
                    default='seq_params/adapters_1.dat')
parser.add_argument('-b', '--base_seqs', help='side1 and side2 of the base sequence helix of tectoRNA',
                    default='seq_params/base_1.dat')
parser.add_argument('-r', '--receptor_seqs', help='side1 and side2 of the receptor sequence', required=True)
parser.add_argument('-s', '--helix_seqs', help='side1 and side2 of the helix sequences', required=True)
parser.add_argument('-out','--out_file', help='file to save output', required=True)

def thread_together(inside, outside):
    return (outside[0] + inside + outside[1])

def load_file(filename):
    ext = os.path.splitext(filename)[1]
    if ext == '.csv':
        return pd.read_csv(filename)
    elif ext == '.dat':
        return pd.read_table(filename)


if __name__ == '__main__':
    args = parser.parse_args()
    
    adapter_seqs = load_file(args.adapter_seqs)
    base_seqs = load_file(args.base_seqs)
    receptor_seqs = load_file(args.receptor_seqs)
    helix_seqs = load_file(args.helix_seqs)
    loop_seqs = load_file(args.loops)
    
    new_seqs = loop_seqs
    
    for add_seqs in [helix_seqs, receptor_seqs, base_seqs, adapter_seqs]:
        all_seqs = []
        for (idx1, row1), (idx2, row2) in itertools.product(new_seqs.iterrows(), add_seqs.iterrows()):
            seq = thread_together(row1.seq, (row2.side1, row2.side2))
            seq_data = pd.concat([pd.Series({'seq':seq}),
                row1.drop('seq'), row2.drop(['side1', 'side2'])])
            all_seqs.append(seq_data)
        all_seqs = pd.concat(all_seqs, axis=1).transpose()
        new_seqs = all_seqs
    
    new_seqs.to_csv(args.out_file, sep='\t', index=False)
    
    sys.exit()
    
