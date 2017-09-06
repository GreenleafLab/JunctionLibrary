#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library


##### IMPORT #####
import numpy as np
import pandas as pd
import os
import sys
import argparse
import itertools

from hjh import processing

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


if __name__ == '__main__':
    args = parser.parse_args()
    
    # load seqs
    adapter_seqs = processing.load_file(args.adapter_seqs)
    base_seqs = processing.load_file(args.base_seqs)
    receptor_seqs = processing.load_file(args.receptor_seqs)
    helix_seqs = processing.load_file(args.helix_seqs)
    loop_seqs = processing.load_file(args.loops)
    
    # starting with the loop, thread together the sequential pieces of the tectoRNA
    new_seqs = loop_seqs
    for add_seqs in [helix_seqs, receptor_seqs, base_seqs]:
        all_seqs = []
        for (idx1, row1), (idx2, row2) in itertools.product(new_seqs.iterrows(), add_seqs.iterrows()):
            seq = processing.thread_together(row1.seq, (row2.side1, row2.side2))
            seq_data = pd.concat([pd.Series({'seq':seq}),
                row1.drop('seq'), row2.drop(['side1', 'side2'])])
            all_seqs.append(seq_data)
        all_seqs = pd.concat(all_seqs, axis=1).transpose()
        new_seqs = all_seqs
    # save seq without adapters
    new_seqs.loc[:, 'tecto_sequence'] = new_seqs.seq
    
    # add adapter_seqs
    all_seqs = []
    for (idx1, row1), (idx2, row2) in itertools.product(new_seqs.iterrows(), adapter_seqs.iterrows()):
        seq = processing.thread_together(row1.seq, (row2.side1, row2.side2))
        seq_data = pd.concat([pd.Series({'seq':seq}),
            row1.drop('seq'), row2.drop(['side1', 'side2'])])
        all_seqs.append(seq_data)
    new_seqs = pd.concat(all_seqs, axis=1).transpose()
    
    # save
    new_seqs.to_csv(args.out_file, sep='\t', index=False)
    
    sys.exit()
    
