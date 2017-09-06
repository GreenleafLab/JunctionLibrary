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
parser.add_argument('-a', '--starting_seq', help='seed sequence, i.e. a loop. Must hve "seq" column.',
                    default='seq_params/loop_1.dat')
parser.add_argument('-b', '--add_seqs', nargs="+", help='list of filenames of sequences to add. '
                    'All should have columns "side1" and "side2" of the adapter sequences')
parser.add_argument('-out','--out_file', help='file to save output', required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    
    # load seqs
    loop_seqs = processing.load_file(args.starting_seq)
    other_seqs = [processing.load_file(filename) for filename in args.add_seqs]
    
    # starting with the loop, thread together the sequential pieces of the tectoRNA
    new_seqs = loop_seqs
    new_seqs.loc[:, 'starting_seq'] = new_seqs.seq
    for add_seqs in other_seqs:
        all_seqs = []
        for (idx1, row1), (idx2, row2) in itertools.product(new_seqs.iterrows(), add_seqs.iterrows()):
            seq = processing.thread_together(row1.seq, (row2.side1, row2.side2))
            seq_data = pd.concat([pd.Series({'seq':seq}),
                row1.drop('seq'), row2.drop(['side1', 'side2'])])
            all_seqs.append(seq_data)
        all_seqs = pd.concat(all_seqs, axis=1).transpose()
        new_seqs = all_seqs

    # save
    new_seqs.to_csv(args.out_file, sep='\t', index=False)
    
    sys.exit()
    
