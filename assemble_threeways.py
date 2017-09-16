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
parser.add_argument('-c', '--helix_seqs', help='filename of helix contexts (side1 and side2)')
parser.add_argument('-r', '--junction_seqs', help='three_way junction sequence (with ', required=True)
parser.add_argument('-out','--out_file', help='file to save output', required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    
    # load seqs
    helix_seqs = processing.load_file(args.helix_seqs)
    loop_seqs = ['GGAA', 'UUCG']
    junction_seqs = processing.load_file(args.junction_seqs)
    
    # for each of the junction sequences, put in one context
    new_seqs = []
    for (idx1, helix_row), (idx2, junction_row) in itertools.product(helix_seqs.iterrows(), junction_seqs.iterrows()):
        # assemble seq
        for loop1, loop2 in [loop_seqs, loop_seqs[::-1]]:
            junction_seq = processing.convert_junction_to_seq(junction_row)
            base_helix, stem_loop_seqs = processing.convert_threeway_helix_to_stemloops(helix_row, [loop1, loop2])
            # find the sequence of the junction and stem loops
            seq = processing.convert_threeway_to_single(junction_seq, stem_loop_seqs)
            # add the base helix
            seq_with_base = base_helix.side1 + seq + base_helix.side2
            # save annotations
            annotations = pd.concat([helix_row.drop(processing.return_threeway_helix_fields()),
                                     junction_row.drop(processing.return_junction_fields()),
                                     pd.Series({'loop_seq':loop1+'_'+loop2})])
            new_seqs.append(pd.concat([pd.Series({'seq':seq_with_base}), annotations]))
    new_seqs = pd.concat(new_seqs, axis=1).transpose()    
    
    # save
    
    
    new_seqs.to_csv(args.out_file, sep='\t', index=False)
    
    sys.exit()
    
