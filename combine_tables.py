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
parser.add_argument('-a', '--add_seqs', nargs="+", help='list of filenames of sequence data to add.')
parser.add_argument('-i', '--libname', help='[optional] the name to prepend to sublibraries, '
                    'i.e. "tertcontact". default is "library"', default='library')
parser.add_argument('-u', '--unique', action="store_true", help='[optional] whether to ensure '
                    '"seq" column is unique. If set, will take first instance of each seq. default=False.')
parser.add_argument('-out','--out_file', help='file to save output', required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    # load
    new_seqs = pd.concat({'%s_%d'%(args.libname, i):processing.load_file(filename)
                          for i, filename in enumerate(args.add_seqs)}, names=['sublibrary', 'index'])
    new_seqs.reset_index(level=0, inplace=True)
    
    # make unique if option given
    if args.unique:
        old_cols = new_seqs.columns.tolist()
        new_seqs = new_seqs.groupby('seq').first().reset_index().loc[:, old_cols]
        
    
    # save
    new_seqs.to_csv(args.out_file, sep='\t', index=False)
    

    
