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
parser.add_argument('-s', '--seqs',  help='filenames with column "seq" to test ss structure of.')
parser.add_argument('-out','--out_file', help='file to save output', required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    seqs = processing.load_file(args.seqs)
    seqs.loc[:, 'ss'] = processing.check_ss_structure_set(seqs)
    
    ## save
    seqs.to_csv(args.out_file, index=False, sep='\t')