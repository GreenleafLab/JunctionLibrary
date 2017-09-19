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
import subprocess
import itertools
import logging
from hjh import processing


#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-s', '--seqs', help='filename of things to mutate with the '
                    'positions in which to mutate (side1, side2, positions), i.e. '
                    '~/JunctionLibrary/seq_params/receptors_expt3_original.dat',
                    required=True)
parser.add_argument('-out','--out_file', help='file to save output', required=True )

if __name__=='__main__':
    args = parser.parse_args()

    receptors = processing.load_file(args.seqs)
    script = 'python ~/JunctionLibrary/mutate_seqs.py -s {in_file} -out {out_file} -p {positions}'
    working_dir = './'
    # groupby side length
    out_filenames = []
    for name, group in receptors.groupby('positions'):
        # make name machine friendly
        filename = working_dir + os.path.basename(os.path.splitext(args.seqs)[0]
                                                  + '_' + name.replace(';', '.').replace(',', ''))
        in_filename, out_filename = filename + '.dat', filename + '_muts.dat', 
        group.drop('positions', axis=1).to_csv(in_filename, index=False, sep='\t')
        call = script.format(in_file=in_filename, out_file=out_filename, positions='"%s"'%name)
        logging.info(call)
        subprocess.call(call, shell=True)
        out_filenames.append(out_filename)
    # join
    final_muts = pd.concat([processing.load_file(filename) for filename in out_filenames])     
    final_muts.to_csv(args.out_file, sep='\t', index=False)
    
    # remove intermediate files
    for filename in in_filenames + out_filenames:
        call = "rm " + filename
        logging.info(call)
        subprocess.call(call, shell=True)        
    

