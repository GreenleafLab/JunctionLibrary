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

from hjh import processing

receptor_filename = '~/JunctionLibrary/seq_params/receptors_expt2_original.dat'
receptors = processing.load_file(receptor_filename)
script = 'python ~/JunctionLibrary/mutate_seqs.py -s {in_file} -out {out_file} -p {positions}'
working_dir = './'
# groupby side length
out_filenames = []
for name, group in receptors.groupby('positions'):
    # make name machine friendly
    filename = working_dir + os.path.basename(os.path.splitext(receptor_filename)[0]
                                              + '_' + name.replace(';', '.').replace(',', ''))
    in_filename, out_filename = filename + '.dat', filename + '_muts.dat', 
    group.drop('positions', axis=1).to_csv(in_filename, index=False, sep='\t')
    call = script.format(in_file=in_filename, out_file=out_filename, positions='"%s"'%name)
    print call
    subprocess.call(call, shell=True)
    out_filenames.append(out_filename)
# join
final_muts = pd.concat([processing.load_file(filename) for filename in out_filenames])     
final_muts.to_csv('~/JunctionLibrary/seq_params/receptors_expt2_muts.dat', sep='\t', index=False)

