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

from hjh import processing, mutations


# load seq_params
seq_params = pd.read_table('/home/sarah/JunctionLibrary/seq_params/seq_params.txt', index_col=(0,1)).loc['receptor', ['side1', 'side2']]

# for receptors that are the same length as the 11nt, find all evolutionary trajectories
receptor_wt = seq_params.loc['11nt']

# starting recpetors, what evolutionary tracks lead to each.
receptor_starts = ['11nt', '11nt', '11nt',  'C7.2', 'IC3']
receptor_ends =   ['C7.2', 'IC3',  'R_Vc2', 'IC3', 'R_Vc2']

new_receptors = {}
for i, (receptor, receptor2) in enumerate(zip(receptor_starts, receptor_ends)):
    steps, graph, paths, paths_parsed = mutations.getPathsFromFormatted(seq_params.loc[receptor].copy(), seq_params.loc[receptor2].copy())
    new_receptors[(receptor, receptor2)] = pd.concat([paths_parsed.loc[:, col].str.replace('_', '') for col in paths_parsed], axis=1)
new_receptors = pd.concat(new_receptors, names=['start', 'end'])

seqs = new_receptors.reset_index().groupby(['side1','side2']).first().reset_index()
