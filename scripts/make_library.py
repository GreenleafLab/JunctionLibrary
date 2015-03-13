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

# load custom libraries
from create_library import threadTogether
from hjh.helix import Helix
from hjh.junction import Junction


#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-map','--map_file', help='file that contains parameters to vary', required=True)
parser.add_argument('-sub','--subset_number', help='maximum number of sequences per junction motif. Default is all possible ')
parser.add_argument('-seq','--junction_sequences', help='overwrite the junction motifs given in map with particular sequenecs in this file')
parser.add_argument('-par','--seq_params', help='location of seq_param file. Default is ~/JunctionLibrary/seq_params/seq_params.txt')
parser.add_argument('-out','--out_file', help='file to save output')

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()
if args.out_file is None: args.out_file = os.path.splitext(args.map_file)[0] + '.library'
if args.seq_params is None: args.seq_params = os.path.expanduser('~/JunctionLibrary/seq_params/seq_params.txt')

# load the parameters
seq_params = pd.read_table(args.seq_params, index_col=[0,1])

# read map file
expt_params = pd.read_table(args.map_file)

# if seq file is defined, replace junction motifs with this
if args.junction_sequences is not None:
    expt_params.loc[:, 'junction'] = ['user_defined'] +[np.nan]*(len(expt_params)-1)

# iterate through all
cols = expt_params.columns.tolist()+['junction_seq', 'helix_seq','tecto_sequence', 'sequence']
numToLog =  np.product((~expt_params.isnull()).sum().values)
allSeqs = pd.DataFrame(columns=cols)
logSeqs = pd.DataFrame(index=np.arange(numToLog), columns=expt_params.columns.tolist()+['number'])

for i, x in enumerate(itertools.product(*[expt_params.loc[:,name].dropna() for name in expt_params])):
    
    x = pd.Series(x, index=expt_params.columns.tolist())
    if x.junction == 'user_defined':
        junction = Junction(sequences=pd.read_table(args.junction_sequences))
    else:                
        junction = Junction(tuple(x.junction.split(',')))
        helix = Helix(seq_params.loc[('helix', x.helix)], junction.length   , x.offset, int(x.length))
    
    allSeqSub = pd.DataFrame(index=junction.sequences.index, columns=cols)
    logSeqs.loc[i] = x
    logSeqs.loc[i, 'number'] = len(junction.sequences)

    
    for loc in junction.sequences.index:
        if x.junction == 'user_defined':    # user defined junctions can have any length- so it depends on sequence, not motif
            helix = Helix(seq_params.loc[('helix', x.helix)],
                          junction.findEffectiveJunctionLength(sequence=junction.sequences.loc[loc]),
                          x.offset, int(x.length))
        x.loc['junction_seq'] = '_'.join(junction.sequences.loc[loc])
        x.loc['helix_seq'] = '&'.join(['_'.join(helix.split.loc['side1']), '_'.join(helix.split.loc['side2'])])
        sequence = seq_params.loc[('loop', x.loop), 'loop']
        sequence = threadTogether(sequence, [helix.split.loc['side1', 'after'], helix.split.loc['side2', 'before']])
        sequence = threadTogether(sequence, junction.sequences.loc[loc])
        sequence = threadTogether(sequence, [helix.split.loc['side1', 'before'], helix.split.loc['side2', 'after']])
        sequence = threadTogether(sequence, seq_params.loc[('receptor', x.receptor), ['side1', 'side2']])
        sequence = threadTogether(sequence, seq_params.loc[('base', x.loc['base']), ['side1', 'side2']])
        x.loc['tecto_sequence'] = sequence
        x.loc['sequence']  = threadTogether(sequence, seq_params.loc[('adapters', x.adapters), ['side1', 'side2']])
        allSeqSub.loc[loc] = x
    
    allSeqs = allSeqs.append(allSeqSub, ignore_index=True)