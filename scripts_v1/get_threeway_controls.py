##### IMPORT #####
import numpy as np
import pandas as pd
import os
import sys
import argparse
import itertools

from hjh import processing


# load sequences
junction_seqs = processing.load_file('~/JunctionLibrary/seq_params/three_way_junctions.dat')

# cut out the 'extra'
junction_seq_controls = {}
for idx, row in junction_seqs.iterrows():
    for loop_context in ['L1','L2']:
        if loop_context == 'L1':
            two_way_seq = pd.Series({'side1':row.j1,
                                     'side2':row.j2[0] + row.j3[-1]})

        elif loop_context == 'L2':
            two_way_seq = pd.Series({'side1':row.j1[0] + row.j2[-1],
                                     'side2':row.j3})
        junction_seq_controls[(idx, loop_context)] = pd.concat([
                two_way_seq, row.drop(['j1', 'j2', 'j3'])])

        
junction_seq_controls = pd.concat(junction_seq_controls, names=['index', 'loop_context']).unstack().swaplevel(0,1).sort_index()
for loop_context in ['L1', 'L2']:
    junction_seq_controls.loc[loop_context].to_csv('~/JunctionLibrary/seq_params/three_way_junction_controls_%s.dat'%loop_context, sep='\t', index=False)