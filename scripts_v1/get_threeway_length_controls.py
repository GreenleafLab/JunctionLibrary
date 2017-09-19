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
                                     'side2':row.j2[:-1] + row.j3[1:],
                                     'j_len_diff':min(len(row.j1), len(row.j2+row.j3)-2)-2})
        elif loop_context == 'L2':
            two_way_seq = pd.Series({'side1':row.j1[:-1] + row.j2[1:],
                                     'side2':row.j3,
                                     'j_len_diff':min(len(row.j1)+len(row.j2)-2, len(row.j3))-2})
        junction_seq_controls[(idx, loop_context)] = pd.concat([
            two_way_seq, row.drop(['j1', 'j2', 'j3'])])

        
junction_seq_controls = pd.concat(junction_seq_controls, names=['index', 'loop_context']).unstack().swaplevel(0,1).sort_index()
junction_seq_controls_sub = pd.concat({name:group for name, group in junction_seq_controls.groupby('j_len_diff') if name <= 3}, names=['j_len_diff'])
# find those with j_len_offset of 1, 2, or 3
# these original and control junctions will be saved separately

for j_len_diff, loop_context in itertools.product([1,2,3], ['L1', 'L2']):
    junction_seq_controls_sub.loc[(j_len_diff, loop_context)].to_csv(
        '~/JunctionLibrary/seq_params/three_way_junction_controls_%s_%d.dat'%(loop_context, j_len_diff), sep='\t', index=False)
    
    junction_seqs.loc[junction_seq_controls_sub.loc[(j_len_diff, loop_context)].index].to_csv(
        '~/JunctionLibrary/seq_params/three_way_junctions_%d.dat'%j_len_diff,  sep='\t', index=False)
    
    
    
    
    
    
    