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
from hjh import mutations
from hjh import processing

#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-s', '--seqs', help='sequences (i.e. side1 and side2) to mutate.', required=True)
parser.add_argument('-p', '--positions', help='positions to mutate. should be a string, ; separated for two sides, comma separated for positions')
parser.add_argument('-out','--out_file', help='file to save output', required=True)



if __name__ == '__main__':
    args = parser.parse_args()
    positions = args.positions
    sequences = processing.load_file(args.seqs)
    mut_sequences_all = []
    for idx, row in sequences.iterrows():
        side1, side2 = row.loc[['side1', 'side2']]
        mutated_seqs = mutations.singlesDF(row)
        
        # subset by certain positions if given
        if positions:
            # collate positions to mutate
            pos_side1, pos_side2 = [[int(i) for i in s.split(',')] for s in positions.split(';')]
            
            # find positions that should be constant
            constant_pos_side1 = [i for i in range(len(side1)) if i not in pos_side1]
            constant_pos_side2 = [i for i in range(len(side2)) if i not in pos_side2]
            
            # go through each position and ensure its the same as the ref seq
            mutated_seqs_edit = {}
            for side in ['side1', 'side2']:
                mutated_seqs_sub = mutated_seqs.loc[side]
                if side=='side1':
                    side1_mut = pd.Series({idx2:seq for idx2, seq in mutated_seqs_sub.side1.iteritems() if all([seq[i]==side1[i] for i in constant_pos_side1])}).rename(side)
                    side2_mut = mutated_seqs_sub.side2.loc[side1_mut.index]
                else:
                    side2_mut = pd.Series({idx2:seq for idx2, seq in mutated_seqs_sub.side2.iteritems() if all([seq[i]==side2[i] for i in constant_pos_side2])}).rename(side)
                    side1_mut = mutated_seqs_sub.side1.loc[side2_mut.index]
                mutated_seqs_edit[side] = pd.concat([side1_mut, side2_mut], axis=1)
            mutated_seqs_edit = pd.concat(mutated_seqs_edit)
            
            mutated_seqs = mutated_seqs_edit
            
        # save data about the mutations/positions
        mut_info = {}
        for side, side_seq in zip(['side1', 'side2'], [side1, side2]):
            mutated_seqs_sub = mutated_seqs.loc[side]
            seq_data = pd.Series({idx:'%s_%s%d%s'%(side, s2, i, s1) for idx, seq in mutated_seqs_sub.loc[:, side].iteritems() for i, (s1, s2) in enumerate(zip(seq, side_seq)) if s1!=s2})
            mut_info[side] = seq_data
        mut_info = pd.concat(mut_info)
        # group any double mutants
        mut_info = mut_info.groupby(level=[0,1]).apply(lambda x: ';'.join(x))
        
        # save info and seqs and other info
        mutated_seqs_info = pd.concat([mutated_seqs, mut_info.rename('mut_info')], axis=1)
        for col in row.drop(['side1', 'side2']).index.tolist():
            # below is a clunky name to prevent column name duplicates
            old_col_name = col
            i = 1
            while old_col_name in mutated_seqs_info:
                old_col_name = '%s_%d'%(col, i)
                i+=1
            mutated_seqs_info.loc[:, old_col_name] = row.loc[col]

        mut_sequences_all.append(mutated_seqs_info.reset_index(drop=True))
        
    # concat
    mut_sequences_all = pd.concat(mut_sequences_all, ignore_index=True)
    
    # reduce duplicates
    mut_sequences_all = mut_sequences_all.groupby(['side1', 'side2']).first().reset_index()
    mut_sequences_all.to_csv(args.out_file, index=False, sep='\t')
    sys.exit()
    
