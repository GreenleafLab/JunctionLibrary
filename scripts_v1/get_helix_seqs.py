import os
import numpy as np
import pandas as pd
import itertools
from fittinglibs import fileio
from tectolibs import tectplots
from hjh.junction import Junction
from hjh.helix import Helix


# to start, load old helix model
model = fileio.loadFile('/home/sarah/JunctionLibrary/seq_params/linear_regression_length10.p')
seqs_per_length = {}
for length in [9, 10, 11]:
    # get all helix seqs of a particular length
    all_seqs = Junction(tuple(['W']*(length-1) + ['G'])).sequences
    
    # get rid of homoplymeric tracts
    max_homo_length = 3
    bases = ['A', 'C', 'G', 'U']
    max_num_GC_or_AU = np.ceil(length/2.)+1
    sub_index = []
    for idx, seq in all_seqs.side1.iteritems():
        too_homopolymeric = any([seq.find(base*(max_homo_length+1))>-1 for base in bases])
        too_gc_rich = len([s for s in seq if s=='C' or s=='G']) > max_num_GC_or_AU
        too_au_rich = len([s for s in seq if s=='A' or s=='U']) > max_num_GC_or_AU
        
        if not too_homopolymeric and not too_gc_rich and not too_au_rich:
            sub_index.append(idx)
    seqs_sub = all_seqs.loc[sub_index].copy()
    
    # predict how they would behave by linear model
    helix_seqs_sep = pd.concat(list(seqs_sub.side1.str), axis=1, ignore_index=True)
    if length == 10:
        helix_seqs_sep = helix_seqs_sep.loc[:, range(length-1)]
    elif length == 9:
        # repeat position 4
        helix_seqs_sep = helix_seqs_sep.loc[:, range(4) + [4] + range(4, length-1)].copy()
        helix_seqs_sep.columns = range(9)
    elif length == 11:
        # delete position 5
        helix_seqs_sep = helix_seqs_sep.loc[:, range(4) + range(5, length-1)].copy()
        helix_seqs_sep.columns = range(9)

    X = pd.get_dummies(helix_seqs_sep)
    dG_predicted = pd.Series(model.predict(X), index=seqs_sub.index)
    
    # choose 100 uniformly
    interval = int(np.ceil(len(dG_predicted)/5000.))
    # choose from whole interval
    dG_sub1 = tectplots.choose_uniform(dG_predicted.iloc[::interval], N=100)
    
    # choose from 95% of interval
    dG_sub2_pre = dG_predicted.iloc[1::interval]
    min_dG, max_dG = dG_predicted.quantile([0.025, 0.975])
    dG_sub2 = tectplots.choose_uniform(dG_sub2_pre.loc[(dG_sub2_pre>min_dG)&(dG_sub2_pre<max_dG)], N=100)
    dG_sub =pd.concat([dG_sub1, dG_sub2]).rename('dG_predicted')
    
    seqs_per_length[length] = pd.concat([seqs_sub, dG_sub], axis=1).dropna(subset=['dG_predicted'])

# save

seqs_per_length = pd.concat(seqs_per_length, names=['length', 'index']).reset_index(level='length')
seqs_per_length.loc[:, ['side1', 'side2', 'length']].to_csv('~/JunctionLibrary/seq_params/helix_sides_random_200.dat', sep='\t', index=False)