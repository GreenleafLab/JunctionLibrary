##### IMPORT #####
import numpy as np
import pandas as pd
import os
import sys
import argparse
import itertools

from hjh import processing


# load sequences
for filename in ['~/JunctionLibrary/seq_params/three_way_helices.dat',
                 '~/JunctionLibrary/seq_params/three_way_helices2.dat',
                 '~/JunctionLibrary/seq_params/three_way_helices_minus1.dat',
                 '~/JunctionLibrary/seq_params/three_way_helices_minus2.dat']:

    helix_seqs = processing.load_file(filename)
    
    # split into base helix, L1 helix, and L2 helix
    keys = ['base_side1', 'base_side2']
    base_helix = helix_seqs.loc[:, keys].rename(columns={key:'h1_'+key.split('_')[-1] for key in keys}).copy()
    
    for loop_context in ['L1', 'L2']:
        if loop_context == 'L1':
            keys = ['h1_side1', 'h1_side2']
        elif loop_context == 'L2':
            keys = ['h2_side1', 'h2_side2']
        else:
            keys=None
        loop_helix = helix_seqs.loc[:, keys].rename(columns={key:'h2_'+key.split('_')[-1] for key in keys}).copy()
    
        predefined_helix =pd.concat([base_helix, loop_helix], axis=1)
        outfilename = os.path.splitext(filename)[0] + '_split_%s.dat'%loop_context
        predefined_helix.to_csv(outfilename, index=False, sep='\t')    
    
        # do the same after subtracting 1, 2 and 3 base pairs from base helix.
    
        
