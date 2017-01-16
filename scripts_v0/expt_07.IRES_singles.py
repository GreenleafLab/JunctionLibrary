from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np
import sys
import hjh.junction_seqs
# functions


## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/ires_singles/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = hjh.junction_seqs.getSingleMutationIres()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files
    lengths = [8, 9, 10, 11, 12]
    sides = ['up', 'down']

    num_lengths = len(lengths)
    expt_map_constant = pd.Series(index = ['junction','receptor', 'helix','loop', 'base',   'adapters', 'offset'],
                                  data  = ['defunct', '11nt',  'wc',  'GGAA', 'normal', 'truseq', 0], dtype=str)
    cols = expt_map_constant.index.tolist()


    expt_map = pd.DataFrame(index = np.arange(num_lengths),
                        columns = cols + ['length'])
    expt_map.loc[0, cols] = expt_map_constant
    expt_map.loc[np.arange(len(sides)), 'side'] = sides
    expt_map.loc[np.arange(len(lengths)), 'length'] = lengths
    filename = os.path.join(saveDir, 'expt.length.map')
    expt_map.to_csv(filename, index=False, sep='\t')

    print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s "%(filename,
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     )
    sys.exit()  # and run make_library commands
    
    filenames = [os.path.join(saveDir, 'all.junctions_to_compare.library')]
    allSeqs = []
    # load allSeqs
    for filename in filenames:
        allSeqSub = pd.read_table(filename+'.txt', index_col=0)
        allSeqSub.loc[:, 'tecto_object'] = np.nan
        with open(filename+'.pkl', 'rb') as input:
            for loc in allSeqSub.index:
                allSeqSub.loc[loc, 'tecto_object'] = pickle.load(input)
        allSeqs.append(allSeqSub)
    allSeqs = pd.concat(allSeqs, ignore_index=True)
    
    # print varna sequences fro moving A bulge and think about it
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'UUUG_UGUAUG')&(allSeqs.loc[:, 'flank'] == 'CC')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand()[0]
        os.system(cmnd.replace('test.png', 'junction_%d.png'%loc))