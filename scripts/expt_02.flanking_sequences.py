from hjh.junction import Junction
import itertools
import os
import pandas as pd
import sys
import hjh.mutations
from hjh.junction import Junction
import numpy as np

# functions
def getAllJunctionSeqs():
    # save junctions : first B1 junctions
    junctionMotifs = []
    junctionSeqs = {}
    
    for motif in ['B1', 'B1,B1', 'B1,B1,B1', 'N', 'B1,M']:

        junctionSeqs[motif] = {}

        junctionMotif = 'G,W,'+motif+',W,C'
        junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
        junctionSeq.loc[:, 'n_flank'] = 2

        if motif == 'B1,B1,B1' or motif == 'B1,M':
            # only take those sequences with an A in the middle
            indices = []
            for index in junctionSeq.index:
                if junctionSeq.loc[index, 'side1'][1] == 'A':
                    indices.append(index)
            junctionSeq = junctionSeq.loc[indices]
        flank = 'GC'
        junctionSeqs[motif][flank] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank'])
    
    return pd.concat(junctionSeqs,  names=['junction'])

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/flanking_sequences/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files
    lengths = [9, 10]
    offsets = [-1]
    
    num_lengths = len(lengths)
    expt_map_constant = pd.Series(index = ['helix','junction','receptor','loop', 'base',   'adapters'],
                                  data  = ['wc',   'defunct', '11nt',    'GGAA', 'normal', 'truseq'], dtype=str)
    cols = expt_map_constant.index.tolist()
    sides = ['up', 'down']

    num = max(len(offsets), len(sides))
    expt_map = pd.DataFrame(index = np.arange(num),
                        columns = cols + ['length','offset','side'])
    expt_map.loc[0, cols] = expt_map_constant
    expt_map.loc[0, 'offset'] = 0
    expt_map.loc[np.arange(len(lengths)), 'length'] = lengths
    expt_map.loc[np.arange(len(sides)),   'side']   = sides
    filename = os.path.join(saveDir, 'expt.length.map')
    expt_map.to_csv(filename, index=False, sep='\t')
    
    print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s "%(filename,
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     )
    
    sys.exit()  # and run make_library commands
    filenames = [os.path.join(saveDir, 'all.junctions_to_compare.library') ]
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
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'length'] == 10)]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'no_flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)
