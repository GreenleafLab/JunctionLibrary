from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np
import sys

def getAllJunctions():
    pdb_junctions = pd.read_table('150311_library_v2/ires_junctions.txt', index_col=0, header=0, names=['junction', 'side1', 'side2'])
    pdb_junctions.loc[:, 'n_flank'] = 1
    return pdb_junctions


## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/ires_junctions/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctions()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files
    lengths = [8, 9, 10, 11, 12]
    offsetsPerLength = {8:[0], 9:[-1,0,1], 10:[-1,0,1], 11:[0], 12:[0] }
    filenames = {}
    num_lengths = len(lengths)
    expt_map_constant = pd.Series(index = ['helix','junction','receptor','loop', 'base',   'adapters'],
                                  data  = ['wc',   'defunct', '11nt',    'GGAA', 'normal', 'truseq'], dtype=str)
    cols = expt_map_constant.index.tolist()
    sides = ['up', 'down']
    for length in lengths:
        offsets = offsetsPerLength[length]
        num = max(len(offsets), len(sides))
        expt_map = pd.DataFrame(index = np.arange(num),
                            columns = cols + ['length','offset','side'])
        expt_map.loc[0, cols] = expt_map_constant
        expt_map.loc[0, 'length'] = length
        expt_map.loc[np.arange(len(offsets)), 'offset'] = offsets
        expt_map.loc[np.arange(len(sides)), 'side'] = sides
        filenames[length] = os.path.join(saveDir, 'expt.length_%d.map'%length)
        expt_map.to_csv(filenames[length], index=False, sep='\t')
    
    filenames = {}
    for length in lengths:
        filenames[length] = os.path.join(saveDir, 'expt.length_%d.map'%length)  
        print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s -out %s"%(filenames[length],
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.length_%d'%length))
    
    sys.exit()  # and run make_library commands
    
    saveDirs = ['150311_library_v2/ires_%s/'%s for s in ['junctions', 'singles', 'doubles']]
    
    filenames = [os.path.join(saveDirs[0], 'all.junctions_to_compare.length_%d'%length) for length in lengths]
    for i in [1,2]:
        filenames += [os.path.join(saveDirs[i], 'all.junctions_to_compare.library')]
        
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
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'flank'] == 'GCGC')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand()[0]
        os.system(cmnd.replace('test.png', 'junction_%d.png'%loc))