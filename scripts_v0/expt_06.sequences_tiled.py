from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np
import sys
import seqfun

# functions
def getAllJunctionSeqs():
    # save junctions : first B1 junctions
    junctionMotifs = []
    junctionSeqs = {}
    flanks = ['']

    for motif in ['_']:
        junctionSeqs[motif] = {}
        
        for flank in flanks:
    
            junctionSeq = pd.DataFrame(index=np.arange(16 + 64), columns=['side1', 'side2', 'n_flank'])
            junctionSeq.loc[:, 'n_flank'] = 0
            count = 0
            for base1, base2 in itertools.product('ACGU', 'ACGU'):
                junctionSeq.loc[count, 'side1'] = ''.join([base1,base2]*3)
                junctionSeq.loc[count, 'side2'] = seqfun.reverseComplement(junctionSeq.loc[count, 'side1'])
                count+=1
            for base1, base2, base3 in itertools.product('ACGU', 'ACGU', 'AGCU'):
                junctionSeq.loc[count, 'side1'] = ''.join([base1,base2, base3]*2)
                junctionSeq.loc[count, 'side2'] = seqfun.reverseComplement(junctionSeq.loc[count, 'side1'])
                count+=1
    
            junctionSeqs[motif][''.join(flank)] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank'])
    
    return pd.concat(junctionSeqs, names=['junction'])

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/sequences_tiled/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files

    expt_map_constant = pd.Series(index = ['junction','receptor','loop', 'base',   'adapters', 'side', 'helix', 'offset'],
                                  data  = ['defunct', '11nt',    'GGAA', 'normal', 'truseq', 'up', 'wc', 0], dtype=str)
    cols = expt_map_constant.index.tolist()
    
    lengths = [8,9,10,11,12]
    num = len(lengths)
    expt_map = pd.DataFrame(index = np.arange(num),
                        columns = cols + ['length'])
    expt_map.loc[0, cols] = expt_map_constant
    expt_map.loc[:, 'length'] = lengths
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
    subset = allSeqs.loc[(allSeqs.loc[:, 'flank'] == 'AAAAAA')]
    # print varna sequences fro moving A bulge and think about it
    #subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'length'] == 10)]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)