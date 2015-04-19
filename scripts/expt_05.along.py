from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np
import sys

# functions
def getAllJunctionSeqs():
    # save junctions : first B1 junctions
    junctionMotifs = []
    junctionSeqs = {}

    flanks = [['G','C'], ['C', 'G']]
    for motif in ['N']:
        junctionSeqs[motif] = {}
        for flank in flanks:
            baseNum = len(flank)/2
            junctionMotif = ','.join(flank[:baseNum] + [motif] + flank[baseNum:])
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = baseNum
            junctionSeqs[motif][''.join(flank)] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])
        
    motif = 'B1'
    junctionSeqs[motif] = {}
    for flank in flanks:
        junctionSeq = pd.DataFrame(index=np.arange(2), columns = ['side1', 'side2', 'n_flank'])
        junctionSeq.loc[:, 'n_flank'] = baseNum
        for i, base in enumerate(['A', 'U']):
            junctionSeq.loc[i, 'side1'] = flank[0] + base + flank[-1]
            junctionSeq.loc[i, 'side2'] = hjh.mutations.complement(flank[-1]) + hjh.mutations.complement(flank[0])
        junctionSeqs[motif][''.join(flank)] = junctionSeq
    junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])

    motif = 'B1,B1'
    junctionSeqs[motif] = {}
    for flank in flanks:
        junctionSeq = pd.DataFrame(index=np.arange(2), columns = ['side1', 'side2', 'n_flank'])
        junctionSeq.loc[:, 'n_flank'] = baseNum
        for i, base in enumerate(['AA', 'UU']):
            junctionSeq.loc[i, 'side1'] = flank[0] + base + flank[-1]
            junctionSeq.loc[i, 'side2'] = hjh.mutations.complement(flank[-1]) + hjh.mutations.complement(flank[0])
        junctionSeqs[motif][''.join(flank)] = junctionSeq
    junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])

    motif = 'B1,B1'
    junctionSeqs[motif] = {}
    for flank in flanks:
        junctionSeq = pd.DataFrame(index=np.arange(2), columns = ['side1', 'side2', 'n_flank'])
        junctionSeq.loc[:, 'n_flank'] = baseNum
        for i, base in enumerate(['AAA', 'UUU']):
            junctionSeq.loc[i, 'side1'] = flank[0] + base + flank[-1]
            junctionSeq.loc[i, 'side2'] = hjh.mutations.complement(flank[-1]) + hjh.mutations.complement(flank[0])
        junctionSeqs[motif][''.join(flank)] = junctionSeq
    junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])


    
    return pd.concat(junctionSeqs, names=['junction'])

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/along_N'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs = {}
    junctionSeqs['up'] = getAllJunctionSeqs()
    junctionSeqs['down'] = junctionSeqs['up'].copy()
    junctionSeqs['down'].loc[:, 'side1'] = junctionSeqs['up'].loc[:, 'side2']
    junctionSeqs['down'].loc[:, 'side2'] = junctionSeqs['up'].loc[:, 'side1']
    junctionSeqs = pd.concat(junctionSeqs)
    # switch side before adding flank
    
    # add a flanker
    junctionSeqs.loc[:, 'side1'] = junctionSeqs.loc[:, 'side1'] + 'U'
    junctionSeqs.loc[:, 'side2'] = 'A' + junctionSeqs.loc[:, 'side2']
    junctionSeqs.loc[:, 'side1'] = 'U' + junctionSeqs.loc[:, 'side1']
    junctionSeqs.loc[:, 'side2'] = junctionSeqs.loc[:, 'side2'] + 'A'
    junctionSeqs.loc[:, 'n_flank'] = junctionSeqs.loc[:, 'n_flank'] + 1
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # load length/offset combos
    lengths = [8, 9, 10, 11]
    offsetsPerLength = {8:range(-2,3), 9:range(-2,4), 10:range(-3,4), 11:range(-3, 5)}
    filenames = {}
    num_lengths = len(lengths)
    expt_map_constant = pd.Series(index = ['helix','junction','receptor','loop', 'base',   'adapters'],
                                  data  = ['wc',   'defunct', '11nt',    'GGAA', 'normal', 'truseq'], dtype=str)
    cols = expt_map_constant.index.tolist()
    sides = ['up']
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
        print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s -out %s"%(filenames[length],
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.length_%d'%length))
    
    
    sys.exit()  # and run make_library commands
    
    
    filenames = [os.path.join(saveDir, 'all.junctions_to_compare.length_%d'%length) for length in lengths]
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
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'flank'] == 'UCGU')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'no_flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)