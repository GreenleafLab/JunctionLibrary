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
    flanks = [['G', 'C', 'G', 'C'],  ['C', 'U', 'A', 'G']]

    for motif in ['B1', '_']:
        junctionSeqs[motif] = {}
        
        for flank in flanks:

            baseNum = len(flank)/2
            junctionMotif = ','.join(flank[:baseNum] + [motif] + flank[baseNum:])
    
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = baseNum
    
            if motif == 'B1':
                # only take those sequences with an A in the middle
                indices = []
                for index in junctionSeq.index:
                    if junctionSeq.loc[index, 'side1'][baseNum] == 'A':
                        indices.append(index)
                junctionSeq = junctionSeq.loc[indices]
    
            junctionSeqs[motif][''.join(flank)] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank'])
    
    return pd.concat(junctionSeqs, names=['junction'])

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/loop_flanking/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # load seq params
    seq_params_file = os.path.expanduser('~/JunctionLibrary/seq_params/seq_params.txt')
    seq_params = pd.read_table(seq_params_file, index_col=[0,1])
    
    # take WC and independently generate seqs from wc
    seq = seq_params.loc[('helix', 'wc'), ['side1', 'side2']]
    
    # make secnond to last bases all and GU, last base any
    lflank = Junction(tuple(['W+'])).sequences
    lflank2 = Junction(tuple(['W'])).sequences
    
    iterables = [lflank.index, lflank2.index]
    new_seq = pd.DataFrame(index=pd.MultiIndex.from_product(iterables, names=['base-2', 'base-1']), columns=['side1', 'side2'])
    for ind1, ind2 in itertools.product(lflank.index, lflank2.index):
        new_seq.loc[(ind1, ind2), 'side1'] = seq.loc['side1'][:-2] + lflank.loc[ind1, 'side1'] + lflank2.loc[ind2, 'side1']
        new_seq.loc[(ind1, ind2), 'side2'] = lflank2.loc[ind2, 'side2'] + lflank.loc[ind1, 'side2'] + seq.loc['side2'][2:]
    
    lindex=['lflank_' + new_seq.loc[(ind1, ind2), 'side1'][-2:]+'_'+new_seq.loc[(ind1, ind2), 'side2'][:2] for ind1, ind2 in itertools.product(lflank.index, lflank2.index)]
    print pd.DataFrame(index=lindex, data=new_seq.values, columns=['side1', 'side2'])
    
    # now do the same for the receptor
    rflank = Junction(tuple(['W'])).sequences
    rflank2 = Junction(tuple(['W+'])).sequences
    iterables = [rflank.index, rflank2.index]
    new_seq = pd.DataFrame(index=pd.MultiIndex.from_product(iterables, names=['base1', 'base2']), columns=['side1', 'side2'])
    for ind1, ind2 in itertools.product(rflank.index, rflank2.index):
        print ind1, ind2
        new_seq.loc[(ind1, ind2), 'side1'] = rflank.loc[ind1, 'side1'] + rflank2.loc[ind2, 'side1'] + seq.loc['side1'][2:]
        new_seq.loc[(ind1, ind2), 'side2'] = seq.loc['side2'][:-2] + rflank2.loc[ind2, 'side2'] + rflank.loc[ind1, 'side2']
    rindex=['rflank_' + new_seq.loc[(ind1, ind2), 'side1'][:2]+'_'+new_seq.loc[(ind1, ind2), 'side2'][-2:] for ind1, ind2 in itertools.product(rflank.index, rflank2.index)]
    print pd.DataFrame(index=rindex, data=new_seq.values, columns=['side1', 'side2'])
                
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files

    expt_map_constant = pd.Series(index = ['junction','receptor','loop', 'base',   'adapters', 'side', 'length', 'offset'],
                                  data  = ['defunct', '11nt',    'GGAA', 'normal', 'truseq', 'up', 10, 0], dtype=str)
    cols = expt_map_constant.index.tolist()

    num = len(rindex+lindex)
    expt_map = pd.DataFrame(index = np.arange(num),
                        columns = cols + ['helix'])
    expt_map.loc[0, cols] = expt_map_constant
    expt_map.loc[:, 'helix'] = rindex+lindex
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
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == '_')]
    # print varna sequences fro moving A bulge and think about it
    #subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'length'] == 10)]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'no_flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)