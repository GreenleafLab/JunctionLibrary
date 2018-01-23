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
    flanks = [['G', 'C', 'G', 'C'], ['C', 'U', 'A', 'G']]
    junctionMotifs = []
    junctionSeqs = {} 
    for motif in ['B1', 'B2', 'B1,B1', 'B1,B1,B1', 'W']:
        junctionSeqs[motif] = {}
        for flank in flanks:
            baseNum = len(flank)/2
            junctionMotif = ','.join(flank[:baseNum] + [motif] + flank[baseNum:])
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = baseNum
    
            if motif == 'B1,B1,B1':
                
                # only take those sequences with an A in the middle
                indices = []
                for index in junctionSeq.index:
                    if junctionSeq.loc[index, 'side1'][baseNum+1] == 'A':
                        indices.append(index)
                
                junctionSeq = junctionSeq.loc[indices]
            junctionSeqs[motif][''.join(flank)] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])
    
    # now 2x2s:
    flanks = [['G','C'], ['C', 'G']]
    for motif in ['N,N']:
        junctionSeqs[motif] = {}
        for flank in flanks:
            baseNum = len(flank)/2
            junctionMotif = ','.join(flank[:baseNum] + [motif] + flank[baseNum:])
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = baseNum
            junctionSeqs[motif][''.join(flank)] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])
    
    # now three by threes: take 16 2x2s and mutat other base pairs
    flanks = [['G','C'], ['C', 'G']]
    motif = 'M,M,M'
    junctionSeqs[motif] = {}
    for flank in flanks:
        
        junctionMotif = 'M,M'
        baseNum = len(flank)/2
        junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
        junctionSeq = junctionSeq.iloc[np.linspace(0, len(junctionSeq)-1, 12)]
        
        # first base, side 1
        junctionSeqsMut = {}
        #junctionSeqsMut['side1,b1'] = []
        for base in hjh.mutations.singles(flank[0]):
            if base != flank[0]:
                junctionSeqsMut['side1:%s._'%base] = pd.DataFrame(index=junctionSeq.index, columns=['side1', 'side2'])
                junctionSeqsMut['side1:%s._'%base].loc[:, 'side1'] = base + junctionSeq.loc[:, 'side1'] + flank[-1]
                junctionSeqsMut['side1:%s._'%base].loc[:, 'side2'] = hjh.mutations.complement(flank[-1]) + junctionSeq.loc[:, 'side2'] + hjh.mutations.complement(flank[0])
                
        for base in hjh.mutations.singles(flank[-1]):
            if base != flank[-1]:
                junctionSeqsMut['side1:_.%s'%base] = pd.DataFrame(index=junctionSeq.index, columns=['side1', 'side2'])
                junctionSeqsMut['side1:_.%s'%base].loc[:, 'side1'] = flank[0] + junctionSeq.loc[:, 'side1'] + base
                junctionSeqsMut['side1:_.%s'%base].loc[:, 'side2'] = hjh.mutations.complement(flank[-1]) + junctionSeq.loc[:, 'side2'] + hjh.mutations.complement(flank[0])

        for base in hjh.mutations.singles(hjh.mutations.complement(flank[-1])):
            if base != hjh.mutations.complement(flank[-1]):
                junctionSeqsMut['side2:%s._'%base] = pd.DataFrame(index=junctionSeq.index, columns=['side1', 'side2'])
                junctionSeqsMut['side2:%s._'%base].loc[:, 'side1'] = flank[0] + junctionSeq.loc[:, 'side1'] + flank[-1]
                junctionSeqsMut['side2:%s._'%base].loc[:, 'side2'] = base + junctionSeq.loc[:, 'side2'] + hjh.mutations.complement(flank[0])

        for base in hjh.mutations.singles(hjh.mutations.complement(flank[0])):
            if base != hjh.mutations.complement(flank[0]):
                junctionSeqsMut['side2:_.%s'%base] = pd.DataFrame(index=junctionSeq.index, columns=['side1', 'side2'])
                junctionSeqsMut['side2:_.%s'%base].loc[:, 'side1'] = flank[0] + junctionSeq.loc[:, 'side1'] + flank[-1]
                junctionSeqsMut['side2:_.%s'%base].loc[:, 'side2'] = hjh.mutations.complement(flank[-1]) + junctionSeq.loc[:, 'side2'] + base
                                                               
        junctionSeq = pd.concat(junctionSeqsMut, ignore_index=True)        
        junctionSeq.loc[:, 'n_flank'] = 0
        junctionSeqs[motif][''.join(flank)] = junctionSeq
    junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])

    # now 2x1s:
    flanks = [['G', 'C', 'G', 'C'], ['C', 'U', 'A', 'G']]
    for motif in ['B1', 'B1,B1']:
        for actual_motif in ['M,' + motif, motif+',M']:
            junctionSeqs[actual_motif] = {}
    
        for flank in flanks:
            # generate single mutants of central region of flank
            baseNum = 1
            junctionMotif = ','.join(flank[:baseNum] + ['M'] + [motif] + flank[baseNum+1:])
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = 1
            
            # take subset of B1,B1's
            if motif == 'B1,B1':
                indices = []
                for index in junctionSeq.index:
                    if junctionSeq.loc[index, 'side1'][2] == 'A':
                        indices.append(index)
                junctionSeq = junctionSeq.loc[indices] 
            
            # now choose subset where the Mismatched base side1 = flank[baseNum] side 1 or Mismatched base side2 = flank[baseNum] Complement 
            indices = []
            for index in junctionSeq.index:
                if junctionSeq.loc[index, 'side1'][baseNum] == flank[baseNum]:
                    indices.append(index)
                if junctionSeq.loc[index, 'side2'][-(baseNum+1)] == hjh.mutations.complement(flank[baseNum]):
                    indices.append(index)
            
            junctionSeqs['M,' + motif][''.join(flank)] = junctionSeq.loc[indices]      
            
            # also other side of junction
            baseNum = 2
            junctionMotif = ','.join(flank[:baseNum] + [motif] + ['M'] + flank[baseNum+1:])
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = 1
            
            # take subset of B1,B1's
            if motif == 'B1,B1':
                indices = []
                for index in junctionSeq.index:
                    if junctionSeq.loc[index, 'side1'][2] == 'A':
                        indices.append(index)
                junctionSeq = junctionSeq.loc[indices]   
                        
            indices = []
            for index in junctionSeq.index:
    
                if junctionSeq.loc[index, 'side2'][baseNum-1] == hjh.mutations.complement(flank[-baseNum]):
                    indices.append(index)
                if junctionSeq.loc[index, 'side1'][-baseNum] == flank[-baseNum]:
                    indices.append(index)
            
            junctionSeqs[motif+',M'][''.join(flank)]  = junctionSeq.loc[indices]
            
        for actual_motif in ['M,' + motif, motif+',M']:
            junctionSeqs[actual_motif] = pd.concat(junctionSeqs[actual_motif], names=['flank', 'junction_num'])
    
    return pd.concat(junctionSeqs, names=['junction'])

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/junction_conformations/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
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
    
    # make map files
    lengths = [8, 9, 10, 11, 12]
    offsetsPerLength = {8:[0], 9:[-1,0,1], 10:[-1,0,1], 11:[0], 12:[0] }
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
    
    filenames = {}
    for length in lengths:
        filenames[length] = os.path.join(saveDir, 'expt.length_%d.map'%length)  
        print "%%run ~/JunctionLibrary/make_library.py -map %s -jun %s -jl 6 -out %s"%(filenames[length],
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
    
    # find numbers again
    for junction in np.unique(allSeqs.loc[:, 'junction']):
        print '%s\t%d'%(junction, (allSeqs.loc[:, 'junction'] == junction).sum())
    
    # print varna sequences fro moving A bulge and think about it
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'flank'] == 'UGCGCU')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(
            name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'no_flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)
        
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'AA_')&(allSeqs.loc[:, 'flank'] == 'GCGC')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.png'%loc))
        os.system(cmnd)

    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'CAG_CA')&(allSeqs.loc[:, 'flank'] == 'UGCU')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.png'%loc))
        os.system(cmnd)
        
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'CAAG_CA')&(allSeqs.loc[:, 'flank'] == 'UGCU')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'no_flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)

    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'AA_AA')&(allSeqs.loc[:, 'flank'] == 'UGCU')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand(name=os.path.join(saveDir, 'junction_%d.%s.length_%d.png'%(loc, allSeqs.loc[loc, 'no_flank'], allSeqs.loc[loc, 'length'])))
        os.system(cmnd)