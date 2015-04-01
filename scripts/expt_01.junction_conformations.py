from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np

# functions
def getAllJunctionSeqs():
    # save junctions : first B1 junctions
    flanks = [['G', 'C', 'G', 'C'], ['C', 'U', 'A', 'G']]
    junctionMotifs = []
    junctionSeqs = {} 
    for motif in ['B1', 'B2', 'B1,B1', 'B1,B1,B1']:
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
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif])
    
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
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif])
            
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
            junctionSeqs[actual_motif] = pd.concat(junctionSeqs[actual_motif])
    
    return pd.concat(junctionSeqs)

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/junction_conformations/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files
    lengths = [8, 9, 10, 11]
    offsetsPerLength = {8:[0], 9:[-1,0,1], 10:[-1,0,1], 11:[0] }
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
        
    for length in lengths:

        print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s -out %s"%(filenames[length],
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.length_%d'%length))
        
