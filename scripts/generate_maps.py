from hjh.junction import Junction
import itertools
## for junction conformation expt, do two flanking base pairs, 7 different positions

saveDir = '150311_library_v2/junction_conformations/'
if not os.path.exists(saveDir): os.mkdir(saveDir)

# make other parameters
lengths = [8, 9, 10, 11]
offsets = [-1, 0, 1]
junctions = ['M,M']
flanking = [['GC', 'CG', 'AU', 'UA'], ['GC', 'CG', 'AU', 'UA']]
junctions_w_flank = [','.join(x) for x in itertools.product(['W'], junctions, ['W'])]

num_junctions = 1
num_offsets = len(offsets)
num_lengths = len(lengths)
expt_map = pd.DataFrame(index = np.arange(np.max([num_offsets, num_junctions, num_lengths])),
                        columns = ['helix','junction','length','receptor','loop','offset','base','adapters', 'side'])
expt_map.loc[0, 'helix'] = 'wc'
expt_map.loc[0, 'adapters'] = 'truseq'
expt_map.loc[0, 'base'] = 'normal'
expt_map.loc[0, 'receptor'] = '11nt'
expt_map.loc[0, 'loop'] = 'GGAA'
expt_map.loc[0, 'junction'] = 'user_defined'
expt_map.loc[np.arange(num_lengths), 'length'] = lengths
expt_map.loc[np.arange(num_offsets), 'offset'] = offsets
expt_map.loc[[0,1], 'side'] = ['up', 'down']

filename = os.path.join(saveDir, 'expt.all.map')
expt_map.to_csv(filename, sep='\t', index=False)

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




junctionSeqs = pd.concat(junctionSeqs)
            
junctionMotifs = ['G,C,B1,G,C', ]
junctionSeqs = pd.DataFrame(columns=['side1', 'side2'])
for motif in junctionMotifs:
    junctionSeqs = junctionSeqs.append(Junction(tuple(motif.split(','))).sequences, ignore_index=True)
junctionSeqs.to_csv(os.path.join(saveDir, 'expt.B1B1.junctions'), sep='\t', index=True)
print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s"%(filename, os.path.join(saveDir, 'expt.B1.junctions'))

junctionMotifs = []
flanking1 = ['UA,GC', 'CG,UA']
flanking2 = ['GC,AU', 'AU,CG']
for x in itertools.product(flanking1, ['B1', 'B1,B1', 'B1,B1']):
    junctionMotifs.append(','.join(x))
    



for length in lengthsDict.keys():
    num_junctions = len(junctions_w_flank)
    num_offsets = len(offsets)

    expt_map = pd.DataFrame(index = np.arange(max(num_offsets, num_junctions)),
                            columns = ['helix','junction','length','receptor','loop','offset','base','adapters', 'side'])
    expt_map.loc[0, 'helix'] = 'wc'
    expt_map.loc[0, 'adapters'] = 'truseq'
    expt_map.loc[0, 'base'] = 'normal'
    expt_map.loc[0, 'receptor'] = '11nt'
    expt_map.loc[0, 'loop'] = 'GGAA'
    expt_map.loc[0, 'length'] = length
    expt_map.loc[np.arange(num_junctions), 'junction'] = junctions_w_flank
    expt_map.loc[np.arange(num_offsets), 'offset'] = offset
    expt_map.loc[[0,1], 'side'] = ['up', 'down']
    print expt_map
    filename = os.path.join(saveDir, 'expt.length_%s.offset_%d.map'%(length, offset))
    filenames.append(filename)
    expt_map.to_csv(filename, sep='\t', index=False)

for filename in filenames:
    print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -swap none"%(filename)
    print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -swap swap"%(filename)
    print "nohup python ~/JunctionLibrary/scripts/make_library.py -map %s >> %s.out &"%(filename, os.path.basename(filename))
    
lengths = [8, 9, 10, 11]
offsets = [-1, 0, 1]
junctions = ['M,M']
flanking = [['GC', 'CG', 'AU', 'UA'], ['GC', 'CG', 'AU', 'UA']]
junctions_w_flank = [','.join(x) for x in itertools.product(flanking[0], junctions, flanking[1])]

num_junctions = len(junctions_w_flank)
num_offsets = len(offsets)
num_lengths = len(lengths)
expt_map = pd.DataFrame(index = np.arange(max(num_offsets, num_junctions)),
                        columns = ['helix','junction','length','receptor','loop','offset','base','adapters', 'side'])
expt_map.loc[0, 'helix'] = 'wc'
expt_map.loc[0, 'adapters'] = 'truseq'
expt_map.loc[0, 'base'] = 'normal'
expt_map.loc[0, 'receptor'] = '11nt'
expt_map.loc[0, 'loop'] = 'GGAA'
expt_map.loc[np.arange(num_lengths), 'length'] = lengths
expt_map.loc[np.arange(num_junctions), 'junction'] = junctions_w_flank
expt_map.loc[np.arange(num_offsets), 'offset'] = offsets
expt_map.loc[[0,1], 'side'] = ['up', 'down']
print expt_map
filename = os.path.join(saveDir, 'expt.all_%s.map'%(''.join(junctions)).replace(',', ''))
expt_map.to_csv(filename, sep='\t', index=False)

# generate sequences
