## for junction conformation expt, do two flanking base pairs, 7 different positions

saveDir = '150311_library_v2/junction_conformations/'
if not os.path.exists(saveDir): os.mkdir(saveDir)


lengthsDict = {8:[-1,0,1], 9:[-1,0, 1], 10:[-1, 0, 1], 11:[-1,0,1]}
junctions = ['B1', 'B1,B1', 'B1,B1,B1', 'M', 'M,M', 'B1,M', 'B1,B1,M']
junctions = ['B1']
flanking = [['GC', 'CG', 'AU', 'UA'], ['GC', 'CG', 'AU', 'UA']]
junctions_w_flank = [','.join(x) for x in itertools.product(flanking[0], junctions, flanking[1])]
filenames = []


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
