import sys
import math
import random
import os
import hjh.tecto_assemble
import pandas as pd
import numpy as np
import hjh.mutations

def calc_end_distro_score(selected_constructs):
    bps = ['GC','CG','AU','UA']

    step_distro = {}
    for bp1 in bps:
        for bp2 in bps:
            step_distro[bp1+"="+bp2] = 0

    for c in selected_constructs:
        for step in c.steps:
            if step not in step_distro:
                continue
            step_distro[step] += 1

    avg = 0
    std = 0
    for v in step_distro.itervalues():
        avg += v
    avg /= len(step_distro)

    for v in step_distro.itervalues():
        std += (v - avg) ** 2
    std /= len(step_distro)
    std = math.sqrt(std)

    return std


def calc_bin_score(selected_constructs):
    bins = [0 for x in range(51)]
    avg = 0
    for c in selected_constructs:
        bins[c.dG_bin] += 1

    std = 0
    avg = 0
    for b in bins:
        avg += b
    avg /= len(bins)
    for b in bins:
        std += (b - avg) ** 2
    std /= len(bins)
    std = math.sqrt(std)
    return std


class Construct(object):
    def __init__(self, seq, dG):
        self.seq = seq
        self.dG = dG
        self.steps = []
        self.dG_bin = 0
        seq1 = seq[10:20]
        seq2 = seq[26:36][::-1]
        bps = []
        for i in range(len(seq1)):
            bps.append(seq1[i]+seq2[i])
        steps = []
        for i in range(1, len(seq1)):
            steps.append(bps[i-1]+"="+bps[i])
        self.steps = steps

f = open(os.path.join(os.path.expanduser('~/JunctionLibrary/seq_params/'), 'exhustive_helices.results'))
lines = f.readlines()
f.close()

constructs = []
selected_constructs = []

for l in lines:
    spl = l.split()
    c = Construct(spl[0], float(spl[-1]))
    constructs.append(c)

min_dG = constructs[-1].dG
max_dG = constructs[0].dG

r = max_dG - min_dG
interval = abs(r) / 50

for c in constructs:
    c.dG_bin = (int)((min_dG - c.dG) / interval)

if len(sys.argv[1]) < 2:
    print "need to set a number for the number of constructs you need"
    sys.exit()

selected_num = int(sys.argv[1])
selected_interval = len(constructs) / selected_num
current = 0
while current < len(constructs):
    selected_constructs.append(constructs[current])
    current += selected_interval

distro_score = calc_end_distro_score(selected_constructs)
bin_score = calc_bin_score(selected_constructs)
score = distro_score + 5*bin_score
new_score = 1000
new_selected_constructs = []

for i in range (100000):

    new_selected_constructs = selected_constructs[:]
    rand_pos = random.randint(0, len(new_selected_constructs)-1)
    new_c = random.choice(constructs)
    while new_c in new_selected_constructs:
        new_c = random.choice(constructs)
    new_selected_constructs[rand_pos] = new_c
    distro_score = calc_end_distro_score(new_selected_constructs)
    bin_score = calc_bin_score(new_selected_constructs)
    new_score = distro_score + 5*bin_score
    if new_score < score:
        selected_constructs = new_selected_constructs
        score = new_score
        print "SCORE", score
        
saveDir = '150311_library_v2/sequences'
if not os.path.exists(saveDir): os.mkdir(saveDir)

seq_params_file = os.path.expanduser('~/JunctionLibrary/seq_params/seq_params.txt')
seq_params    = pd.read_table(seq_params_file , index_col=[0,1])

params = pd.Series({'helix':'joe',
'junction':'_',
'receptor':       '11nt',
'loop':           'GGAA',
'base':         'normal',
'adapters':     'truseq',
'length':           10,
'offset':            0,
'side':             'up'})
tectoSeq = hjh.tecto_assemble.TectoSeq(seq_params, params)
allSeqs_10bp = pd.DataFrame(index=np.arange(selected_num), columns=tectoSeq.params.index)
for i, c in enumerate(selected_constructs):
    
    allSeqs_10bp.loc[i, params.index] = params
    allSeqs_10bp.loc[i, 'helix'] += '_%d'%i
    allSeqs_10bp.loc[i, 'tecto_sequence'] = c.seq.replace('T', 'U')
    seq = seq_params.loc[('adapters', params.loc['adapters'])]
    allSeqs_10bp.loc[i, 'sequence'] = seq.loc['side1'] +  c.seq.replace('U', 'T') + seq.loc['side2']
allSeqs = allSeqs_10bp

num_9bp = 100
allSeqs_9bp = pd.DataFrame(index=allSeqs_10bp.index[np.linspace(0, len(allSeqs_10bp)-1,num_9bp).astype(int)] , columns=tectoSeq.params.index)
for i in allSeqs_9bp.index:
    allSeqs_9bp.loc[i, params.index] = params
    allSeqs_9bp.loc[i, 'length'] = 9
    allSeqs_9bp.loc[i, 'helix'] += '_%d'%i
    old_seq = list(allSeqs_10bp.loc[i, 'tecto_sequence'])
    split_inds = [15, 30]
    new_seq = old_seq[:split_inds[0]] + old_seq[split_inds[0]+1:split_inds[1]] + old_seq[split_inds[1]+1:]
    allSeqs_9bp.loc[i, 'tecto_sequence'] = ''.join(new_seq).replace('T', 'U')
    seq = seq_params.loc[('adapters', params.loc['adapters'])]
    allSeqs_9bp.loc[i, 'sequence'] = seq.loc['side1'] +  allSeqs_9bp.loc[i, 'tecto_sequence'].replace('U', 'T') + seq.loc['side2']
allSeqs = pd.concat([allSeqs, allSeqs_9bp])

num_11bp = 100
new_bp = 'U'
allSeqs_11bp = pd.DataFrame(index=allSeqs_10bp.index[np.linspace(0, len(allSeqs_10bp)-1,num_11bp).astype(int)] , columns=tectoSeq.params.index)
for i in allSeqs_11bp.index:
    allSeqs_11bp.loc[i, params.index] = params
    allSeqs_11bp.loc[i, 'length'] = 11
    allSeqs_11bp.loc[i, 'helix'] += '_%d'%i
    old_seq = list(allSeqs_10bp.loc[i, 'tecto_sequence'])
    split_inds = [15, 30]
    new_seq = old_seq[:split_inds[0]] + [new_bp] + old_seq[split_inds[0]:split_inds[1]] + [hjh.mutations.complement(new_bp)] + old_seq[split_inds[1]:]
    allSeqs_11bp.loc[i, 'tecto_sequence'] = ''.join(new_seq).replace('T', 'U')
    seq = seq_params.loc[('adapters', params.loc['adapters'])]
    allSeqs_11bp.loc[i, 'sequence'] = seq.loc['side1'] + allSeqs_11bp.loc[i, 'tecto_sequence'].replace('U', 'T') + seq.loc['side2']
allSeqs = pd.concat([allSeqs, allSeqs_11bp])

num_11bp = 100
new_bp = 'G'
allSeqs_11bp = pd.DataFrame(index=allSeqs_10bp.index[np.linspace(0, len(allSeqs_10bp)-1,num_11bp).astype(int)] , columns=tectoSeq.params.index)
for i in allSeqs_11bp.index:
    allSeqs_11bp.loc[i, params.index] = params
    allSeqs_11bp.loc[i, 'length'] = 11
    allSeqs_11bp.loc[i, 'helix'] += '_%d'%i
    old_seq = list(allSeqs_10bp.loc[i, 'tecto_sequence'])
    split_inds = [15, 30]
    new_seq = old_seq[:split_inds[0]] + [new_bp] + old_seq[split_inds[0]:split_inds[1]] + [hjh.mutations.complement(new_bp)] + old_seq[split_inds[1]:]
    allSeqs_11bp.loc[i, 'tecto_sequence'] = ''.join(new_seq).replace('T', 'U')
    seq = seq_params.loc[('adapters', params.loc['adapters'])]
    allSeqs_11bp.loc[i, 'sequence'] = seq.loc['side1'] + allSeqs_11bp.loc[i, 'tecto_sequence'].replace('U', 'T') + seq.loc['side2']
allSeqs = pd.concat([allSeqs, allSeqs_11bp])

out_file = os.path.join(saveDir, 'all.sequences.library')
   
neworder = [ 'junction', 'length', 'offset', 'helix', 'receptor', 'loop',  'side',
            'flank', 'no_flank', 'n_flank', 'junction_seq',  'junction_SS',  'ss_correct',
            'adapters', 'base', 'helix_one_length', 'effective_length', 'helix_seq',
            'tecto_sequence', 'sequence', 'ss']
allSeqs.loc[:, neworder].to_csv(out_file+'.txt', sep='\t')



