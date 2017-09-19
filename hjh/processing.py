import os
import pandas as pd
import subprocess

def thread_together(inside, outside):
    return (outside[0] + inside + outside[1])

def load_file(filename):
    ext = os.path.splitext(filename)[1]
    if ext == '.csv':
        return pd.read_csv(filename)
    elif ext == '.dat':
        return pd.read_table(filename)
  


#### THREEWAYJUNCTION FUNCTIONS ####
def convert_threeway_to_twoway(junction_seq, stem_loop_seq, position):
    """Take a three way junction and add a stem loop to one of the junciton positions to make a two-way junction.
    
    'junction_seq' gives the sequence of the 3way junction, e.g.: ['GAAC','GU','AC'] (2 A's inserted at position 1; flanking base pairs GC;CG;UA
    'stem-loop_seq' gives the sequence of the stem loop, i.e. 'CCUCGGGAACGAGG' (5 bp stem helix and GGAA tetraloop)
    'position' indicates between which two parts of the junction shoudl have the stem loop inserted (i.e. (2_0, 0_1, or 1_2))
    
    Returns:
    the sequence of the resulting two-way junction, oriented so that the stem-loop helix is on side 1. 
    """
    
    seqs = junction_seq
    
    if position=='2_0':
        two_way_seq = seqs[2] + stem_loop_seq + seqs[0] + '_' +  seqs[1]
    elif position == '0_1':
        two_way_seq = seqs[0] + stem_loop_seq + seqs[1] + '_' + seqs[2]
    elif position == '1_2':
        two_way_seq = seqs[1] + stem_loop_seq + seqs[2] + '_' + seqs[0]
    
    return two_way_seq

def convert_threeway_to_single(junction_seqs, stem_loops):
    """Take a three way junction and add a stem loop to one of the junciton positions to make a two-way junction.
    
    'junction_seq' gives the sequence of the 3way junction, e.g.: ['GAAC','GU','AC'] (2 A's inserted at position 1; flanking base pairs GC;CG;UA
    'stem-loop_seq' gives the sequence of the two stem loops, i.e. ['CCUCGGGAACGAGG', 'CCUCGUUCGCGAGG'] (5 bp stem helix and GGAA tetraloop)

    Returns:
    
    """
    return junction_seqs[0] + stem_loops[0] + junction_seqs[1] + stem_loops[1] + junction_seqs[2]


def convert_threeway_helix_to_stemloops(helix_row, loop_seqs):
    """Take a helix Series with fields ['base_side1', 'base_side2', 'h1_side1',
    'h1_side2', 'h2_side1', 'h2_side2'] and a list of two loop seqs.
    Returns:
    base_helix, stem_loop1, stem_loop2
    """
    
    base_helix = pd.Series({'side1':helix_row.base_side1, 'side2':helix_row.base_side2})
    stem_loop1 = helix_row.loc['h1_side1'] + loop_seqs[0] +  helix_row.loc['h1_side2']
    
    stem_loop2 = helix_row.loc['h2_side1'] + loop_seqs[1] +  helix_row.loc['h2_side2']
    
    return base_helix, [stem_loop1, stem_loop2]

def convert_junction_to_seq(junction_row, permutation=0):
    """Take a junction Series with fields ['j1', 'j2', 'j3'] and return a list
    of the sequences. Optional, circularly permutate the sequence. (permute=[1,2])"""
    
    if permutation == 0:
        fields = ['j1', 'j2', 'j3']
    elif permutation == 1:
        fields = ['j3', 'j1', 'j2']
    elif permutation == 2:
        fields = ['j2', 'j3', 'j1']
    return [junction_row.loc[i] for i in fields]

def return_junction_fields():
    return ['j1', 'j2', 'j3']

def return_threeway_helix_fields():
    return ['base_side1', 'base_side2', 'h1_side1',
    'h1_side2', 'h2_side1', 'h2_side2']


def check_ss_structure_set(seqs):
    """Get the dotbracket ss structure of a set of sequences."""
    seqs.seq.to_csv('temp.dat', sep='\t', header=False)
    s = subprocess.check_output('cat temp.dat | awk  -F "\\t" \'{print ">"$1"\\n"$2}\' | RNAfold --noPS', shell=True).strip().split('\n')
    ss_all = {}
    for idx_almost, ss_almost in zip(s[::3], s[2::3]):
        if idx_almost[0]!= '>':
            print "error: every third entry is not the fasta header. "
            sys.exit()
        idx = idx_almost[1:]
        try:
            idx = int(idx)
        except ValueError:
            pass
        ss = ss_almost.split()[0]
        ss_all[idx] = ss
    return pd.Series(ss_all).loc[seqs.index]
        
