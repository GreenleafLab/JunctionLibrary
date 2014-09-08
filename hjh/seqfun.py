def reverseComplement(seq):
    for base in seq:
        if base not in 'ATCGatcg':
            print "TypeError: NOT a DNA sequence"
            return None
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])