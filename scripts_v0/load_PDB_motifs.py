import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import *
import rnamake.motif_library_sqlite
import rnamake.motif_type
import rnamake.motif_tree
import rnamake.cluster
import rnamake.motif
mlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(rnamake.motif_type.TWOWAY)
tlib = rnamake.motif_library_sqlite.MotifLibrarySqlite(rnamake.motif_type.TCONTACT)
mlib.load_all()
tlib.load_all()

# Notes:
# ./a.out -cseq (chip sequence) -css (secondary structure "dot bracket") -f for flow

motifs = []
m = [m for m in mlib.mdict.itervalues() if m.secondary_structure() == '(.(&).)'][0]
motifs.append(mlib.get_motif(m.name))

for i, m in enumerate(mlib.mdict.itervalues()):
	#print m.sequence()
	#print m.secondary_structure()
	motif_type = 'other'
	if m.secondary_structure() == '(.(&))':
		motif_type = '1x0'
	if m.secondary_structure() == '((&).)':
		motif_type = '0x1'
	if m.secondary_structure() == '(..(&))':
		motif_type = '2x0'
	if m.secondary_structure() == '((&)..)':
		motif_type = '0x2'
	if m.secondary_structure() == '((&)...)':
		motif_type = '0x3'
	if m.secondary_structure() == '(...(&))':
		motif_type = '3x0'
	if m.secondary_structure() == '(.(&).)':
		motif_type = '1x1'
	if m.secondary_structure() == '(..(&)..)':
		motif_type = '2x2'
	if m.secondary_structure() == '(.(&)..)':
		motif_type = '1x2'
	if m.secondary_structure() == '(..(&).)':
		motif_type = '2x1'
	if m.secondary_structure() == '(...(&).)':
		motif_type = '3x1'
	if m.secondary_structure() == '(.(&)...)':
		motif_type = '1x3'
	if m.secondary_structure() == '(...(&)..)':
		motif_type = '3x2'
	if m.secondary_structure() == '(..(&)...)':
		motif_type = '2x3'
	if m.secondary_structure() == '(...(&)...)':
		motif_type = '3x3'
		
	#mt = rnamake.motif_tree.MotifTree()
	#mt.add_motif(m)
	#mt.nodes[1].motif.to_pdb('%s.%s.pdb'%(m.name, motif_type))
	#if motif_type == '1x1':
	#	motifs.append(mt.nodes[1].motif)
	#mt.remove_node_level()

	print '\t'.join([m.name] + m.sequence().split('&') + [motif_type])