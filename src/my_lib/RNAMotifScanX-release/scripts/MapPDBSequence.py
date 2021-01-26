# Module for parsing MC-Annotate and RNAVIEW outputs

import re
from Bio import pairwise2


# The first column is the combined ID, the second column is the chain,
# the third column is the index, the fourth column is the nucleotide
def GetRecognizedNucleotides(mc_annotate_file): 
	mca_fh = open(mc_annotate_file, 'r')
	# skip the first line	
	temp = mca_fh.readline()
	residues = []
	for line in mca_fh:
		if re.search('^Adjacent\sstackings', line):
			break
		decom = re.split('\s+', line)
		if decom[2] == 'A' or decom[2] == 'C' or decom[2] == 'G' or decom[2] == 'U':
			if decom[0][0] == '\'':
				chain_and_pos = re.search('\'(\d)\'(\d+)', decom[0])
			else: 
				chain_and_pos = re.search('(\D)(\d+)', decom[0])			
			residues.append([decom[0], chain_and_pos.group(1), chain_and_pos.group(2), decom[2]])
	return residues

# Colapse the chains from the recognized residues
def GetRecognizedChains(recognized_residues):
	chain_hash = {}
	for k in recognized_residues:
		chain_hash[k[1]] = 1
	presented_chains = []
	for key in chain_hash:	
		presented_chains.append(key)		
	return presented_chains


# return the mappings between the input sequence and the recognized residues
def AlignToReference(chain, ref_sequence, recognized_residues):
	# concaternate the sequence, assume the residues are in order
	rec_sequence = ""
	residue_id = []	
	for residue_info in recognized_residues:
		if residue_info[1] == chain:
			rec_sequence += residue_info[3]
			residue_id.append(residue_info[0])	
		# ENDIF
	# ENDFOR
	# align the rec_sequence to the reference sequence
	alignment = pairwise2.align.globalms(ref_sequence, rec_sequence, 3, -100, -10, -2)
	# parse the alignment result
	matched_id = {}
	ref_accumulate = 0
	rec_accumulate = 0
	for i in range(0, len(alignment[0][0])):
		if alignment[0][0][i] != '-' and alignment[0][1][i] != '-':
			matched_id[residue_id[rec_accumulate]] = ref_accumulate
			ref_accumulate += 1
			rec_accumulate += 1
		elif alignment[0][0][i] == '-' and alignment[0][1][i] != '-':
			rec_accumulate += 1
		elif alignment[0][0][i] != '-' and alignment[0][1][i] == '-':
			ref_accumulate += 1
		# ENDIF
	# ENDFOR		
	return matched_id
