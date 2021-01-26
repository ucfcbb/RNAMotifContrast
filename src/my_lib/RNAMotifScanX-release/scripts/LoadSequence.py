# Module for loading the PDB sequence

import re

def LoadSequenceWithTag(seq_file_name, tag):
	seq_mh = open(seq_file_name, 'r');
	# assume it is a multi-fasta file
	loaded_seqs = []	
	tag = tag.upper()
	to_include_seq = False
	for line in seq_mh:
		line = line.rstrip()
		line = line.upper()
		if(to_include_seq):
			loaded_seqs.append(line)
			to_include_seq = False
			continue
		if re.search('^>', line) and re.search(tag, line):
			loaded_seqs.append(line)
			to_include_seq = True
	return loaded_seqs
