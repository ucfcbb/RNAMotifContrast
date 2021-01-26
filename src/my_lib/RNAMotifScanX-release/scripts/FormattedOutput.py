import re

def WriteInteractions(pdb_code, chain, ref_seq, interactions, nucleotide_hash, out_file):
	out_fh = open(out_file, 'w')
	out_fh.write(">%s_%s\n" %(pdb_code, chain))
	out_fh.write("%s\n" %ref_seq)
	out_fh.write("#info=basepair\n")
	# print basepairs
	for it in interactions:
		if len(it) == 4 and it[0] in nucleotide_hash and it[1] in nucleotide_hash:
			i = nucleotide_hash[it[0]]
			j = nucleotide_hash[it[1]]
			out_fh.write("%s-%s,%s,%s,%s-%s\n" %(i, j, it[2], it[3], it[0], it[1]))
	# print base stackings
	out_fh.write("#info=stacking\n")
	for it in interactions:
		if len(it) == 3 and it[0] in nucleotide_hash and it[1] in nucleotide_hash:
			i = nucleotide_hash[it[0]]
			j = nucleotide_hash[it[1]]
			out_fh.write("%s-%s,%s,%s-%s\n" %(i, j, it[2], it[0], it[1]))

def WriteNucleotideHash(nucleotide_hash, out_file):
	out_fh = open(out_file, 'w')
	for it in nucleotide_hash:
		out_fh.write("%s	%s\n" %(it, nucleotide_hash[it]))
