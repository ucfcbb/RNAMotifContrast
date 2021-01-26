#!/usr/bin/python
#!python

# The script is used to generate input file for RNAMotifScanX
# The script accepts a PDB file and use MC-Annotate and RNAVIEW to parse
#	the file into base-interaction patterns 
# by Cuncong Zhong (send comments or bug reports to czhong@jvci.org)

import MapPDBSequence
import LoadSequence
import ParseStructureAnnotation
import FormattedOutput

import argparse
import os
import sys
import re


parser = argparse.ArgumentParser()

# the input PDB structure
parser.add_argument("pdb_file", type = str, help = "The PDB file to be annotated")
# the reference sequence file
parser.add_argument("seq_ref", type = str, help = "The reference sequence file (in (multi)FASTA format)")

path = os.environ["RNAMOTIFSCANX_PATH"]
mc_exe = path + "/thirdparty/MC-Annotate"
rv_exe = path + "/thirdparty/RNAVIEW/bin/rnaview"
# the MC-Annotate executable
parser.add_argument("-m", "--mc_annotate_exe", type = str, help = "The MC-Annotate executable location", default = mc_exe)
# the RNAVIEW exetutable
parser.add_argument("-r", "--rnaview_exe", type = str, help = "The RNAVIEW executable location", default=rv_exe)
# whether incorporate RNAVIEW annotation
parser.add_argument("--incorporate_rnaview", action = "store_true", help = "Incoporate RNAVIEW annotation", default = False)
# whether to enforce the entire process disregarding whether annotation files already exist
parser.add_argument("-f", "--force", action = "store_true", help = "Force to redo the input preparation and overwrite existing files", default = False)

# specify use which annotation
#group = parser.add_mutually_exclusive_group()
#group.add_argument("--use_mc_annotate", action = "store_true", help = "Use MC-Annotate annotation only", default = False)
#group.add_argument("--use_rnaview", action = "store_true", help = "Use RNAVIEW annotation only", default = False)
#group.add_argument("--use_both", action = "store_true", help = "Use (the union of) both annotations", default = True)

# verbosity
parser.add_argument("-v", "--verbose", action = "store_true", help = "Enable progress report", default = True)

# passing in arguments
args = parser.parse_args()
structure_to_parse = args.pdb_file
reference_sequence_file = args.seq_ref
mc_annotate_executable = args.mc_annotate_exe
rnaview_executable = args.rnaview_exe
annotation_both = True
verbose_enabled = args.verbose
force_enabled = args.force

# Guessing the PDB code and the directory	
structure_prefix = os.path.realpath(structure_to_parse)
path_search_obj = re.search('.*\/', structure_prefix)
pdb_path = path_search_obj.group(0);
pdb_code = structure_prefix[len(pdb_path) : len(pdb_path) + 4]


if not os.path.exists(mc_annotate_executable):
	sys.exit("MC-Annotate executable does not exist!")
if annotation_both and not os.path.exists(rnaview_executable):
	sys.exit("RNAVIEW executable does not exist!")

if not re.search('[0-9A-Z]{4,4}', pdb_code):
	print "The name of valid PDB files should have the following format: XXXX.pdb ('X' is either a single digit or a CAPTALIZED letter)."
	print "If you are confortable that the input file is correct, please rename it and try again."
	sys.exit()
if not os.path.exists(structure_to_parse):
	print "The PDB file does not exist in local archive. Trying to download it from PDB."
	os.system("wget http://www.rcsb.org/pdb/files/%s.pdb.gz -o %s.wget.log -O %s.gz" %(pdb_code, structure_to_parse, structure_to_parse))
	os.system("gzip -d %s.gz" %structure_to_parse)
	if not os.path.exists(structure_to_parse):
		print "No such structure can be downloaded from PDB. Please double check you input."
		sys.exit()


# running MC-Annotate and RNAVIEW
result_file_mca = ""
result_file_rvw = ""
result_file_rvw_notNMR = ""
result_file_rvw_NMR = ""

result_file_mca = structure_prefix + ".mca"	
if not os.path.exists(result_file_mca) or force_enabled:	
	if verbose_enabled:
		print "Running MC-Annotate, please wait..."
	os.system("%s %s >%s.mca" %(mc_annotate_executable, structure_to_parse, structure_to_parse))
else:	
	if verbose_enabled:
		print "MC-Annotate file %s already exists, skip annotation." %result_file_mca

if annotation_both:
	result_file_rvw_notNMR = structure_prefix + ".out"
	result_file_rvw_NMR = structure_prefix + "_nmr.pdb.out"	
	if (not os.path.exists(result_file_rvw_notNMR) and not os.path.exists(result_file_rvw_NMR)) or force_enabled:
		if verbose_enabled:
			print "Running RNAVIEW, please wait..."			
		os.system("%s %s" %(rnaview_executable, structure_to_parse))	
	elif os.path.exists(result_file_rvw_notNMR):
		if verbose_enabled:
			print "RNAVIEW annotation file %s already exists, skip annotation." %result_file_rvw_notNMR
	elif os.path.exists(result_file_rvw_NMR):			
		if verbose_enabled:
			print "RNAVIEW annotation file %s already exists, skip annotation." %result_file_rvw_NMR

	if os.path.exists(result_file_rvw_notNMR):
		result_file_rvw = result_file_rvw_notNMR		
	elif os.path.exists(result_file_rvw_NMR):			
		result_file_rvw = result_file_rvw_NMR


# Parse MC-Annotate and RNAVIEW files
mca_interactions = ParseStructureAnnotation.GetMCAnnotateInteractions(result_file_mca)
rvw_interactions = []
if annotation_both:
  rvw_interactions = ParseStructureAnnotation.GetRNAVIEWInteractions(result_file_rvw)
merged_interactions = ParseStructureAnnotation.MergeInteractions(mca_interactions, rvw_interactions)


# Producing the mapping between sequence and PDBids based on sequence alignment
nucleotides_in_structure = MapPDBSequence.GetRecognizedNucleotides(result_file_mca);
chains_in_structure = MapPDBSequence.GetRecognizedChains(nucleotides_in_structure);

for chain in chains_in_structure:
	seq_tag = pdb_code + '_' + chain
	# loads the sequence from the reference sequence set
	if verbose_enabled:
		print "Building nucleotide hash for %s..." %seq_tag
	interaction_out_file = pdb_path + seq_tag + '.rmsx.in'
	nuchash_out_file = pdb_path + seq_tag + '.rmsx.nch'
	if not os.path.exists(interaction_out_file) or not os.path.exists(nuchash_out_file) or force_enabled:
		chain_ref_seq_info = LoadSequence.LoadSequenceWithTag(reference_sequence_file, seq_tag);
		chain_mapping = MapPDBSequence.AlignToReference(chain, chain_ref_seq_info[1], nucleotides_in_structure);
		if verbose_enabled:
			print "Formatted annotation and nucleotide mappings have been written to %s and %s." %(interaction_out_file, nuchash_out_file)
		FormattedOutput.WriteInteractions(pdb_code, chain, chain_ref_seq_info[1], merged_interactions, chain_mapping, interaction_out_file)
		FormattedOutput.WriteNucleotideHash(chain_mapping, nuchash_out_file)
	else:
		if verbose_enabled:
			print "RNAMotifScanX input files %s and %s already exist, skip parsing annotations." %(interaction_out_file, nuchash_out_file)
	

