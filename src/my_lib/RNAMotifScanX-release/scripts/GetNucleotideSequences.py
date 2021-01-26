#!/usr/bin.python2.6
#!python

import argparse
import os
import sys
import re

parser = argparse.ArgumentParser()

parser.add_argument("input", type = str, help = "The input multi-fasta file (should named pdb_seqres.txt)")
parser.add_argument("output", type = str, help = "The output multi-fasta file containing only nucleotide sequences")

args = parser.parse_args()

in_file = args.input
out_file = args.output

in_fh = open(in_file, 'r')
out_fh = open(out_file, 'w')

record_tag = False

for line in in_fh:
	if re.search('^>.*mol:na.*', line):
		record_tag = True
		out_fh.write(line)
	elif record_tag == True:
		record_tag = False
		out_fh.write(line)
