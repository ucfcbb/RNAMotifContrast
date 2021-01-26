#!/usr/bin/env python

import os, glob, sys, operator, math, re

if len(sys.argv) < 3:
	print "Error in executing the script: Argument missing!"
	print "Usage: Reconsolidate_Results.py result_file target_DIR"
	
result_file_name = os.path.abspath(sys.argv[1])
target_dir = os.path.abspath(sys.argv[2])

RIN = open(result_file_name, 'r');
content = RIN.readlines();
RIN.close();

hits_info = []

x_bar = 0
y_bar = 0
xsquare_bar = 0
xy_bar = 0

if(len(content) < 2):
	for i in range(len(content)):
		print "%s" % (content[i].rstrip())
	exit(0)

#	reads in the information in the results file and sort it

errorpattern = re.compile('ERROR')

for i in range(len(content)):
	if (not errorpattern.search(content[i])):
		#print "%s" % (content[i])
		line_info = content[i].rstrip().split("\t");
		hits_info.append(line_info)
		pindex = len(hits_info) - 1
		hits_info[pindex][2] = float(hits_info[pindex][2])	
		hits_info[pindex][3] = float(hits_info[pindex][3])
		x_bar += math.sqrt(1 / hits_info[pindex][3])
		y_bar += hits_info[pindex][2]
		xsquare_bar += math.sqrt(1 / hits_info[pindex][3]) * math.sqrt(1 / hits_info[pindex][3])
		xy_bar += hits_info[pindex][2] * math.sqrt(1 / hits_info[pindex][3])
	
	
hits_info.sort(key=operator.itemgetter(2), reverse=True);	
	
#	re-estimate the mean and variance
x_bar = x_bar / len(content)
y_bar = y_bar / len(content)
xsquare_bar = xsquare_bar / len(content)
xy_bar = xy_bar / len(content)
	
reestimated_std = (xy_bar - x_bar * y_bar) / (xsquare_bar - x_bar * x_bar)
reestimated_mean = y_bar - reestimated_std * x_bar

for i in range(len(hits_info)):
	hits_info[i][3] = (reestimated_std / (hits_info[i][2] - reestimated_mean))**2
	

#	detect overlapping of the hits
for i in range(len(hits_info)):
	HIN = open(os.path.join(target_dir, "%s.rmf" % (hits_info[i][1])));
	mf_contents = HIN.readlines()
	HIN.close()
	starts = hits_info[i][1].split('-')
	starts = starts[0].split('_')
	seq_lens = mf_contents[1].rstrip().split('...')
	#	add label to indicate the hits in not traversed or retained
	hits_info[i].append('Y')
	for j in range(1, len(starts)):
		hits_info[i].append(int(starts[j]))
		hits_info[i].append(int(starts[j]) + len(seq_lens[j - 1]) - 1);

def len_overlap(a_start, a_end, b_start, b_end):
	if b_end >= a_start and b_end <= a_end:
		if b_start < a_start:
			return b_end - a_start + 1
		else:
			return b_end - b_start + 1
	elif a_end >= b_start and a_end <= b_end:
		if a_start < b_start:
			return a_end - b_start + 1
		else:
			return a_end - a_start + 1
	else:
		return 0
		
for i in range(len(hits_info) - 1):
	for j in range(i + 1, len(hits_info)):
		if hits_info[i][4] == 'N' or hits_info[j][4] == 'N':
			continue
		significant_overlap = False
		for k in range(5, len(hits_info[i]), 2):
			for l in range(5, len(hits_info[j]), 2):
				if len_overlap(hits_info[i][k], hits_info[i][k + 1], hits_info[j][l], hits_info[j][l + 1]) > 2:
					significant_overlap = True
		if significant_overlap:
			hits_info[j][4] = 'N'

				
for i in range(len(hits_info)):
	if hits_info[i][4] == 'Y':
		print "%s	%s	%s	%s" % (hits_info[i][0], hits_info[i][1], hits_info[i][2], hits_info[i][3])
