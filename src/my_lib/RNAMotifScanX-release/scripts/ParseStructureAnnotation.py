# module for parsing MC-Annotate and RNAVIEW inputs

import re

def GetMCAnnotateInteractions(mc_annotate_file):
	mca_fh = open(mc_annotate_file, 'r')
	contents = []
	start_recording = False	
	for line in mca_fh:
		if re.search('^Adjacent\sstackings', line):
			start_recording = True
			contents.append(line)
		if start_recording:
			contents.append(line)
	# parsing the information
	interactions = []
	non_adjacent_stacking_info = False
	adjacent_stacking_info = False
	pairing_info = False	
	for info in contents:
		if re.search('Non-Adjacent stackings', info):
			non_adjacent_stacking_info = True
			adjacent_stacking_info = False
			pairing_info = False
		elif re.search('Adjacent stackings', info):
			non_adjacent_stacking_info = False
			adjacent_stacking_info = True
			pairing_info = False
		elif re.search('Base\-pairs', info):
			non_adjacent_stacking_info = False
			adjacent_stacking_info = False
			pairing_info = True
		# parsing interactions
		if re.search('\:', info):
			decom = re.split('\s+', info)
			residues = re.split('\-', decom[0])
			if non_adjacent_stacking_info:
				if decom[2] == 'inward' or decom[2] == 'outward' or decom[2] == 'upward' or decom[2] == 'downward':	
					interactions.append([residues[0], residues[1], decom[2]])
			if adjacent_stacking_info:
			  if decom[3] == 'inward' or decom[3] == 'outward' or decom[3] == 'upward' or decom[3] == 'downward':	
					interactions.append([residues[0], residues[1], decom[3]])
			if pairing_info:
				m = re.search('([HWS]).*\/([HWS]).*', decom[3])
				if m:
					edge = m.group(1) + '/' + m.group(2)
					if re.search('cis', info):
						interactions.append([residues[0], residues[1], edge, 'cis'])
					elif re.search('trans', info):	
						interactions.append([residues[0], residues[1], edge, 'trans'])
				#elif re.search('O2', decom[3]):
				#  interactions.append([residues[0], residues[1], 'S/S', 'trans'])
				else:
					interactions.append([residues[0], residues[1], '-/-', 'hbond'])				
	return interactions

def GetRNAVIEWInteractions(rnaview_file):
	rvw_fh = open(rnaview_file, 'r');
	start_parsing = False
	interactions = []
	for line in rvw_fh:
		line = ' ' + line
		if re.search('BEGIN\_base\-pair', line):
			start_parsing = True
		elif re.search('END\_base\-pair', line):
			start_parsing = False
			break
		# annotating base interactions if indicated to start parsing
		if start_parsing:
			decom = re.split('\s+', line)
			if (len(decom) >= 9 and re.search('[HWShws+-]\/[HWShws+-]', decom[7]) and (decom[8] == 'cis' or decom[8] == 'tran')):
				decom[2] = decom[2].rstrip(':')
				decom[6] = decom[6].rstrip(':')
				if not decom[2] == decom[6]:
					continue
				if re.search('\d', decom[2]):
					decom[2] = '\'' + decom[2] + '\''
				if re.search('\d', decom[6]):
					decom[6] = '\'' + decom[6] + '\''
				decom[2] = decom[2] + decom[3]
				decom[6] = decom[6] + decom[5]
				if decom[7] == '+/+' or decom[7] == '-/-':
					decom[7] = 'W/W'				
				decom[7] = decom[7].upper()
				if decom[8] == 'tran':
					decom[8] = 'trans'
				#if (re.search('!', line)):
				#	decom[7] = '-/-'
				#	decom[8] = 'hbond'
				interactions.append([decom[2], decom[6], decom[7], decom[8]])
	return interactions


def MergeInteractions(mca_interactions, rvw_interactions):
	merged_interactions = []
	interaction_hash = {}
	for single_interaction in mca_interactions:
		merged_interactions.append(single_interaction)		
		key = single_interaction[0] + '_' + single_interaction[1]
		interaction_hash[key] = 1
	for single_interaction in rvw_interactions:
		key = single_interaction[0] + '_' + single_interaction[1]
		if not key in interaction_hash:
			merged_interactions.append(single_interaction)
	return merged_interactions
