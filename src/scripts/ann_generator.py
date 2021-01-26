import sys
import os
import logging

from future.standard_library import install_aliases
install_aliases()
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request, urlretrieve
from urllib.error import HTTPError

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from ann_parser import *

def generate_dssr_annotation_for_single_pdb(selected_annotation_dir, pdb_id):

	if output_env == 'global':
		if annotation_source == 'dssr':
			logger.error('DSSR annotation file not exists for ' + pdb_id + '.')
			sys.exit()
		else:
			# logger.warning('DSSR annotation file not exists for ' + pdb_id + '.')
			return

	pdb_fname = os.path.join(pdbx_dir, pdb_id + '.cif')
	t_dir = os.path.join(temp_dir, pdb_id)
	create_directory(t_dir)
	os.chdir(t_dir)
	# os.system("./x3dna-dssr -i=%s -o=%s --non-pair" % (pdb_fname, os.path.join(selected_annotation_dir, pdb_id + '.dssr')))
	os.system("%s -i=%s -o=%s --non-pair" % (os.path.join(dssr_dir, "x3dna-dssr"), pdb_fname, os.path.join(selected_annotation_dir, pdb_id + '.dssr')))

def download_fr3d_annotation_for_single_pdb(selected_annotation_dir, pdb_id):
	try:
		csvfile = urlopen(fr3d_url.replace('XXXX', pdb_id))
		fp = open(os.path.join(selected_annotation_dir, pdb_id+'.fr3d'), 'wb')
		fp.write(csvfile.read())
		fp.close()
		# logger.info('File downloaded successfully: ' + os.path.join(selected_annotation_dir, pdb_id+'.fr3d'))
	except HTTPError as e:
		if annotation_source == 'fr3d':
			logger.error('Error downloading fr3d annotation file for ' + pdb_id + '. ' + e.Read())
			sys.exit()
		else:
			# logger.warning('Error downloading fr3d annotation file for ' + pdb_id + '. ' + e.Read())
			return

def get_highest_freq_interactions(interactions):
	interact_dict = {}
	for interact in interactions:
		if interact not in interact_dict:
			interact_dict[interact] = 0
		interact_dict[interact] += 1

	max_freq = 0
	for interact in interact_dict:
		max_freq = max(interact_dict[interact], max_freq)

	max_freq_interactions = []
	for interact in interact_dict:
		if interact_dict[interact] == max_freq:
			max_freq_interactions.append(interact)

	return max_freq_interactions

def get_resolved_dict(pdb_id, pairwise_dict, interaction_category_rank):
	resolved_dict = {}
	for ind_pair in pairwise_dict:
		interactions = get_highest_freq_interactions(pairwise_dict[ind_pair])
		if len(interactions) > 1:
			best_rank = 1000
			best_rank_interaction = ("", "")
			for bp, interact in interactions:
				if bp not in interaction_category_rank:
					logger.error('ERROR: Problem in ' + pdb_id + '. Check bp ' + bp + '.')
					sys.exit()
				if interact not in interaction_category_rank[bp]:
					rank = len(interaction_category_rank[bp])
				else:
					rank = interaction_category_rank[bp].index(interact)

				if rank < best_rank:
					best_rank = rank
					best_rank_interaction = (bp, interact)
			interactions = []
			interactions.append((best_rank_interaction[0], best_rank_interaction[1]))

		resolved_dict[ind_pair] = interactions

	return resolved_dict

def write_merged_annotation_to_file(pdb_id, ann_list, merged_annotation_dir, detailed):
	if detailed:
		bp_item_len = 5
		stk_item_len = 4
	else:
		logger.error('Detailed BP info not available.')
		sys.exit()

	# for pdb_id in ann_dict:
	fp = open(os.path.join(merged_annotation_dir, pdb_id + ".merged"), "w")
	for annotation in ann_list:
		if len(annotation) == bp_item_len:
			index1, index2, edges, orientation, bp = annotation
			fp.write(str(index1) + '\t' + str(index2) + '\t' + bp + '\t' + edges + '\t' + orientation + '\n')

		elif len(annotation) == stk_item_len:
			index1, index2, direction, bp = annotation
			fp.write(str(index1) + '\t' + str(index2) + '\t' + bp + '\t' + direction + '\n')

		else:
			logger.error('ERROR: invalid interact_info length in ' + pdb_id)
			sys.exit()
	fp.close()

def resolve_pairwise_annotation_conflict_helper(pdb_id, ann_list, bp_interaction_category_rank, stk_interaction_category_rank, merged_annotation_dir, detailed):
	if detailed:
		bp_item_len = 5
		stk_item_len = 4
	else:
		logger.error('Detailed BP info not available.')
		sys.exit()

	pairwise_bp_dict = {}
	pairwise_stk_dict = {}

	logger.info('Resolving ' + pdb_id)
	for interact_info in ann_list:
		if len(interact_info) == bp_item_len:
			index1, index2, edges, orientation, bp = interact_info
			interaction = orientation[0] + edges[0] + edges[2]
			if (index1, index2) not in pairwise_bp_dict:
				pairwise_bp_dict[(index1, index2)] = []
			pairwise_bp_dict[(index1, index2)].append((bp, interaction))

		elif len(interact_info) == stk_item_len:
			index1, index2, direction, bp = interact_info
			if (index1, index2) not in pairwise_stk_dict:
				pairwise_stk_dict[(index1, index2)] = []
			pairwise_stk_dict[(index1, index2)].append((bp, direction))

		else:
			logger.error('ERROR: invalid interact_info length in ' + pdb_id)
			sys.exit()

	pairwise_bp_dict = get_resolved_dict(pdb_id, pairwise_bp_dict, bp_interaction_category_rank)
	pairwise_stk_dict = get_resolved_dict(pdb_id, pairwise_stk_dict, stk_interaction_category_rank)

	resolved_annotation_list = []
	for (index1, index2) in pairwise_bp_dict:
		(bp, interaction) = pairwise_bp_dict[(index1, index2)][0]
		orientation = 'cis' if interaction[0] == 'c' else 'trans'
		edges = interaction[1] + '/' + interaction[2]
		resolved_annotation_list.append((index1, index2, edges, orientation, bp))

	for index1, index2 in pairwise_stk_dict:
		bp, direction = pairwise_stk_dict[(index1, index2)][0]
		resolved_annotation_list.append((index1, index2, direction, bp))

	resolved_annotation_list = sorted(resolved_annotation_list)
	write_merged_annotation_to_file(pdb_id, resolved_annotation_list, merged_annotation_dir, detailed)

def load_category_rank_from_file(fname):
	category_rank = {}

	fp = open(fname)
	rank_list = csv_to_list(fp.readlines())
	fp.close()

	for item in rank_list:
		category_rank[item[0]] = item[1:]

	return category_rank

def get_dssr_and_fr3d_merged_annotation_for_single_pdb(selected_annotation_dir, pdb_id):
	dssr_ann_dir = os.path.join(annotation_dir, 'dssr')
	fr3d_ann_dir = os.path.join(annotation_dir, 'fr3d')
	dssr_fname = os.path.join(dssr_ann_dir, pdb_id + '.dssr')
	fr3d_fname = os.path.join(fr3d_ann_dir, pdb_id + '.fr3d')
	if not os.path.isfile(dssr_fname):
		generate_dssr_annotation_for_single_pdb(dssr_ann_dir, pdb_id)
	if not os.path.isfile(fr3d_fname):
		download_fr3d_annotation_for_single_pdb(fr3d_ann_dir, pdb_id)

	detailed_info = True
	ann_list = parseDSSR(dssr_fname, detailed_info)
	
	if os.path.isfile(fr3d_fname):
		ann_list.extend(parseFR3D(fr3d_fname, detailed_info))

	if len(ann_list) == 0:
		logger.error('None of the DSSR and FR3D annotation files were found. Cannot generate the merged annotation file for ' + pdb_id + '.')
		sys.exit()

	bp_interaction_category_rank = load_category_rank_from_file(os.path.join(lib_dir, 'bp_category_rank.csv'))
	stk_interaction_category_rank = load_category_rank_from_file(os.path.join(lib_dir, 'stk_category_rank.csv'))

	resolve_pairwise_annotation_conflict_helper(pdb_id, ann_list, bp_interaction_category_rank, stk_interaction_category_rank, selected_annotation_dir, detailed_info)
