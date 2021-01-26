import os
import sys
import logging
import pickle
sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from my_log import *
from utils import *

def get_node_list_from_rmsd_data(rmsd_data_dict):
    node_list_dict = {}

    for c_id in rmsd_data_dict:

        _, rmsd_data_list_dict = rmsd_data_dict[c_id]
        key_list = list(rmsd_data_list_dict.keys())

        if len(key_list) > 0:
            (i, r1) = key_list[0]
            node1 = strToNode(r1)
            _, fit_ret = rmsd_data_list_dict[(i, r1)]
            node2_list = list(map(lambda x: strToNode(x[1]), fit_ret))
            node_list_dict[c_id] = [node1] + node2_list

    return node_list_dict

# def get_formatted_rmsd_data(rmsd_data_dict):
#     rmsd_formatted_data = {}
#     for c_id in rmsd_data_dict:
#         rmsd_formatted_data[c_id] = {}
#         _, rmsd_data_list_dict = rmsd_data_dict[c_id]
#         for i, r1 in rmsd_data_list_dict:
#             node1 = strToNode(r1)
#             if node1 not in rmsd_formatted_data[c_id]:
#                 rmsd_formatted_data[c_id][node1] = []
#             _, fit_ret = rmsd_data_list_dict[(i, r1)]
#             for _, r2, _, _ in fit_ret:
#                 node2 = strToNode(r2)
#                 if node2 not in rmsd_formatted_data[c_id][node1]:
#                     rmsd_formatted_data[c_id][node1].append(node2)
#     return rmsd_formatted_data

def is_valid_graph(clusters, graph_fname):
    if not os.path.isfile(graph_fname):
        return False

    fp = open(graph_fname)
    lines = fp.readlines()
    fp.close()

    graph_data = {}
    for line in lines:
        r1, r2, _, _, _ = line.strip().split(' ')
        if r1 not in graph_data:
            graph_data[r1] = []
        if r2 not in graph_data:
            graph_data[r2] = []

        graph_data[r1].append(r2)
        graph_data[r2].append(r1)

    if align_all_pair == False:
        for c_id in clusters:
            r1 = clusters[c_id][0]
            if r1 not in graph_data:
                return False

            r1_graph_data = [r1]
            r1_graph_data += graph_data[r1]
            set_a = set(clusters[c_id])
            set_b = set(r1_graph_data)
            if not(len(set_a) == len(set_b) and len(set_a.intersection(set_b)) == len(set_a)):
                return False
        return True
    else:
        all_loops = []
        for c_id in clusters:
            all_loops += clusters[c_id]

        r1 = all_loops[0]
        if r1 not in graph_data:
            return False

        r1_graph_data = [r1]
        r1_graph_data += graph_data[r1]
        set_a = set(all_loops)
        set_b = set(r1_graph_data)
        if len(set_a) == len(set_b) and len(set_a.intersection(set_b)) == len(set_a):
            return True

    return False

def is_valid_pickle(pickle_fname, clusters):

    if os.path.basename(pickle_fname).startswith('alignment'):
        # check alignment pickle
        alignment_data_fname = pickle_fname

        if not os.path.isfile(alignment_data_fname):
            return False

        f = open(alignment_data_fname, 'rb')
        cluster_alignment_data = pickle.load(f)
        f.close()
        
        for c_id in clusters:
            if c_id not in cluster_alignment_data:
                return False

        for c_id in clusters:
            loop_nodes = list(map(lambda x: strToNode(x), clusters[c_id]))
            node1 = loop_nodes[0]
            set_a = set(loop_nodes)
            set_b = set(list(cluster_alignment_data[c_id][node1].keys()) + [node1])

            if not(len(set_a) == len(set_b) and len(set_a.intersection(set_b)) == len(set_a)):
                return False
            
    elif os.path.basename(pickle_fname).startswith('rmsd'):
        # hard to detect if this is invalid or not
        # check rmsd pickle
        rmsd_data_fname = pickle_fname

        if not os.path.isfile(rmsd_data_fname):
            return False

        f = open(rmsd_data_fname, 'rb')
        rmsd_data_dict = pickle.load(f)
        f.close()

        for c_id in clusters:
            if c_id not in rmsd_data_dict:
                return False

        node_list_dict = get_node_list_from_rmsd_data(rmsd_data_dict)

        for c_id in clusters:
            loop_nodes = list(map(lambda x: strToNode(x), clusters[c_id]))
            set_a = set(loop_nodes)
            set_b = set(node_list_dict[c_id])

            if not(len(set_a) == len(set_b) and len(set_a.intersection(set_b)) == len(set_a)):
                return False

    else:
        return False

    return True

def validate_all(input_fname, draw_figures):
    if not os.path.isfile(input_fname):
        logger.error('Input file does not exists.')
        sys.exit()

    if draw_figures == True:
        if sys.version_info > (3, 0):
            # if len(pymol_py3_dir) == 0:
            if not os.path.exists(pymol_py3_dir):
                logger.error('Please configure pymol, and set pymol directory in ' + os.path.join(root_dir, 'config.py')[base_path_len:] + ' file to generate images. (see instructions in README file)')
                sys.exit()
        else:
            try:
                import pymol
            except Exception as e:
                logger.error('Please install pymol to generate images. (see instructions in README file)')
                sys.exit()

	# if draw_figures == True and (sys.version_info > (3, 0)) and len(pymol_py3_dir) == 0:
	# 	logger.error('Please set value for \'pymol_py3_dir\' in \'config.py\' file.')
	# 	sys.exit()
