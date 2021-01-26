import sys
import os
import glob
import logging
import operator
import time
import pickle
import multiprocessing as mp
import numpy as np
from datetime import datetime
from Bio.PDB import *
import copy
import gc

# python 3 compatibility
from functools import reduce
from past.builtins import map

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from utils import *
from pymol_helper import *
from classes import *
from validators import *

def rmsd(V, W):
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def kabsch(P, Q):
    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)
    return U

def kabsch_rmsd(P, Q):
    P = rotate(P, Q)
    return rmsd(P, Q)


def rotate(P, Q):
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P

def convert_array(c1, c2):
    ret1 = []
    ret2 = []
    for i, j in zip(c1, c2):
        if i == 0 or j == 0:
            continue
        ret1.append(i)
        ret2.append(j)
    return np.array(ret1), np.array(ret2)

def load_graph_data(fn):
    graph_data = {}
    fp = open(fn)
    cnt = 0
    for line in fp.readlines():
        pieces = line.split()

        node1 = strToNode(pieces[0])
        node2 = strToNode(pieces[1])

        if node1 not in graph_data:
            graph_data[node1] = {}
        if node2 not in graph_data:
            graph_data[node2] = {}

        # t1 t2 is for best alignment order
        #(t1, t2, z-score, local_aln_indx, local_aln_indx, seq1, seq2, aln_score)
        if node1 == strToNode(pieces[3]):
            graph_data[node1][node2] = (pieces[3], pieces[4], float(pieces[2]), '', '', '', '', 0.)
            graph_data[node2][node1] = (pieces[4], pieces[3], float(pieces[2]), '', '', '', '', 0.)
        else:
            graph_data[node1][node2] = (pieces[4], pieces[3], float(pieces[2]), '', '', '', '', 0.)
            graph_data[node2][node1] = (pieces[3], pieces[4], float(pieces[2]), '', '', '', '', 0.)

        cnt += 1

    fp.close()

    # print('Maximum number of edges in cluster graph: ' + str(cnt))

    return graph_data

def load_cluster_alignment_data(alignment_data, clusters):
    cluster_alignment_data = {}

    for c_id in clusters:
        cluster_alignment_data[c_id] = {}

        for i in range(len(clusters[c_id])):

            node1 = strToNode(clusters[c_id][i])

            if node1 not in cluster_alignment_data[c_id]:
                cluster_alignment_data[c_id][node1] = {}

            for j in range(len(clusters[c_id])):
                if i == j:
                    continue

                node2 = strToNode(clusters[c_id][j])

                if node2 not in cluster_alignment_data[c_id]:
                    cluster_alignment_data[c_id][node2] = {}

                if node1 not in alignment_data or node2 not in alignment_data[node1]:
                    logger.error('ERROR: (' + str(node1) + ', ' + str(node2) + ') pair not found in graph file!')
                    sys.exit()

                cluster_alignment_data[c_id][node1][node2] = alignment_data[node1][node2]
                cluster_alignment_data[c_id][node2][node1] = alignment_data[node2][node1]

    return cluster_alignment_data, get_loops_in_cluster(clusters)

def load_alignment_data(input_fname_base, alignment_dir, graphs_and_pickles_dir, alignment_fname, clusters, previous_graph_file_reused):

    alignment_data_fname = os.path.join(graphs_and_pickles_dir, 'alignment_data_' + input_fname_base + '.pickle2')
    if (sys.version_info >= (3, 0)):
        alignment_data_fname = os.path.join(graphs_and_pickles_dir, 'alignment_data_' + input_fname_base + '.pickle3')

    if previous_graph_file_reused == True:
        if is_valid_pickle(alignment_data_fname, clusters) == True:
            logger.info('Loading saved alignment data from previous run (In ' + alignment_data_fname[base_path_len:] + ') ...\n')
            f = open(alignment_data_fname, 'rb')
            cluster_alignment_data = pickle.load(f)
            f.close()
            return cluster_alignment_data

    logger.info('Loading alignment data from alignment files.')
    start_time = time.time()

    alignment_data = load_graph_data(alignment_fname)
    cluster_alignment_data, loops_node_list = load_cluster_alignment_data(alignment_data, clusters)

    # print('Number of loops in cluster file: ' + str(len(loops_in_cluster)))

    file_counter = 0
    for node in loops_node_list:
        for r1 in get_all_loop_combination(str(node)):
            node1 = strToNode(r1)
    # for fn in glob.glob(os.path.join(alignment_dir, '*.aln')):
            fn = os.path.join(alignment_dir, r1 + '.aln')

            if not os.path.isfile(fn):
                return None

            # stime = datetime.now()
            file_counter += 1
            print('Processsing ' + fn[base_path_len:] + ' ... (' + str(file_counter) + ')')
            # sys.stdout.flush()
            # r1 = os.path.basename(fn)[:-4]
            # node1 = strToNode(r1)

            # if node1 not in loops_in_cluster:
            #     continue

            cid_nodelist_pair = find_nodes_in_cluster(node1, cluster_alignment_data)

            fp = open(fn)
            lines = fp.readlines()
            fp.close()
            # flag = 0
            # print('No. of lines: ' + str(len(lines)))
            # print('Connected nodes: ' + str(len(node_dict)))
            line_index = 0
            while line_index < len(lines):
                # print 'Reading line ' + str(line_index)
                # sys.exit()
                if lines[line_index].startswith('#  Aligning'):
                    test_r1 = lines[line_index].split('::')[1].split(' and ')[0].strip().strip(':')
                    if r1 != test_r1:
                        logger.error('ERROR: filename and loop mismatch. r1(filename): ' + r1 + ', r1(loop): ' + test_r1)
                        sys.exit()

                    r2 = lines[line_index].split('::')[1].split(' and ')[1].strip().strip(':')
                    node2 = strToNode(r2)

                    for c_id, node_dict in cid_nodelist_pair:
                        if node2 in node_dict:

                            t1 = cluster_alignment_data[c_id][node1][node2][0]
                            t2 = cluster_alignment_data[c_id][node1][node2][1]
                            # make sure the best alignment ordering is same as the ordering of loop in current file
                            if not ((r1 == t1 and r2 == t2) or (r1 == t2 and r2 == t1)):
                                # line_index += 12
                                continue

                            #safety condition for reverse order (might be redundant)
                            if (r1 == t2 and r2 == t1):
                                t1 = r2
                                t2 = r1

                            score_text = lines[line_index+1].split(':')[1].strip()
                            if score_text == '':
                                score = -50.
                                logger.error('ERROR: No alignment score found for: ' + r1 + ' and ' + r2)
                                sys.exit()
                            else:
                                score = float(score_text)

                            text = lines[line_index+3].split(':')[1].strip()
                            cr1 = get_local_alignment_index(r1, text)

                            text = lines[line_index+4].split(':')[1].strip()
                            cr2 = get_local_alignment_index(r2, text)

                            aln1 = lines[line_index+6].strip()
                            aln2 = lines[line_index+7].strip()

                            if len(aln1) == 0 or len(aln2) == 0:
                                # set dummy alignment
                                aln1 = 'A'
                                aln2 = 'A'
                                temp_cr1 = cr1.split(':')
                                dummy_index = temp_cr1[1].split('_')[0].split('-')[0]
                                cr1 = temp_cr1[0] + ':' + dummy_index + '-' + dummy_index
                                temp_cr2 = cr2.split(':')
                                dummy_index = temp_cr2[1].split('_')[0].split('-')[0]
                                cr2 = temp_cr2[0] + ':' + dummy_index + '-' + dummy_index

                            zscore = cluster_alignment_data[c_id][node1][node2][2]
                            cluster_alignment_data[c_id][node1][node2] = (t1, t2, zscore, cr1, cr2, aln1, aln2, score)
                            cluster_alignment_data[c_id][node2][node1] = (t2, t1, zscore, cr2, cr1, aln2, aln1, score)
                            # break
                    # we can skip at least 11 lines from current alignment data
                    line_index += 11
                line_index += 1
            # while end
            
            # etime = datetime.now()
            # difftime = (etime - stime).total_seconds()
            # print('Time taken: ')
            # print(difftime)
            # break

    f = open(alignment_data_fname,"wb")
    pickle.dump(cluster_alignment_data, f)
    f.close()

    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

    return cluster_alignment_data

def extract_atom_coordinate(loop_coord, pdb_pm, pdb_id):
    loop_coord_backbone, loop_coord_sugar = loop_coord

    aligned_segment_coord = []
    for chain_id, index, icode in pdb_pm:
        if chain_id == '':
            continue
        if (pdb_id, chain_id, index, icode) in loop_coord_backbone and loop_coord_backbone[(pdb_id, chain_id, index, icode)] != 0.:
            aligned_segment_coord.append(loop_coord_backbone[(pdb_id, chain_id, index, icode)])
        if (pdb_id, chain_id, index, icode) in loop_coord_sugar and loop_coord_sugar[(pdb_id, chain_id, index, icode)] != 0.:
            aligned_segment_coord.append(loop_coord_sugar[(pdb_id, chain_id, index, icode)])

    return aligned_segment_coord

def generate_rmsd_data(input_fname_base, partial_pdbx_dir, graphs_and_pickles_dir, alignment_data, clusters, loop_list, previous_graph_file_reused):

    rmsd_data_fname = os.path.join(graphs_and_pickles_dir, 'rmsd_data_' + input_fname_base + '.pickle2')
    if (sys.version_info >= (3, 0)):
        rmsd_data_fname = os.path.join(graphs_and_pickles_dir, 'rmsd_data_' + input_fname_base + '.pickle3')
     
    if previous_graph_file_reused == True:
        if is_valid_pickle(rmsd_data_fname, clusters) == True:
            logger.info('Loading saved RMSD data from previous run (In ' + rmsd_data_fname[base_path_len:] + ') ...\n')
            f = open(rmsd_data_fname, 'rb')
            rmsd_data_dict = pickle.load(f)
            f.close()
            return rmsd_data_dict

    # print(alignment_data['GNGA'][strToNode('4V9F_0:2621-2626')][strToNode('4V9F_0:460-465')])
    # print(alignment_data['GNAA'][strToNode('4V9F_0:2621-2626')][strToNode('4V9F_0:460-465')])
    # sys.exit()
    logger.info('Generating RMSD data.')
    start_time = time.time()

    rmsd_data_dict = {}
    coord_dict = {}
    pdb_structure = None
    prev_pdb_chain = ''
    structure_counter = 0

    # for lp in loop_list:
    #     pdb_chain, regions = lp.split(':')
    #     pdb = pdb_chain.split('_')[0]
    #     pdb_pm = get_pdb_index_list(lp)
    #     if prev_pdb_chain != pdb_chain:
    #         pdb_structure = None
    #     else:
    #         structure_counter += 1
    #         if structure_counter % 500 == 0:
    #             gc.collect()
    #     coord_backbone, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(pdbx_dir, pdb+'.cif'), pdb_pm, pdb_structure)
    #     coord_dict[lp] = (coord_backbone, coord_sugar)
    #     prev_pdb_chain = pdb_chain

    for lp in loop_list:
        pdb_pm = get_pdb_index_list(lp)
        coord_backbone, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(partial_pdbx_dir, lp + '.cif'), pdb_pm)
        coord_dict[lp] = (coord_backbone, coord_sugar)

    # time_align_residue = 0
    # time_get_coordinate = 0
    # time_rmsd = 0
    # time_start = time.time()
    for cluster_id in alignment_data:
        # if cluster_id not in clusters:
        #     continue
        sum_of_avg_rmsd_for_c = 0.
        rmsd_data_list_dict = {}
        index_dict = {}
        i = 0
        # print(loop_list)
        for l1 in alignment_data[cluster_id]:
            index_dict[l1] = i
            i += 1

        # pdb_chain_dict = {}
        # pdb_res_mapping_dict = {}
        # fasta_seq_dict = {}
        # for l1 in alignment_data[cluster_id]:
        #     pdb_chain, _ = str(l1).strip().split(':')
            
        #     pdb_id, chain_id = pdb_chain.strip().split('_')
        #     if pdb_id not in pdb_chain_dict:
        #         pdb_chain_dict[pdb_id] = []
        #     pdb_chain_dict[pdb_id].append(chain_id)

        #     if pdb_chain not in pdb_res_mapping_dict:
        #         pdb_res_mapping_dict[pdb_chain] = load_pdb_res_map(pdb_chain)

        # for pdb_id in pdb_chain_dict:
        #     fasta_seq_dict.update(load_fasta_seq(pdb_id, pdb_chain_dict[pdb_id]))

        pdb_res_mapping_dict, fasta_seq_dict = load_pdb_fasta_mapping_and_fasta_seq_dict(cluster_id, alignment_data)

        # i = 0
        for l1 in alignment_data[cluster_id]:
            # if str(l1) not in loop_list:
            #     continue
            fit_ret = []
            sum_of_rmsd_for_l1 = 0.
            # rmsd_data_list_item = {}
            # j = 0
            for l2 in alignment_data[cluster_id][l1]:
                # if str(l2) not in loop_list:
                #     continue
                if l1 != l2:
                    # print(cluster_id)
                    r1, r2, zscore, cr1, cr2, aln_1, aln_2, score = alignment_data[cluster_id][l1][l2]
                    # print('r1, r2: ')
                    # print(r1, r2)
                    # print(r1, r2, zscore, cr1, cr2, aln_1, aln_2, score)
                    # if output_env == 'local':
                    # time_s = time.time()
                    pdb1_pm, pdb2_pm, i1_pm, i2_pm = aln_residue_temp(pdb_res_mapping_dict, fasta_seq_dict, r1, r2, cr1, cr2, aln_1, aln_2, 0, len(aln_1)-1, 0)
                    # time_align_residue += time.time() - time_s
                    # else:
                    #     pdb1_pm, pdb2_pm, i1_pm, i2_pm = aln_residue(r1, r2, cr1, cr2, aln_1, aln_2, 0, len(aln_1)-1, 0)

                    pdb_chain1, _ = r1.split(':')
                    pdb_chain2, _ = r2.split(':')
                    pdb1 = pdb_chain1.split('_')[0]
                    pdb2 = pdb_chain2.split('_')[0]

                    # structures = {}
                    #returns centroid of backbone atoms
                    # time_s = time.time()
                    coord1 = extract_atom_coordinate(coord_dict[str(l1)], pdb1_pm, pdb1)
                    coord2 = extract_atom_coordinate(coord_dict[str(l2)], pdb2_pm, pdb2)
                    # time_get_coordinate += time.time() - time_s
                    # print(coord1, coord2)

                    X, Y = convert_array(coord1, coord2)

                    if len(X) != len(Y):
                        logger.warning('WARNING: Corresponding co-ordinates for alignments not found! rmsd = 20 assigned.')
                        rmsd = 20.
                    elif len(X) == 0:
                        logger.warning('WARNING: Co-ordinates for alignments not found! rmsd = 20 assigned.')
                        rmsd = 20.
                    else:
                        XC = sum(X)/len(X)
                        YC = sum(Y)/len(Y)
                        # calculating relative co-ordinate using mean as reference
                        X -= XC
                        Y -= YC
                        # time_s = time.time()
                        rmsd = kabsch_rmsd(X, Y)
                        # time_rmsd += time.time() - time_s

                    sum_of_rmsd_for_l1 += rmsd
                    # X.shape[0] represents the number of aligned nucleotides
                    # fit_ret.append((index_dict[l2], str(l2), rmsd, X.shape[0]))
                    fit_ret.append((index_dict[l2], str(l2), rmsd, len(pdb1_pm)))
                    # end of if
                # j += 1
                # end of l2 for loop

            # time_diff = time.time() - time_start
            # if time_diff > 30:
            #     print(cluster_id, l1)
            #     print('time_align_residue', time_align_residue)
            #     print('time_get_coordinate', time_get_coordinate)
            #     print('time_rmsd', time_rmsd)
            #     time_start = time.time()
            #     print('')

            avg_of_rmsd_for_l1 = 0.0
            if (len(clusters[cluster_id]) - 1) > 0:
                avg_of_rmsd_for_l1 = sum_of_rmsd_for_l1 / (len(clusters[cluster_id]) - 1)
            sum_of_avg_rmsd_for_c += avg_of_rmsd_for_l1
            # print str(i)+','+l1+'\t'+'|'.join(map(lambda x: str(x[0])+','+x[1]+','+str(x[2])+','+str(x[3]), sorted(fit_ret, key=lambda x: x[2])))
            rmsd_data_list_dict[(index_dict[l1], str(l1))] = (avg_of_rmsd_for_l1, sorted(fit_ret, key=lambda x: x[2]))
            # i += 1
            # end of l1 for loop

        avg_of_avg_rmsd_for_c = 0.0
        if len(clusters[cluster_id]) > 0:
            avg_of_avg_rmsd_for_c = sum_of_avg_rmsd_for_c / len(clusters[cluster_id])

        rmsd_data_dict[cluster_id] = (avg_of_avg_rmsd_for_c, rmsd_data_list_dict) # sorted(rmsd_data_list_dict, key=lambda x: x[0][1]) #order by loop average

    f = open(rmsd_data_fname,"wb")
    pickle.dump(rmsd_data_dict, f)
    f.close()

    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

    return rmsd_data_dict

def read_pdb_chain_organism_details(fname):
    fp = open(fname)
    first_line = True
    pdb_organism_details = {}
    for line in fp.readlines():
        if first_line:
            first_line = False
            continue
        pieces = line.strip('\n').strip('\r').split('\t')
        # For each chain, store RNA Types, Organism, Class, Type (Manually Defined), Source
        if len(pieces) > 0:
            pdb_organism_details[pieces[0].strip()] = pieces[1:]
    fp.close()
    return pdb_organism_details

def remove_based_on_zscore(zscore, rmsd, align_length, is_alignment_from_user):
    # if cluster_source == 'DeNovo':
    if is_alignment_from_user == False:
        # zscore < 0.0
        if get_zscore_rank(zscore) == 100:
            return True

        # zscore 0 to 0.5, rmsd >= 1.0
        if get_zscore_rank(zscore) >= 5 and get_rmsd_rank(rmsd, align_length, is_length_adjusted_score) >= 3:
            return True

        # zscore 0.5 to 1, rmsd >= 2.0
        if get_zscore_rank(zscore) == 4 and get_rmsd_rank(rmsd, align_length, is_length_adjusted_score) >= 4:
            return True

    return False

def remove_based_on_rmsd(rmsd, align_length):
    # rmsd >= 4.0
    if get_rmsd_rank(rmsd, align_length, is_length_adjusted_score) >= 5:
        return True
    return False

def extreme_filtering_based_on_rmsd(rmsd, align_length):
    if extreme_filtering == True and rmsd > rmsd_threshold_for_merging:
        return True
    return False

def get_loop_string_for_filtering_log(r):
    loop_str = ''
    r_pdb_ind = convert_a_loop_from_FASTA_to_PDB(r)
    loop_str = r_pdb_ind
    if input_index_type == 'fasta':
        loop_str += ' (PDB) ' + r + ' (FASTA)'

    return loop_str

def filter_loops_in_cluster(clusters, rmsd_data_dict, alignment_data, is_alignment_from_user):

    current_rmsd_data_dict = copy.deepcopy(rmsd_data_dict)
    removed_loops = 0
    while True:
        max_cid_len = max([len(x) for x in clusters])
        filtered_cluster = {}
        for cluster_id in sorted(clusters):
            loops = clusters[cluster_id]
            filtered_loops = []
            for i in range(len(loops)):
                loops[i] = str(strToNode(loops[i]))

            _, cluster_pairwise_align_details = current_rmsd_data_dict[cluster_id]
            align_len_threshold = generate_align_length_threshold(cluster_pairwise_align_details)
            for (i, r1) in cluster_pairwise_align_details:
                if r1 in loops:
                    _, pairwise_align_details = cluster_pairwise_align_details[(i, r1)]

                    j, r2, rmsd, align_length = find_best_aligned_pair(pairwise_align_details, align_len_threshold)

                    # (t1, t2, zscore, cr1, cr2, aln1, aln2, score) = alignment_data[cluster_id.strip().split('_')[0]][strToNode(r1)][strToNode(r2)]
                    (t1, t2, zscore, cr1, cr2, aln1, aln2, score) = alignment_data[cluster_id][strToNode(r1)][strToNode(r2)]

                    if extreme_filtering == True and not is_acceptable_align_len(align_length, align_len_threshold):
                        logger.info('Filtering ' + get_loop_string_for_filtering_log(r1).ljust(75) + ' from ' + str(cluster_id).ljust(max_cid_len) + ' based on length \t(align_length: ' + str(align_length) + ') [Extreme filtering].')
                        removed_loops += 1
                    
                    elif extreme_filtering_based_on_rmsd(rmsd, align_length):
                        logger.info('Filtering ' + get_loop_string_for_filtering_log(r1).ljust(75) + ' from ' + str(cluster_id).ljust(max_cid_len) + ' based on rmsd \t(zscore: ' + "{:.3f}".format(round(zscore, 3)) + ', rmsd: ' + "{:.3f}".format(round(rmsd, 3)) + ') [Extreme filtering].')
                        removed_loops += 1
                    
                    elif remove_based_on_zscore(zscore, rmsd, align_length, is_alignment_from_user):
                        logger.info('Filtering ' + get_loop_string_for_filtering_log(r1).ljust(75) + ' from ' + str(cluster_id).ljust(max_cid_len) + ' based on zscore\t(zscore: ' + "{:.3f}".format(round(zscore, 3)) + ', rmsd: ' + "{:.3f}".format(round(rmsd, 3)) + ').')
                        removed_loops += 1

                    elif remove_based_on_rmsd(rmsd, align_length):
                        logger.info('Filtering ' + get_loop_string_for_filtering_log(r1).ljust(75) + ' from ' + str(cluster_id).ljust(max_cid_len) + ' based on rmsd \t(zscore: ' + "{:.3f}".format(round(zscore, 3)) + ', rmsd: ' + "{:.3f}".format(round(rmsd, 3)) + ').')
                        removed_loops += 1

                    else:
                        filtered_loops.append(r1)

            # filtered_loops = list(set(filtered_loops))
            if len(filtered_loops) > 1:
                filtered_cluster[cluster_id] = filtered_loops
            else:
                removed_loops += len(filtered_loops)
        
        if len(get_loops_in_cluster(filtered_cluster)) == len(get_loops_in_cluster(clusters)):
            break

        clusters = copy.deepcopy(filtered_cluster)

        for cluster_id in clusters:
            loops = clusters[cluster_id]
            current_rmsd_data_dict[cluster_id] = extract_current_rmsd_data_dict(rmsd_data_dict, cluster_id, loops)
        

    if removed_loops > 0:
        logger.info('Removed ' + str(removed_loops) + ' loop' + ('s' if removed_loops > 1 else '') + ' from input data through filtering based on zvalue and rmsd.\n')
    return filtered_cluster

# def filter_loops_in_cluster_denovo_result_analysis(clusters, rmsd_data_dict, alignment_data):

#     current_rmsd_data_dict = copy.deepcopy(rmsd_data_dict)
#     removed_loops = 0
#     while True:
#         filtered_cluster = {}        
#         for c_id in sorted(clusters, key=lambda x: x):
#             loops = clusters[c_id]
#             filtered_loops = []
#             for i in range(len(loops)):
#                 loops[i] = str(strToNode(loops[i]))

#             for cluster_id in sorted(current_rmsd_data_dict, key= lambda x: x):
#                 _, cluster_pairwise_align_details = current_rmsd_data_dict[cluster_id]
#                 align_len_threshold = generate_align_length_threshold(cluster_pairwise_align_details)
#                 for (i, r1) in cluster_pairwise_align_details:
#                     if r1 in loops:
#                         _, pairwise_align_details = cluster_pairwise_align_details[(i, r1)]

#                     j, r2, rmsd, align_length = find_best_aligned_pair(pairwise_align_details, align_len_threshold)
#                         # (t1, t2, zscore, cr1, cr2, aln1, aln2, score) = alignment_data[cluster_id.strip().split('_')[0]][strToNode(r1)][strToNode(r2)]
#                         (t1, t2, zscore, cr1, cr2, aln1, aln2, score) = alignment_data[cluster_id][strToNode(r1)][strToNode(r2)]
#                         if remove_based_on_zscore(zscore, rmsd, align_length):
#                             logger.info('removing ' + r1 + ' from ' + cluster_id + ' based on zscore (zscore: ' + str(zscore) + ', rmsd: ' + str(rmsd) + ').')
#                             removed_loops += 1
#                         elif remove_based_on_rmsd(rmsd, align_length):
#                             logger.info('removing ' + r1 + ' from ' + cluster_id + ' based on rmsd (zscore: ' + str(zscore) + ', rmsd: ' + str(rmsd) + ').')
#                             removed_loops += 1
#                         else:
#                             filtered_loops.append(r1)

#             filtered_loops = list(set(filtered_loops))
#             if len(filtered_loops) > 1:
#                 filtered_cluster[c_id] = filtered_loops
#             else:
#                 removed_loops += len(filtered_loops)
        
#         if len(get_loops_in_cluster(filtered_cluster)) == len(get_loops_in_cluster(clusters)):
#             break

#         clusters = copy.deepcopy(filtered_cluster)

#         for cluster_id in clusters:
#             loops = clusters[cluster_id]
#             current_rmsd_data_dict[cluster_id] = extract_current_rmsd_data_dict(rmsd_data_dict, cluster_id, loops)
        

#     logger.info('Removed Loops from clusters through filtering zvalue and rmsd: ' + str(removed_loops))
#     return filtered_cluster

def generate_length_adjusted_rmsd_score(rmsd_data_dict):
    length_adjusted_rmsd_score_dict = {}
    for cluster_id in rmsd_data_dict:
        _, rmsd_data_list_dict = rmsd_data_dict[cluster_id]
        # total_rmsd = 0.0
        new_rmsd_data_list_dict = {}
        rmsd_align_len_list = []
        for (i, r1) in rmsd_data_list_dict:
            _, pairwise_align_details = rmsd_data_list_dict[(i, r1)]
            rmsd_align_len_list_for_r1 = []
            # total_rmsd_for_r1 = 0.0
            new_pairwise_align_details = []
            for (j, r2, rmsd, align_length) in pairwise_align_details:
                adjusted_score = rmsd / math.sqrt(align_length)
                new_pairwise_align_details.append((j, r2, adjusted_score, align_length))
                # total_rmsd_for_r1 += adjusted_score
                rmsd_align_len_list_for_r1.append((adjusted_score, align_length))

            # avg_rmsd_for_r1 = total_rmsd_for_r1 / len(new_pairwise_align_details)
            avg_rmsd_for_r1, total_align_len_for_r1 = get_weighted_avg_rmsd(rmsd_align_len_list_for_r1)
            new_rmsd_data_list_dict[(i, r1)] = (avg_rmsd_for_r1, sorted(new_pairwise_align_details, key=lambda x: x[2]))
            # total_rmsd += avg_rmsd_for_r1
            rmsd_align_len_list.append((avg_rmsd_for_r1, total_align_len_for_r1))

        # avg_rmsd = total_rmsd / len(new_rmsd_data_list_dict)
        avg_rmsd, total_align_len = get_weighted_avg_rmsd(rmsd_align_len_list)
        length_adjusted_rmsd_score_dict[cluster_id] = (avg_rmsd, new_rmsd_data_list_dict)
    return length_adjusted_rmsd_score_dict

# def extract_current_rmsd_data_dict_denovo_analysis(rmsd_data_dict, loops):

#     for i in range(len(loops)):
#         loops[i] = str(strToNode(loops[i]))

#     new_rmsd_data_list_dict = {}
#     # total_rmsd = 0.0
#     rmsd_align_len_list = []
#     for cluster_id in rmsd_data_dict:
#         _, rmsd_data_list_dict = rmsd_data_dict[cluster_id]
#         not_found_count = 0
#         for (i, r1) in rmsd_data_list_dict:
#             if r1 in loops:
#                 _, pairwise_align_details = rmsd_data_list_dict[(i, r1)]
#                 # total_rmsd_for_r1 = 0.0
#                 rmsd_align_len_list_for_r1 = []
#                 new_pairwise_align_details = []
#                 for (j, r2, rmsd, align_length) in pairwise_align_details:
#                     if r2 in loops:
#                         new_pairwise_align_details.append((j, r2, rmsd, align_length))
#                         rmsd_align_len_list_for_r1.append((rmsd, align_length))
#                         # total_rmsd_for_r1 += rmsd

#                 avg_rmsd_for_r1, total_align_len_for_r1 = get_weighted_avg_rmsd(rmsd_align_len_list_for_r1)
#                 # avg_rmsd_for_r1 = 0.0
#                 # if len(new_pairwise_align_details) > 0:
#                 #     avg_rmsd_for_r1 = total_rmsd_for_r1 / len(new_pairwise_align_details)
#                 new_rmsd_data_list_dict[(i, r1)] = (avg_rmsd_for_r1, sorted(new_pairwise_align_details, key=lambda x: x[2]))
#                 rmsd_align_len_list.append((avg_rmsd_for_r1, total_align_len_for_r1))
#                 # total_rmsd += avg_rmsd_for_r1
#             else:
#                 # print r1 + ' not found'
#                 not_found_count += 1
#                 #sys.exit()

#         # print str(not_fount_count) + ' loops not found.'
#         avg_rmsd, total_align_len = get_weighted_avg_rmsd(rmsd_align_len_list)
#         # avg_rmsd = 0.0
#         # if len(new_rmsd_data_list_dict) > 0:
#         #     avg_rmsd = total_rmsd / len(new_rmsd_data_list_dict)

#     return (avg_rmsd, new_rmsd_data_list_dict)

def extract_current_rmsd_data_dict(rmsd_data_dict, cluster_id, loops):

    for i in range(len(loops)):
        loops[i] = str(strToNode(loops[i]))

    new_rmsd_data_list_dict = {}
    # total_rmsd = 0.0
    rmsd_align_len_list = []
    _, rmsd_data_list_dict = rmsd_data_dict[cluster_id]
    not_found_count = 0
    for (i, r1) in rmsd_data_list_dict:
        if r1 in loops:
            _, pairwise_align_details = rmsd_data_list_dict[(i, r1)]
            # total_rmsd_for_r1 = 0.0
            rmsd_align_len_list_for_r1 = []
            new_pairwise_align_details = []
            for (j, r2, rmsd, align_length) in pairwise_align_details:
                if r2 in loops:
                    new_pairwise_align_details.append((j, r2, rmsd, align_length))
                    rmsd_align_len_list_for_r1.append((rmsd, align_length))
                    # total_rmsd_for_r1 += rmsd

            avg_rmsd_for_r1, total_align_len_for_r1 = get_weighted_avg_rmsd(rmsd_align_len_list_for_r1)
            # avg_rmsd_for_r1 = 0.0
            # if len(new_pairwise_align_details) > 0:
            #     avg_rmsd_for_r1 = total_rmsd_for_r1 / len(new_pairwise_align_details)
            new_rmsd_data_list_dict[(i, r1)] = (avg_rmsd_for_r1, sorted(new_pairwise_align_details, key=lambda x: x[2]))
            rmsd_align_len_list.append((avg_rmsd_for_r1, total_align_len_for_r1))
            # total_rmsd += avg_rmsd_for_r1
        else:
            # print r1 + ' not found'
            not_found_count += 1
            #sys.exit()

    # print str(not_fount_count) + ' loops not found.'
    avg_rmsd, total_align_len = get_weighted_avg_rmsd(rmsd_align_len_list)
    # avg_rmsd = 0.0
    # if len(new_rmsd_data_list_dict) > 0:
    #     avg_rmsd = total_rmsd / len(new_rmsd_data_list_dict)

    return (avg_rmsd, new_rmsd_data_list_dict)

def load_alignment_and_rmsd_data(clusters, loop_list, input_fname_base, partial_pdbx_dir, alignment_dir, graphs_and_pickles_dir, previous_graph_file_reused):
    graph_fname = os.path.join(graphs_and_pickles_dir, input_fname_base + '.z.graph')
    alignment_data = load_alignment_data(input_fname_base, alignment_dir, graphs_and_pickles_dir, graph_fname, clusters, previous_graph_file_reused)
    
    if alignment_data == None:
        return None, None

    rmsd_data_dict = generate_rmsd_data(input_fname_base, partial_pdbx_dir, graphs_and_pickles_dir, alignment_data, clusters, loop_list, previous_graph_file_reused)

    return alignment_data, rmsd_data_dict

# def generate_superimposition_images(removable_text_file_list, partial_pdbx_dir, alignment_dir, superimposition_output_dir, subfamily_details_dir, loop_type, summary_dir, subfamilies_dir, superimposition_details_dir, representative_dir, progressive_dir, graphs_and_pickles_dir, pymol_session_dir, user_input_fname, clusters, loop_list, previous_graph_file_reused, is_alignment_from_user, draw_figures, filter_cluster, set_view_manually, show_extended_loop):
def generate_superimposition_images(clusters, loop_list, alignment_data, rmsd_data_dict, draw_figures, filter_cluster, set_view_manually, show_extended_loop, is_alignment_from_user, user_input_fname, removable_text_file_list, directories, loop_type):

    (partial_pdbx_dir, summary_dir, subfamilies_dir, subfamily_details_dir, representative_dir, superimposition_output_dir, superimposition_details_dir, progressive_dir, pymol_session_dir) = directories

    if draw_figures == True:
        logger.info('Generating superimposition image and output files ...\n')
    else:
        logger.info('Generating superimposition files ...\n')
    
    start_time = time.time()    

    length_adjusted_rmsd_score_dict = rmsd_data_dict
    if is_length_adjusted_score:
        length_adjusted_rmsd_score_dict = generate_length_adjusted_rmsd_score(rmsd_data_dict)

    pdb_organism_details = read_pdb_chain_organism_details(os.path.join(lib_dir, 'PDB_Chain_Organism_Details.tsv'))
    pdb_organism_details_scrapped = read_pdb_chain_organism_details(os.path.join(lib_dir, 'PDB_Chain_Organism_Details_scrapped.tsv'))

    for pdb_chain in pdb_organism_details_scrapped:
        if pdb_chain not in pdb_organism_details:
            pdb_organism_details[pdb_chain] = pdb_organism_details_scrapped[pdb_chain]

    if filter_cluster:
        clusters = filter_loops_in_cluster(clusters, length_adjusted_rmsd_score_dict, alignment_data, is_alignment_from_user)

    include_organism_info = True
    current_rmsd_data_dict = {}
    for cluster_id in clusters:
        loops = clusters[cluster_id]
        for lp in loops:
            pdb_chain = lp.strip().split(':')[0]
            if pdb_chain not in pdb_organism_details:
                include_organism_info = False
                break
        # print(cluster_id,len(loops))
        logger.info('Extracting rmsd data dict for ' + cluster_id)
        current_rmsd_data_dict[cluster_id] = extract_current_rmsd_data_dict(length_adjusted_rmsd_score_dict, cluster_id, loops)
        _, cluster_pairwise_alignment_details = current_rmsd_data_dict[cluster_id]
        logger.info('Completed extracting rmsd data dict for ' + cluster_id)
        # print(len(cluster_pairwise_alignment_details))

    #     print(current_rmsd_data_dict[cluster_id])
    #     print('\n')

    # sys.exit()

    if include_organism_info == False and output_env == 'global':
        pdb_organism_details = {}
    # generate_pymol_images(current_rmsd_data_dict, loop_type, pymol_image_dir, alignment_data, mapping_dir, pdb_dir, aligned_dir, is_cif, is_normalized_score, merge_components, pdb_organism_details, log_file_list, is_length_adjusted_score, draw_pymol_figure)
    time_in_distance_calc = generate_pymol_images(0, removable_text_file_list, partial_pdbx_dir, summary_dir, superimposition_output_dir, subfamily_details_dir, superimposition_details_dir, representative_dir, pymol_session_dir, current_rmsd_data_dict, alignment_data, pdb_organism_details, loop_type, set_view_manually, draw_figures, show_extended_loop)

    if draw_figures == True:
        logger.info('Superimposition image and output file generation complete.')
    else:
        if set_view_manually == True:
            logger.info('View files for the first loop(s) set successfully.')
            return
        else:
            logger.info('Superimposition file generation complete.')

    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.')

    print('\nProcessed input file: ' + os.path.join(data_dir, user_input_fname)[base_path_len:])

    print('Basic configurations:')
    print('Input index type: ' + input_index_type)
    print('Annotation source: ' + annotation_source)
    print('Traversal algorithm: ' + traversal_algorithm)
    print('\nFor generated text outputs, please check the following directories')
    print('==================================================================')
    print('Superimposition details: '.ljust(60) + superimposition_details_dir[base_path_len:])
    print('Annotation of representative motifs: '.ljust(60) + representative_dir[base_path_len:])
    print('Subfamilywise annotations of all motifs: '.ljust(60) + os.path.join(superimposition_output_dir, 'subfamilywise_bp_ann')[base_path_len:])
    print('Subfamily summary and familywise align length threshold: '.ljust(60) + summary_dir[base_path_len:])
    if draw_figures == True:
        print('\nFor generated image outputs, please check the following directories')
        print('===================================================================')
        print('Superimposition outputs: '.ljust(60) + subfamilies_dir[base_path_len:])
        print('Representative motifs: '.ljust(60) + representative_dir[base_path_len:])
        print('Progressive superimposition images: '.ljust(60) + progressive_dir[base_path_len:])
        if output_env == 'local':
            print('')
            print('Time in distance calculation: ' + str(time_in_distance_calc) + ' seconds.')

# # Rotate the first loop of the cluster to define the orientation
def rotate_first_loop_alignto(load_name, rotation_matrix):
    pdb_data = get_pdb_coordinates(load_name)
    pdb_translated = pdb_data# - centroid
    pdb_rotated = numpy.dot(pdb_translated, rotation_matrix)
    alter_structure(pdb_rotated, load_name)

def generate_superimposition_images_using_alignto(superimposition_output_dir, partial_pdbx_dir, clusters, draw_figures):

    if draw_figures == False:
        return

    try:
        import pymol
        from pymol import stored
    except Exception as e:
        try:
            sys.path.append(pymol_py3_dir)
            import pymol
            from pymol import stored
        except Exception as e:
            logger.error('PyMOL not found.')
            sys.exit()

    logger.info('Generating superimposition image files using pymol default alignto.\n')
    start_time = time.time()

    pymol.finish_launching(['pymol', '-cq'])
    # pymol.finish_launching()
    # alignto_methods = ['align', 'super', 'cealign']
    alignto_methods = ['align', 'super']
    # alignto_methods = ['align']
    # alignto_methods = ['super']
    # alignto_methods = ['cealign']   # edit /usr/lib/python2.7/dist-packages/pymol/fitting.py line 30: window=8 => window=5 or 4

    output_dir = os.path.join(superimposition_output_dir, 'alignto_output')
    create_directory(output_dir)

    for cluster_id in clusters:
        loops = clusters[cluster_id]
        # loops = list(set(loops))
        converted_loops = []
        reset_pymol()
        r1 = ''
        for i, loop in enumerate(loops):
            load_color = 'red' if i == 0 else 'green'
            converted_loops.append(convert_a_loop_from_FASTA_to_PDB(loop))
            load_name = 'loop_'+str(i)
            if i == 0:
                r1 = loop
                # print('setting r1 to ' + r1)
                # sys.exit()
            pymol.cmd.load(os.path.join(partial_pdbx_dir, loop+'.cif'), load_name)
            pymol.cmd.hide('everything', load_name)
            pymol.cmd.show('cartoon', load_name)
            pymol.cmd.color(load_color, load_name)

        # pymol.commanding.sync()
        pymol.cmd.sync()
        for method in alignto_methods:
            rotation_matrices = get_multiple_orientation_rotation_matrices()
            for v, rotation_matrix in enumerate(rotation_matrices):
                rotation_version = 'v' + str(v + 1)
                image_fname = os.path.join(superimposition_output_dir, 'alignto_output', cluster_id + '_' + method + '_' + rotation_version + '.png')
                align_to_target = 'loop_0'
                view_fn = os.path.join(views_dir, str(r1) + '.view')
                if os.path.isfile(view_fn):
                    logger.info('View file found for ' + r1 + '. Setting view of this loop for all loops in this cluster.')
                    fv = open(view_fn)
                    view_lines = fv.readlines()
                    fv.close()
                    pymol.cmd.set_view(view_lines[0].strip())
                    # sys.exit()
                # rotate_first_loop_alignto(align_to_target, rotation_matrix)
                # pymol.cmd.show('cartoon', align_to_target)
                pymol.cmd.alignto(align_to_target, method)
                # if os.path.isfile(view_fn):
                #     logger.info('View file found for ' + r1 + '. Setting view of this loop for all loops in this cluster.')
                #     fv = open(view_fn)
                #     view_lines = fv.readlines()
                #     fv.close()
                #     pymol.cmd.set_view(view_lines[0].strip())
                # print(v)
                # pymol.cmd._do('zoom')
                pymol.cmd.zoom()
                pymol.cmd.sync()
                # pymol.cmd._do('set ray_opaque_background, 0')
                # pymol.cmd.set(name='ray_opaque_background',value=0,quiet=1)
                pymol.cmd.png(image_fname, 1200, 1200, dpi=300, ray=1, quiet=1)
                pymol.cmd.sync()

    logger.info('Superimposition image generation (using PyMol \'alignto\') complete.')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.')
