import sys
import logging
import heapq
import numpy
import random
import time
import re
import math
from Bio.PDB import *
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# python 3 compatibility
from functools import reduce

from builtins import input

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from utils import *
from image_helper import *

# if draw_figures == True:
try:
    import pymol
    from pymol import stored
except Exception as e:
    try:
        sys.path.append(pymol_py3_dir)
        import pymol
        from pymol import stored
    except Exception as e:
        pass
        # logger.error('PyMOL not found.')
    
# using z-value for all rmsd
# def generate_merging_rmsd_threshold(cluster_pairwise_alignment_details, rmsd_zscore_threshold=-0.75):
#     rmsd_list = []
#     for (i, r1) in cluster_pairwise_alignment_details:
#         average_rmsd, pairwise_scores = cluster_pairwise_alignment_details[(i, r1)]
#         for j, r2, rmsd, align_len in pairwise_scores:
#             rmsd_list.append(rmsd)
#     mean = get_mean(rmsd_list)
#     std = numpy.std(rmsd_list)
#     merging_rmsd_threshold = round(mean + std * rmsd_zscore_threshold, 4)
#     logger.info('Merging threshold is set to: ' + str(merging_rmsd_threshold))
#     # mean = get_mean(rmsd_list, True)
#     # std = numpy.std(rmsd_list)
#     # merging_rmsd_threshold_temp = round(mean + std * rmsd_zscore_threshold, 4)
#     # logger.info('Merging threshold could be set to (for median): ' + str(merging_rmsd_threshold_temp))
#     return merging_rmsd_threshold

# # previous stable version with condition: align_len_threshold = min(max_align_length * x%, 10)
# def generate_align_length_threshold(cluster_pairwise_alignment_details, threshold_defining_coefficient=0.75):
#     align_len_threshold_upper_limit = 10
#     max_align_len = 0
#     for (i, r1) in cluster_pairwise_alignment_details:
#         average_rmsd, pairwise_scores = cluster_pairwise_alignment_details[(i, r1)]
#         for j, r2, rmsd, align_len in pairwise_scores:
#             if align_len > max_align_len:
#                 max_align_len = align_len
#     # consider doing it using z-value of align lenth distribution
#     align_len_threshold = min(round(max_align_len * threshold_defining_coefficient),align_len_threshold_upper_limit)
#     return align_len_threshold

# using z-value for all align_len
def generate_align_length_threshold(cluster_pairwise_alignment_details):
    if align_len_threshold_type == 'length':
        return align_len_threshold_value

    elif align_len_threshold_type == 'z-score':
        align_len_list = []
        max_align_len = 0
        for (i, r1) in cluster_pairwise_alignment_details:
            average_rmsd, pairwise_scores = cluster_pairwise_alignment_details[(i, r1)]
            for j, r2, rmsd, align_len in pairwise_scores:
                align_len_list.append(align_len)
                if align_len > max_align_len:
                    max_align_len = align_len

        mean = get_mean(align_len_list)
        std = numpy.std(align_len_list)
        align_len_threshold = int(round(mean + std * align_len_zscore_threshold))
        if use_max_align_len_in_equation == True:
            align_len_threshold = min(align_len_threshold, int(round(max_align_len * 0.66)))
        
        align_len_threshold = max(align_len_threshold, min_align_len_threshold)
        # logger.info('Align len threshold is set to: ' + str(align_len_threshold))

        return align_len_threshold

    else:
        logger.error('Invalid align_len_threshold_type. Exiting...')
        sys.exit()

# # using z-value for all align_len and loop length ratio
# def generate_align_length_threshold(cluster_pairwise_alignment_details, align_len_zscore_threshold=0.0):
#     ratio_list = []
#     loop_length_list = []
#     loop_list = []
#     for (i, r1) in cluster_pairwise_alignment_details:
#         loop_len1 = get_loop_length(r1.strip().split(':')[1].strip().split('_'))
#         if r1 not in loop_list:
#             loop_list.append(r1)
#             loop_length_list.append(loop_len1)
#         average_rmsd, pairwise_scores = cluster_pairwise_alignment_details[(i, r1)]
#         for j, r2, rmsd, align_len in pairwise_scores:
#             loop_len2 = get_loop_length(r2.strip().split(':')[1].strip().split('_'))
#             if r2 not in loop_list:
#                 loop_list.append(r2)
#                 loop_length_list.append(loop_len2)
#             ratio1 = align_len / float(loop_len1)
#             ratio2 = align_len / float(loop_len2)
#             ratio_list.append(ratio1)
#             ratio_list.append(ratio2)
#     # print(get_z_scores(align_len_list))
#     # print(get_z_scores(align_len_list, True))
#     # sys.exit()
#     mean = get_mean(ratio_list)
#     # sys.exit()
#     std = numpy.std(ratio_list)
#     align_len_threshold = int(round(get_mean(loop_length_list) * mean))
#     logger.info('Align len threshold is set to: ' + str(align_len_threshold))
#     return align_len_threshold

def is_better_rmsd(rmsd_a, rmsd_b, is_normalized_score = False):
    # Normalized score denotes the higher the better score (reverse of RMSD property)
    if is_normalized_score:
        if compare_rmsd(rmsd_a, rmsd_b) < 0:
            return True
        return False
    if compare_rmsd(rmsd_a, rmsd_b) > 0:
        return True
    return False

#For RMSD, the lower, the better
#Returns 1 if better RMSD, -1 for worse and 0 for equal
def compare_rmsd(rmsd_a, rmsd_b):
    prec = 0.000001
    # definitively better RMSD
    if rmsd_a < rmsd_b - prec:
        return 1

    # definitively worse RMSD
    if rmsd_a > rmsd_b + prec:
        return -1

    return 0

def is_better_alignment_score_length_priority(rmsd_a, alignment_length_a, rmsd_b, alignment_length_b, is_normalized_score = False):
    if alignment_length_a > alignment_length_b:
        return True
    if alignment_length_a < alignment_length_b:
        return False
    return is_better_rmsd(rmsd_a, rmsd_b, is_normalized_score)

def is_better_alignment_score_RMSD_priority(rmsd_a, alignment_length_a, rmsd_b, alignment_length_b, is_normalized_score):
    rmsd_comp = compare_rmsd(rmsd_a, rmsd_b)

    # Normalized score denotes the higher, the better score (reverse of RMSD property)
    if is_normalized_score:
        if rmsd_comp < 0:
            return True
        if rmsd_comp > 0:
            return False
    else:
        if rmsd_comp > 0:
            return True
        if rmsd_comp < 0:
            return False
    return alignment_length_a > alignment_length_b

def is_acceptable_align_len(align_len, align_len_threshold):
    if align_len >= align_len_threshold:
        return True
    return False

def is_acceptable_rmsd(rmsd, align_len, is_length_adjusted_score, merging_rmsd_threshold):
    if is_length_adjusted_score:
        rmsd = rmsd * math.sqrt(align_len)

    if rmsd <= merging_rmsd_threshold:
        return True
    return False

def extract_alignment_scores(scores):
    alignment_length = alignment_score = 0
    rmsd = zscore = 0.0
    if len(scores) == 1:
        rmsd = scores
    elif len(scores) == 2:
        rmsd, alignment_length = scores
    elif len(scores) == 3:
        rmsd, alignment_length, zscore = scores
    elif len(scores) == 4:
        rmsd, alignment_length, zscore, alignment_score = scores
    else:
        # print(alignment)
        logger.error('Invalid "scores" value.')
    
    return rmsd, alignment_length, zscore, alignment_score

#scores can be (rmsd, [align_length, [zscore, [alignment_score]]])
def is_better_alignment_score(scores_a, scores_b, align_len_threshold, is_normalized_score, priority = 'balance'):
    rmsd_a, align_length_a, zscore_a, alignment_score_a = extract_alignment_scores(scores_a)
    rmsd_b, align_length_b, zscore_b, alignment_score_b = extract_alignment_scores(scores_b)

    # if priority == 'balance':
    if is_acceptable_align_len(align_length_a, align_len_threshold) and is_acceptable_align_len(align_length_b, align_len_threshold):
        return is_better_alignment_score_RMSD_priority(rmsd_a, align_length_a, rmsd_b, align_length_b, is_normalized_score)
    else:
        return is_better_alignment_score_length_priority(rmsd_a, align_length_a, rmsd_b, align_length_b, is_normalized_score)
    
    # Other options can be added (not considered here yet)

def find_best_aligned_pair(pairwise_align_details, align_len_threshold):
    # max align loop

    max_align_length = max_loop_index = 0
    max_length_loop = ''
    max_loop_rmsd = 1000.0

    for (j, r2, rmsd, align_length) in pairwise_align_details:
        if is_better_alignment_score((rmsd, align_length), (max_loop_rmsd, max_align_length), align_len_threshold, is_normalized_score):
            max_align_length = align_length
            max_length_loop = r2
            max_loop_index = j
            max_loop_rmsd = rmsd
    return max_loop_index, max_length_loop, max_loop_rmsd, max_align_length

def get_weighted_avg_rmsd(rmsd_align_len_list):
    total_alignment_length = 0
    total_rmsd = 0.0
    for rmsd, align_len in rmsd_align_len_list:
        total_alignment_length += align_len

        if is_length_adjusted_score == True:
            rmsd = rmsd * math.sqrt(align_len)

        total_rmsd += rmsd * rmsd * align_len

    avg_rmsd = 0.0
    if total_alignment_length > 0:
        avg_rmsd = math.sqrt(total_rmsd / total_alignment_length)
        if is_length_adjusted_score == True:
            avg_rmsd /= math.sqrt(total_alignment_length)

    return avg_rmsd, total_alignment_length

def centroid(coord_list):
    if len(coord_list) > 0:
        return list(map(lambda z: 1.*z/len(coord_list), reduce(lambda x, y: (x[0]+y[0], x[1]+y[1], x[2]+y[2]), coord_list)))
    return None

def get_atom_coordinate(pdb_fn, residue_list):

    backbone_atoms, sugar_atoms = get_backbone_and_sugar_atoms()
    pdb_id = os.path.basename(pdb_fn)[:4]

    parser = FastMMCIFParser()
    structure = parser.get_structure('struct', pdb_fn)

    backbone = {}
    sugar = {}

    for chain_id, index, icd in residue_list:
        # if chain_id == 'n' and index == 'a':
        if chain_id == '':
            continue
        chain = structure[0][chain_id]
        residues = chain.get_residues()

        my_residues = {}

        for r in residues:
            hetflag, resseq, icode = r.get_id()
            my_residues[(resseq, icode)] = r

        i = int(index)
        icode = icd if len(icd) > 0 else ' '

        if (i, icode) not in my_residues:
            # ret.append(0)
            backbone[(pdb_id, chain_id, index, icd)] = 0.
            sugar[(pdb_id, chain_id, index, icd)] = 0.
        else:
            atom_coord = []
            for atom in backbone_atoms:
                if atom in my_residues[(i, icode)]:
                    atom_coord.append(my_residues[(i, icode)][atom].get_vector())

            backbone[(pdb_id, chain_id, index, icd)] = centroid(atom_coord)

            atom_coord = []
            for atom in sugar_atoms:
                if atom in my_residues[(i, icode)]:
                    atom_coord.append(my_residues[(i, icode)][atom].get_vector())

            sugar[(pdb_id, chain_id, index, icd)] = centroid(atom_coord)

    return backbone, sugar, structure

def pdb_pos_map(pdb_res_map, m):
    """position in pdb index alignment"""
    ret = []
    for i in m:
        if i in pdb_res_map:
            ret.append(pdb_res_map[i])
        # if 
        else:
            # ret.append('na')
            ret.append(('', '', ''))
            logger.warning('!!!!!!!!!!!!!!!!!!!!!ALERT: APPENDING EMPTY TUPLE (NA) !!!!!!!!!!!!!!!!!!!!')

    return ret

def get_pdb_index_list(lp):
    pdb_chain, regions = lp.split(':')
    # segments = regions.strip().split('_')
    # index_list = []
    r = list(map(lambda x: x.split('-'), regions.split('_')))
    index_list = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0]), int(x[1])+1)), r)))
    pdb_res_map = load_pdb_res_map(pdb_chain)
    pdb_index_list = pdb_pos_map(pdb_res_map, index_list)

    return pdb_index_list

def aln_map(aln1, aln2, aln_start, aln_end):
    """the aligned nucleotides in two regions"""
    """return index in the aligned regions"""
    if len(aln1) != len(aln2):
        return None

    aln1 = aln1.replace('.', '')
    aln2 = aln2.replace('.', '')
    aln1 = aln1.replace('~', '')
    aln2 = aln2.replace('~', '')
    # print 'aln1: ' + aln1
    # print 'aln2: ' + aln2
    ret1 = []
    ret2 = []
    i = j = 0
    for index, (c1, c2) in enumerate(zip(aln1, aln2)):
        if c1 == '-' and c2 != '-':
            j += 1
        elif c1 != '-' and c2 == '-':
            i += 1
        elif c1 != '-' and c2 != '-':
            if index in range(aln_start, aln_end+1):
                ret1.append(i)
                ret2.append(j)
            i += 1
            j += 1
        else:
            return None
    # print 'ret1: ' + ret1
    # print 'ret2: ' + ret2
    return ret1, ret2

def pos_map(region, m, extend):
    """region is based on ref_index"""
    """m is the output of aln_map"""
    """position in sequence alignment"""
    r = list(map(lambda x: x.split('-'), region.split('_')))

    i = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0]), int(x[1])+1)), r)))
    # i = reduce(lambda x, y: range(int(x[0]), int(x[1])+1)+range(int(y[0]), int(y[1])+1), r)
    i_e = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0])-extend, int(x[1])+1+extend)), r)))
    # i_e = reduce(lambda x, y: range(int(x[0])-extend, int(x[1])+1+extend)+range(int(y[0])-extend, int(y[1])+1+extend), r)
    ret = []
    for j in m:
        ret.append(i[j])
    return ret, i_e

def pos_map_is_safe(region, m, extend):
    """region is based on ref_index"""
    """m is the output of aln_map"""
    """position in sequence alignment"""
    r = list(map(lambda x: x.split('-'), region.split('_')))

    i = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0]), int(x[1])+1)), r)))
    # i = reduce(lambda x, y: range(int(x[0]), int(x[1])+1)+range(int(y[0]), int(y[1])+1), r)
    i_e = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0])-extend, int(x[1])+1+extend)), r)))
    # i_e = reduce(lambda x, y: range(int(x[0])-extend, int(x[1])+1+extend)+range(int(y[0])-extend, int(y[1])+1+extend), r)
    ret = []
    # print i
    # print m
    if m[len(m)-1] >= len(i):
    # if len(i) != len(m):
        return False
    return True

def aln_residue(r1, r2, loop1, loop2, aln1, aln2, aln_start, aln_end, extend):
    chain1, region1 = loop1.split(':')
    chain2, region2 = loop2.split(':')

    # load the ref_index to pdb_index mapping
    pdb1_res_map = load_pdb_res_map(chain1)
    pdb2_res_map = load_pdb_res_map(chain2)

    # find the index of mapping nucleotides
    am1, am2 = aln_map(aln1, aln2, aln_start, aln_end)

    (pm1, i1) = pos_map(region1, am1, extend)
    (pm2, i2) = pos_map(region2, am2, extend)

    # get the ref_index for aligned nucleotides
    (pm1, i1) = pos_map(region1, am1, extend)
    (pm2, i2) = pos_map(region2, am2, extend)

    # get the pdb_index for aligned nucleotides
    pdb1_pm = pdb_pos_map(pdb1_res_map, pm1)
    pdb2_pm = pdb_pos_map(pdb2_res_map, pm2)
    i1_pm = pdb_pos_map(pdb1_res_map, i1)
    i2_pm = pdb_pos_map(pdb2_res_map, i2)

    return pdb1_pm, pdb2_pm, i1_pm, i2_pm

def write_alignment_fixing_log(fw, fix_cnt, r1, r2, loop1, loop2, aln1, aln2, status):
    if output_env == 'local':
        if status == 'before':
            fw.write('fix count: ' + str(fix_cnt) + '\n')
            fw.write(r1 + '\n')
            fw.write(r2 + '\n')
            fw.write('Aligned region(s):\n')
            fw.write(loop1 + '\n')
            fw.write(loop2 + '\n')
            fw.write('Alignment string BEFORE fixing:' + '\n')
            fw.write(aln1 + '\n')
            fw.write(aln2 + '\n')
        else:
            fw.write('Alignment string AFTER fixing:' + '\n')
            fw.write(aln1 + '\n')
            fw.write(aln2 + '\n\n')

def remove_hiphen_and_corresponding_character(aln1, aln2):
    new_aln1 = ''
    new_aln2 = ''
    for i in range(len(aln1)):
        if aln1[i] != '-' and aln2[i] != '-':
            new_aln1 += aln1[i]
            new_aln2 += aln2[i]

    return new_aln1, new_aln2

def get_segment_fasta_seq(fasta_seq, regions):
    seq = ''
    regions = regions.strip().split('_')
    for region in regions:
        s, e = region.strip().split('-')
        s = int(s)
        e = int(e)
        for i in range(s, e+1):
            seq += fasta_seq[i]

    return seq

def fix_alignment(fasta_seq, aln1, aln2):
    new_aln1 = list(aln1)
    new_aln2 = list(aln2)

    last_hyphen_ind = -1
    for i in range(len(aln1)):
        
        if aln2[i] == '-':
            last_hyphen_ind = i

        if fasta_seq[i] != aln1[i]:
            j = last_hyphen_ind
            new_aln1.pop(j)
            new_aln2.pop(j)
            break

    return ''.join(new_aln1), ''.join(new_aln2)

def fix_alignment_pair(fasta_seq1, fasta_seq2, aln1, aln2):

    if '#' in fasta_seq1:
        new_aln1, new_aln2 = fix_alignment(fasta_seq1, aln1, aln2)

    if '#' in fasta_seq2:
        new_aln2, new_aln1 = fix_alignment(fasta_seq2, aln2, aln1)

    return new_aln1, new_aln2

def aln_residue_temp(pdb_res_mapping_dict, fasta_seq_dict, r1, r2, loop1, loop2, aln1, aln2, aln_start, aln_end, extend):

    if output_env == 'local':
        fw = open('alignment_fixing_log.log', 'a')
    # print('loops: ')
    # print(loop1)
    # print(loop2)
    chain1, region1 = loop1.split(':')
    chain2, region2 = loop2.split(':')

    # load the ref_index to pdb_index mapping
    # pdb1_res_map = load_pdb_res_map(chain1)
    # pdb2_res_map = load_pdb_res_map(chain2)
    pdb1_res_map = pdb_res_mapping_dict[chain1]
    pdb2_res_map = pdb_res_mapping_dict[chain2]

    # find the index of mapping nucleotides
    am1, am2 = aln_map(aln1, aln2, aln_start, aln_end)
    # print am1, am2
    # sys.exit()

    fix_cnt = 0
    while pos_map_is_safe(region1, am1, extend) == False or pos_map_is_safe(region2, am2, extend) == False:
        loop1_expected_length = get_loop_length(region1.strip().split('_'))
        loop2_expected_length = get_loop_length(region2.strip().split('_'))
        # fix alignment
        fix_cnt += 1
        
        write_alignment_fixing_log(fw, fix_cnt, r1, r2, loop1, loop2, aln1, aln2, 'before')

        fasta_seq1 = get_segment_fasta_seq(fasta_seq_dict[chain1], region1)
        fasta_seq2 = get_segment_fasta_seq(fasta_seq_dict[chain2], region2)

        hyphen_free_seq1 = aln1.replace('-', '')
        hyphen_free_seq2 = aln2.replace('-', '')

        fasta_seq1 = fasta_seq1.ljust(len(hyphen_free_seq1), '#')
        fasta_seq2 = fasta_seq2.ljust(len(hyphen_free_seq2), '#')

        aln1, aln2 = fix_alignment_pair(fasta_seq1, fasta_seq2, aln1, aln2)
        
        write_alignment_fixing_log(fw, fix_cnt, r1, r2, loop1, loop2, aln1, aln2, 'after')

        am1, am2 = aln_map(aln1, aln2, aln_start, aln_end)

        if fix_cnt > max(len(aln1), len(aln2)):
            break

    if pos_map_is_safe(region1, am1, extend) == False or pos_map_is_safe(region2, am2, extend) == False:
        logger.error('Alignment length missmatch fixing failed. r1: ' + r1 + ', r2: ' + r2)
        logger.error(aln1)
        logger.error(aln2)
        sys.exit()
        
    (pm1, i1) = pos_map(region1, am1, extend)
    (pm2, i2) = pos_map(region2, am2, extend)

    # get the ref_index for aligned nucleotides
    (pm1, i1) = pos_map(region1, am1, extend)
    (pm2, i2) = pos_map(region2, am2, extend)

    # if loop1 == '3J7Q_5:3277-3280' and loop2 == '4V9F_0:1699-1703':
    #    sys.exit()
    # get the pdb_index for aligned nucleotides
    pdb1_pm = pdb_pos_map(pdb1_res_map, pm1)
    pdb2_pm = pdb_pos_map(pdb2_res_map, pm2)
    i1_pm = pdb_pos_map(pdb1_res_map, i1)
    i2_pm = pdb_pos_map(pdb2_res_map, i2)

    if output_env == 'local':
        fw.close()

    return pdb1_pm, pdb2_pm, i1_pm, i2_pm

def load_fasta_seq(pdb_id, chains):
    fasta_seq_dict = {}
    fasta_fn = os.path.join(fasta_dir, pdb_id + '.fasta')
    for record in SeqIO.parse(fasta_fn, 'fasta'):
        # fasta_seq_dict[pdb_id + '_' + record.id.strip().split('|')[0].strip().split(':')[1]] = str(record.seq)
        chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1].strip().split(',')
        for chain_id in chain_ids:
            fasta_seq_dict[chain_id] = str(record.seq)

    return fasta_seq_dict

def load_pdb_res_map(chain):
    """load sequence index->pdb index"""
    """{ref_index: (chain_id, pdb_index)}"""
    ret = {}
    # map_dir = '../nrPDBs_Old' # the directory for the mapping data
    fp = open(os.path.join(pdb_fasta_mapping_dir, chain+'.rmsx.nch'))
    for line in fp.readlines():
        decom = line.strip().split('\t')
        ##################### for PDB #####################
        # if decom[0][0] == "'":
        #     ret[int(decom[1])] = (decom[0][1], decom[0][3:].replace('.', ''))
        # else:
        #     ret[int(decom[1])] = (decom[0][0], decom[0][1:].replace('.', ''))
        ##################### for PDB #####################
        ##################### for PDBx ####################
        if decom[0][0] == "'":
            chain_id = decom[0][1:].strip().split("'")[0]
            i = len(chain_id)+2
        else:
            chain_id = re.split('-?(\d+)',decom[0])[0]
            i = len(chain_id)

        if decom[0][-1].isalpha():
            icode = decom[0][-1]
            j = len(decom[0])-2
        else:
            icode = ''
            j = len(decom[0])

        seqnum = decom[0][i:j]
        ret[int(decom[1])] = (chain_id, seqnum, icode)
        ##################### for PDBx ####################
    return ret

def create_best_alignment_graph(cluster_pairwise_alignment_details, align_len_threshold):
    adjacency_list = {}
    directed_adjacency_list = {}
    transpose_directed_adjacency_list = {}

    # Create undirected, weighted directed, transpose-directed graph from best alignment edge among loops
    for (i, r1) in sorted(cluster_pairwise_alignment_details, key=lambda x: cluster_pairwise_alignment_details[x][0]):
        adjacency_list[(i, r1)] = []
        directed_adjacency_list[(i, r1)] = {}
        transpose_directed_adjacency_list[(i, r1)] = []
    for (i, r1) in sorted(cluster_pairwise_alignment_details, key=lambda x: cluster_pairwise_alignment_details[x][0]):
        avg_rmsd, pairwise_align_details = cluster_pairwise_alignment_details[(i, r1)]
        if len(pairwise_align_details) > 0:
            j, r2, max_loop_rmsd, max_align_length = find_best_aligned_pair(pairwise_align_details, align_len_threshold)
            adjacency_list[(i, r1)].append((j, r2))
            adjacency_list[(j, r2)].append((i, r1))
            directed_adjacency_list[(j, r2)][(i, r1)] = (max_loop_rmsd, max_align_length)
            transpose_directed_adjacency_list[(i, r1)].append((j, r2))

    return adjacency_list, directed_adjacency_list, transpose_directed_adjacency_list

def create_weighted_complete_directed_graph(cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score):
    # Create weighted directed graph from alignment edge among loops
    weighted_directed_adjacency_matrix = {}

    for (i, r1) in sorted(cluster_pairwise_alignment_details, key=lambda x: cluster_pairwise_alignment_details[x][0]):
        weighted_directed_adjacency_matrix[(i, r1)] = {}

    for (i, r1) in sorted(cluster_pairwise_alignment_details, key=lambda x: cluster_pairwise_alignment_details[x][0]):
        avg_rmsd, pairwise_align_details = cluster_pairwise_alignment_details[(i, r1)]
        for (j, r2, rmsd, align_length) in pairwise_align_details:
            weighted_directed_adjacency_matrix[(j, r2)][(i, r1)] = (rmsd, align_length)

    return weighted_directed_adjacency_matrix

def get_connected_components(adjacency_list):
    # Find the connected components of the graph
    connected_components = []
    all_visited = []
    for (i, r1) in adjacency_list:
        if (i, r1) not in all_visited:
            current_visted = []
            dfs(adjacency_list, (i,r1), current_visted, None, [])
            connected_components.append(current_visted)
            all_visited.extend(current_visted)
    return connected_components

def extract_cycle(traverse_list_with_parent, node, initial_node, cycle):
    if node == initial_node:
        return
    cycle.append(node)
    extract_cycle(traverse_list_with_parent, traverse_list_with_parent[node], initial_node, cycle)

# Take the transpose of the directed connected graph and find dependency order (on the best alignment component)
# Feature 1: each node will have exactly one out-edge in the reverse graph
# Feature 1 implies one node can be part of at most one cycle
# Feature 2: There will always be exactly 1 cycle
# Feature 2 can be proved by contradiction
# Contradiction Case 1 (no cycle): n nodes, n edge cannot for a tree, where (n-1) can be accomodated
# Contradiction Case 2 (2 or more cycles): cycles can have shared node, can't have connected edge
def bfs_find_cycle(graph, start):
    visited = []
    traverse_list_with_parent = {}
    parent = None
    queue = [(start, parent)]
    while queue:
        (node, parent) = queue.pop(0)
        if node not in visited:
            visited.append(node)
            traverse_list_with_parent[node] = parent
            for n in sorted(graph[node], key = lambda x:len(graph[x])):
                queue.append((n, node))
        else:
            cycle = [node]
            extract_cycle(traverse_list_with_parent, parent, node, cycle)
            return cycle
    return []

def get_central_node_in_the_component(component, directed_adjacency_list, cluster_pairwise_alignment_details):
    # There will always be one in-degree, but 0+ outdegree
    # Find the best node based on the out degree, and then avg rmsd
    central_node = ''
    central_node_out_degree = 0
    central_node_avg_rmsd = 100.0

    for loop in component:
        out_degree = len(directed_adjacency_list[loop])

        rmsd_align_len_list = []
        for loop2 in directed_adjacency_list[loop]:
            rmsd_align_len_list.append(directed_adjacency_list[loop][loop2])
        avg_rmsd, total_align_len = get_weighted_avg_rmsd(rmsd_align_len_list)
        # avg_rmsd = cluster_pairwise_alignment_details[loop]

        if out_degree > central_node_out_degree or (out_degree == central_node_out_degree and avg_rmsd < central_node_avg_rmsd):
            central_node = loop
            central_node_out_degree = out_degree
            central_node_avg_rmsd = avg_rmsd

    return central_node

def get_component_features(cluster_pairwise_alignment_details, directed_adjacency_list, transpose_directed_adjacency_list, connected_components):
    ### Find features (central_node, cycle_nodes_of_the_component, component_nodes) of components
    connected_components_features = {}

    for component_id, component_nodes in enumerate(sorted(connected_components, key=lambda x:len(connected_components), reverse = True)):
        ## Find the nodes that can be drawn without dependency
        start_node = component_nodes[0]
        # Returns the cycle in reverse order
        cycle_nodes_of_the_component = bfs_find_cycle(transpose_directed_adjacency_list, start_node)

        # extract the adjaceny list of the current component
        component_directed_adjacency_list = {}
        for (i, r1) in directed_adjacency_list:
            if (i, r1) in component_nodes:
                component_directed_adjacency_list[(i, r1)] = directed_adjacency_list[(i, r1)]

        central_node = start_node
        ## If the graph has more than one node, there will the a cycle and we can find a sutable node_to_start there
        if len(cycle_nodes_of_the_component) > 0:
            
            if len(cycle_nodes_of_the_component) > 2 and output_env == 'local':
                logger.info('### CYCLE FOUND WITH MORE THAN TWO MOTIFS ###')
                print(component_id)
                print('Cycle nodes: ')
                print(cycle_nodes_of_the_component)
                sys.exit()

            # Find the most important loop in the cycle of the current connected component
            # Returns the node that has most outgoing edge, then best average rmsd
            central_node = get_central_node_in_the_component(cycle_nodes_of_the_component, component_directed_adjacency_list, cluster_pairwise_alignment_details)
        connected_components_features[component_id] = (central_node, cycle_nodes_of_the_component, component_nodes, component_directed_adjacency_list)
    return connected_components_features

def find_best_edge_between_components(cycle_nodes_of_the_component2, component1, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score):
    best_align_length = 0
    best_align_rmsd = 100.0
    best_align_loop_source = ''
    best_align_loop_source_index = 0
    best_align_loop_dest = ''
    best_align_loop_dest_index = 0

    for (i, r1) in cycle_nodes_of_the_component2:
        _, pairwise_align_details = cluster_pairwise_alignment_details[(i,r1)]
        new_pairwise_align_details = []
        for (j, r2, rmsd, align_length) in pairwise_align_details:
            if (j, r2) in component1:
                new_pairwise_align_details.append((j, r2, rmsd, align_length))
        best_j, best_r2, rmsd, align_length = find_best_aligned_pair(new_pairwise_align_details, align_len_threshold)

        if is_better_alignment_score((rmsd, align_length), (best_align_rmsd, best_align_length), align_len_threshold, is_normalized_score):
            best_align_length = align_length
            best_align_rmsd = rmsd
            best_align_loop_source = best_r2
            best_align_loop_source_index = best_j
            best_align_loop_dest = r1
            best_align_loop_dest_index = i

    return ((best_align_loop_source_index, best_align_loop_source), (best_align_loop_dest_index, best_align_loop_dest), best_align_length, best_align_rmsd)

def find_best_edges_to_another_component(component1, component2, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score):
    best_edges = {}
    for (i, r1) in component1:
        _, pairwise_align_details = cluster_pairwise_alignment_details[(i, r1)]
        new_pairwise_align_details = []
        for (j, r2, rmsd, align_length) in pairwise_align_details:
            if (j, r2) in component2:
                new_pairwise_align_details.append((j, r2, rmsd, align_length))
        
        best_j, best_r2, best_align_rmsd, best_align_length = find_best_aligned_pair(new_pairwise_align_details, align_len_threshold)
        best_edges[(i, r1)] = (best_j, best_r2), best_align_rmsd, best_align_length

    return best_edges

# Implements Prim's MST approach
def get_component_ordering(component_adj_matrix, start, align_len_threshold, is_normalized_score):
    visited = []
    component_ordering = []
    best_next_node = start

    best_edge_loop = None
    best_edge_parent = None
    best_align_length = 0
    best_align_loop_rmsd = 0.0

    while best_next_node != None:
        visited.append(best_next_node)
        component_ordering.append((best_next_node, best_edge_loop, best_edge_parent, best_align_loop_rmsd, best_align_length))
        best_next_node = None
        best_align_length = 0
        best_align_loop_rmsd = 100.0

        #Find which edge is the best among all the outgoing edge from the already visited nodes
        #Time efficiency was not considered as the number of nodes (component) is small
        for i in visited:
            for j in component_adj_matrix[i]:
                if j not in visited:
                    (loop_source, loop_dest, align_length, rmsd) = component_adj_matrix[i][j]
                    # if align_len > best_align_len or (align_len == best_align_len and is_better_rmsd(rmsd, best_align_loop_rsmd, is_normalized_score)):
                    if is_better_alignment_score((rmsd, align_length), (best_align_loop_rmsd, best_align_length), align_len_threshold, is_normalized_score):
                        best_align_length = align_length
                        best_align_loop_rmsd = rmsd
                        best_next_node = j
                        best_edge_loop = loop_dest
                        best_edge_parent = loop_source

    return component_ordering

# Test if at least given percentage/count of the best edges conecting components passes the merging threshold
def assess_merging_two_components(best_edges, align_len_threshold):
    acceptable_edge_rmsd = []
    for (i, r1) in best_edges:
        (j, r2), rmsd, align_len = best_edges[(i,r1)]
        if merge_components and acceptable_edge_for_merging(rmsd, align_len, align_len_threshold, is_length_adjusted_score):
            acceptable_edge_rmsd.append(rmsd)
    
    # For the case of connectivity_test_type = "count", consider the given value or # of loops
    match_count = min(connectivity_test_threshold, len(best_edges))

    if connectivity_test_type == "percent":
        match_count = math.ceil(len(best_edges)*(connectivity_test_threshold/100.0))

    if len(acceptable_edge_rmsd) < match_count:
        # Returning really bad RMSD
        return False, 1000.0
    else:        
        RMSD_sum = sum(sorted(acceptable_edge_rmsd)[:match_count])
        if match_count > 0:
            return True, RMSD_sum/match_count
        else:
            return True, 0.0

# Apply guided-tree based approach to merge the components
def merge_and_order_connected_components(connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score):

    ### Find relational features among components
    component_count = len(connected_components_features)
    component_adj_matrix = {}
    for i in connected_components_features:
        component_adj_matrix[i] = {}
        for j in connected_components_features:
            component_adj_matrix[i][j] = []

    component_nodes_dict = {}
    ## Build graph with pair-wise edge of the best loop-alignment choice among the components
    for i in connected_components_features:
        _, _, component_nodes1, _ = connected_components_features[i]
        component_nodes_dict[i] = copy.deepcopy(component_nodes1)
        for j in connected_components_features:
            if i != j:
                _, _, component_nodes2, _ = connected_components_features[j]
                # Find the best edge from all nodes in the component1 to any node in the component 2
                best_edges = find_best_edges_to_another_component(component_nodes1, component_nodes2, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)
                component_adj_matrix[i][j] = best_edges

    component_merged = {}
    while(True):
        got_mergable_components = False
        mergable_component_id_min = None
        mergable_component_id_max = None
        best_component_merge_avg_RMSD = 1000.0
        for i in connected_components_features:
            if i not in component_merged:
                for j in range(i+1, len(connected_components_features)):
                    if j not in component_merged:
                        is_mergable1, component_merge_avg_RMSD1 = assess_merging_two_components(component_adj_matrix[i][j], align_len_threshold)
                        is_mergable2, component_merge_avg_RMSD2 = assess_merging_two_components(component_adj_matrix[j][i], align_len_threshold)

                        if (is_mergable1 and is_mergable2) or (connectivity_direction == "one-way" and (is_mergable1 or is_mergable2)):
                            got_mergable_components = True
                            component_merge_avg_RMSD = (component_merge_avg_RMSD1 + component_merge_avg_RMSD2) / 2.0
                            if connectivity_direction == "one-way":
                                component_merge_avg_RMSD = min(component_merge_avg_RMSD1, component_merge_avg_RMSD2)
                                if not (is_mergable1 and is_mergable2):
                                    component_merge_avg_RMSD = component_merge_avg_RMSD1 if is_mergable1 else component_merge_avg_RMSD2

                            if component_merge_avg_RMSD <  best_component_merge_avg_RMSD:
                                best_component_merge_avg_RMSD = component_merge_avg_RMSD
                                mergable_component_id_min = i
                                mergable_component_id_max = j

        if got_mergable_components == False:
            break
        else:
            # Merge two components
            component_merged[mergable_component_id_max] = mergable_component_id_min
            component_nodes_dict[mergable_component_id_min] += component_nodes_dict[mergable_component_id_max]            
            component1 = component_nodes_dict[mergable_component_id_min]

            # Update the matrix
            for i in connected_components_features:
                if i != mergable_component_id_min and (i not in component_merged):
                    best_edges1 = find_best_edges_to_another_component(component1, component_nodes_dict[i], cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)
                    component_adj_matrix[mergable_component_id_min][i] = best_edges1
                    best_edges2 = find_best_edges_to_another_component(component_nodes_dict[i], component1, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)
                    component_adj_matrix[i][mergable_component_id_min] = best_edges2

    merged_components_dict = {}
    
    # Include the components that was not merged
    for i in sorted(connected_components_features):
        if i not in component_merged:
            merged_components_dict[i] = []
            merged_components_dict[i].append(i)

    # Complete the list with the merged components     
    for i in sorted(component_merged):
        key = component_merged[i]
        # Find the smallest component id that it is connected to (the component that was not merged)
        while key in component_merged:
            key = component_merged[key]
        merged_components_dict[key].append(i)

    merged_components_features = {}
    for i in sorted(merged_components_dict):
        merged_components_nodes = []
        merged_cycle_nodes_of_the_components = []
        for j in sorted(merged_components_dict[i]):
            _, cycle_nodes_of_the_component, component_nodes, _ = connected_components_features[j]
            merged_cycle_nodes_of_the_components += cycle_nodes_of_the_component
            merged_components_nodes += component_nodes            

        merged_components_features[i] = (merged_cycle_nodes_of_the_components, merged_components_nodes)

    start_node_list = None
    if start_traversal_with_largest_subfamily == True:
        start_node_list = get_largest_subfamily(merged_components_features, cluster_pairwise_alignment_details)

    merged_component_ordering = get_best_component_ordering(merged_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score, start_node_list)

    merged_component_ordering_list = []
    
    first_component = True
    for a_component in merged_component_ordering:
        i, best_edge_next_loop, best_edge_parent, best_align_loop_rmsd, best_align_length = a_component
        if len(merged_components_dict[i]) > 1:
            next_loop_component = None
            next_loop_component_list = None
            filtered_connected_components_features = {}
            for j in merged_components_dict[i]:
                _, cycle_nodes_of_the_component, component_nodes, _ = connected_components_features[j]
                # Find which is the next component to traverse given the next loop we know
                if first_component == False and best_edge_next_loop in component_nodes:
                    next_loop_component = j
                    next_loop_component_list = [j]
                filtered_connected_components_features[j] = cycle_nodes_of_the_component, component_nodes

            if start_traversal_with_largest_subfamily == True and first_component == True:
                next_loop_component_list = get_largest_subfamily(filtered_connected_components_features, cluster_pairwise_alignment_details)
                
            component_ordering = get_best_component_ordering(filtered_connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score, next_loop_component_list)
            if first_component == False:
                component_ordering[0] = (next_loop_component, best_edge_next_loop, best_edge_parent, best_align_loop_rmsd, best_align_length)
            merged_component_ordering_list.append(component_ordering)
        else:
            merged_component_ordering_list.append([a_component])    
        first_component = False

    return merged_component_ordering_list

def get_largest_subfamily(connected_components_features, cluster_pairwise_alignment_details):
    largest_subfamily_ind_list = None
    largest_node_count = -1
    for i in connected_components_features:
        _, component_nodes = connected_components_features[i]
        if len(component_nodes) > largest_node_count:
            largest_node_count = len(component_nodes)
            largest_subfamily_ind_list = [i]
        elif len(component_nodes) == largest_node_count:
            largest_subfamily_ind_list.append(i)

    return largest_subfamily_ind_list

def get_best_component_ordering(connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score, start_node_list = None):
    ### Find relational features among components
    component_count = len(connected_components_features)
    component_adj_matrix = {}
    for i in connected_components_features:
        component_adj_matrix[i] = {}
        for j in connected_components_features:
            component_adj_matrix[i][j] = ((0,''),(0,''),0,0.0)
    best_edge_index_i = 0
    overall_best_align_length = 0
    overall_best_align_rmsd = 100.0

    ## Build graph with pair-wise edge of the best loop-alignment choice among the components
    for i in connected_components_features:
        cycle_nodes_of_the_component1, component_nodes1 = connected_components_features[i]
        best_align_length = 0
        best_align_rmsd = 100.0
        best_align__j = 0
        for j in connected_components_features:
            if i != j:
                cycle_nodes_of_the_component2, component_nodes2 = connected_components_features[j]
                # Find the best edge from any node in the component1 to any node in the cycle of component 2 (best align from component2 to component 1)
                (source_loop, dest_loop, align_length, rmsd) = find_best_edge_between_components(cycle_nodes_of_the_component2, component_nodes1, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)
                component_adj_matrix[i][j] = (source_loop, dest_loop, align_length, rmsd)
                # Find the component that has the best incoming edge
                # if is_better_alignment_score((rmsd, align_length), (best_align_rmsd, best_align_length), align_len_threshold, is_normalized_score):
                #     best_align_rmsd = rmsd
                #     best_align_length = align_length
                #     best_align__j = j
        # Find the component that has the best outgoing edge
    #     if is_better_alignment_score((best_align_rmsd, best_align_length), (overall_best_align_rmsd, overall_best_align_length), align_len_threshold, is_normalized_score):
    #         overall_best_align_rmsd = best_align_rmsd
    #         overall_best_align_length = best_align_length
    #         best_edge_index_i = i

    # print best_edge_index_i

    best_traversal_rmsd = 1000000.0
    best_component_ordering = None
    # Find the traversal that reduces the overall traversal rmsd
    if start_node_list == None:
        start_node_list = connected_components_features.keys()

    for i in start_node_list:
        # Order components with traversal order of Prim's MST  algorithm
        component_ordering = get_component_ordering(component_adj_matrix, i, align_len_threshold, is_normalized_score)
        sum = 0.0
        for _,_,_,rmsd,_ in component_ordering:
            sum += rmsd
        
        if sum < best_traversal_rmsd:
            best_traversal_rmsd = sum
            best_component_ordering = copy.deepcopy(component_ordering)
    # print('Best traversal RMSD ' + str(best_traversal_rmsd))
    # else:
    #     best_component_ordering = get_component_ordering(component_adj_matrix, start_node, align_len_threshold, is_normalized_score)
    return best_component_ordering

# def get_best_subfamily_and_component_ordering(merged_components_dict, connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score):
#     ### Find relational features among components
#     merged_components_adj_matrix = {}
#     for i in merged_components_dict:
#         merged_components_adj_matrix[i] = {}
#         for j in merged_components_dict:
#             merged_components_adj_matrix[i][j] = ((0,''),(0,''),0,0.0)

#     max_edge_index_i = 0
#     overall_best_align_length = 0
#     overall_best_align_rmsd = 100.0

#     merged_components_features = {}
#     for k in sorted(merged_components_dict):
#         merged_components_nodes = []
#         merged_cycle_nodes_of_the_components = []
#         for i in sorted(merged_components_dict[k]):
#             _, cycle_nodes_of_the_component, component_nodes, _ = connected_components_features[i]
#             merged_cycle_nodes_of_the_components += cycle_nodes_of_the_component
#             merged_components_nodes += component_nodes            

#         merged_components_features[k] = (merged_cycle_nodes_of_the_components, merged_components_nodes)

#     ## Build graph with pair-wise edge of the best loop-alignment choice among the merged components (subfamily)
#     for i in sorted(merged_components_dict):
#         merged_cycle_nodes_of_the_components1, merged_components_nodes1 = merged_components_features[i]
#         best_align_length = 0
#         best_align_rmsd = 100.0
#         max_j = 0
#         for j in sorted(merged_components_dict):
#             if i != j:
#                 merged_cycle_nodes_of_the_components2, merged_components_nodes2 = merged_components_features[j]
#                 # Find the best edge from any node in the component1 to any node in the cycle of component 2 (best align from component2 to component 1)
#                 (source_loop, dest_loop, align_length, rmsd) = find_best_edge_between_components(merged_cycle_nodes_of_the_components2, merged_components_nodes1, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)
#                 merged_components_adj_matrix[i][j] = (source_loop, dest_loop, align_length, rmsd)
#                 # Find the merged component that has the best incoming edge
#                 if is_better_alignment_score((rmsd, align_length), (best_align_rmsd, best_align_length), align_len_threshold, is_normalized_score):
#                     best_align_rmsd = rmsd
#                     best_align_length = align_length
#                     max_j = j
#         # Find the component that has the best outgoing edge
#         if is_better_alignment_score((best_align_rmsd, best_align_length), (overall_best_align_rmsd, overall_best_align_length), align_len_threshold, is_normalized_score):
#             overall_best_align_rmsd = best_align_rmsd
#             overall_best_align_length = best_align_length
#             max_edge_index_i = i

#     ## Build graph with pair-wise edge of the best loop-alignment choice among the components inside merged components
#     for k in sorted(merged_components_dict):
#         for k in sorted(merged_components_dict):
#             central_node1, cycle_nodes_of_the_component1, component_nodes1, component_directed_adjacency_list1 = connected_components_features[i]
#             best_align_length = 0
#             best_align_rmsd = 100.0
#             max_j = 0
#             for j in merged_components_dict[k]:
#                 if i != j:
#                     central_node2, cycle_nodes_of_the_component2, component_nodes2, component_directed_adjacency_list2 = connected_components_features[j]
#                     # Find the best edge from any node in the component1 to any node in the cycle of component 2 (best align from component2 to component 1)
#                     (source_loop, dest_loop, align_length, rmsd) = find_best_edge_between_components(cycle_nodes_of_the_component2, component_nodes1, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)
#                     component_adj_matrix[i][j] = (source_loop, dest_loop, align_length, rmsd)
#                     # Find the component that has the best incoming edge
#                     if is_better_alignment_score((rmsd, align_length), (best_align_rmsd, best_align_length), align_len_threshold, is_normalized_score):
#                         best_align_rmsd = rmsd
#                         best_align_length = align_length
#                         max_j = j
#             # Find the component that has the best outgoing edge
#             if is_better_alignment_score((best_align_rmsd, best_align_length), (overall_best_align_rmsd, overall_best_align_length), align_len_threshold, is_normalized_score):
#                 overall_best_align_rmsd = best_align_rmsd
#                 overall_best_align_length = best_align_length
#                 max_edge_index_i = i

#     best_traversal_rmsd = 1000000.0
#     best_component_ordering = None
#     for i in connected_components_features:
#         # Order components with traversal order of Prim's MST  algorithm
#         component_ordering = get_component_ordering(component_adj_matrix, i, align_len_threshold, is_normalized_score)
#         sum = 0.0
#         for _,_,_,rmsd,_ in component_ordering:
#             sum += rmsd
#         # print sum
#         if sum < best_traversal_rmsd:
#             best_traversal_rmsd = sum
#             best_component_ordering = copy.deepcopy(component_ordering)

#     print('Best traversal RMSD ' + str(best_traversal_rmsd))
#     return best_component_ordering

def dfs(graph, node, visited, parent, traverse_list_with_parent):
    if node not in visited:
        visited.append(node)
        traverse_list_with_parent.append((node, parent))
        for n in sorted(graph[node], key = lambda x:len(graph[x])):
            dfs(graph, n, visited, node, traverse_list_with_parent)

def bfs(graph, start, visited, parent, traverse_list_with_parent):
    queue = [(start, parent)]
    while queue:
        (node, parent) = queue.pop(0)
        if node not in visited:
            visited.append(node)
            traverse_list_with_parent.append((node, parent))            
            for n in sorted(graph[node], key = lambda x:len(graph[x])):
                queue.append((n, node))

def sort_value(edge, align_len_threshold, is_length_adjusted_score, max_threshold = 10000.0):
    rmsd, align_len = edge    
    # if is_length_adjusted_score:
        # return rmsd
    inverted_val_of_align_len = (max_threshold - align_len)
    if is_acceptable_align_len(align_len, align_len_threshold):
        return (rmsd * 1000 * max_threshold) + inverted_val_of_align_len
    else:
        # Returning a pseudo RMSD, considering 100000 as max possible RMSD
        # Subtracting align_len as lower RMSD is better, while larger align_len is better
        return inverted_val_of_align_len

def dijkstra(graph, weighted_directed_adjacency_matrix, start, parent, align_len_threshold, is_length_adjusted_score):
    visited = []
    traverse_list_with_parent = []
    traverse_list_with_parent_and_weight = []
    heap = []
    
    heapq.heappush(heap, (0.0, (0.0, 0), start, parent))

    while len(heap) > 0:
        (_, edge_weight, node, parent) = heapq.heappop(heap)
        visited.append(node)
        traverse_list_with_parent.append((node, parent))
        traverse_list_with_parent_and_weight.append((node, parent, edge_weight))
        for n in graph[node]:
            if n not in visited:
                edge_weight = graph[node][n]
                traversal_weight = None
                if traversal_algorithm == "dijkstra":
                    traversal_weight = sort_value(edge_weight, align_len_threshold, is_length_adjusted_score)
                elif traversal_algorithm == "root-oriented":
                    traversal_weight = sort_value(weighted_directed_adjacency_matrix[start][n], align_len_threshold, is_length_adjusted_score)
                
                heapq.heappush(heap, (traversal_weight, edge_weight, n, node))

    return traverse_list_with_parent, traverse_list_with_parent_and_weight

def acceptable_edge_for_merging(rmsd, align_len, align_len_threshold, is_length_adjusted_score):
    if is_acceptable_align_len(align_len, align_len_threshold):

        # need to adjust rmsd rank
        if is_length_adjusted_score:
            rmsd = rmsd * math.sqrt(align_length)

        if rmsd <= rmsd_threshold_for_merging:
            return True


        # rmsd_rank = get_rmsd_rank(rmsd, align_len, is_length_adjusted_score)
        # # if rmsd_rank <= 2:
        # if rmsd_rank < 2:
        # # if rmsd <= acceptable_RMSD_threshold_for_merging:
        #     return True

    return False

# def acceptable_edge_for_merging(rmsd, align_len, align_len_threshold, merging_rmsd_threshold, is_length_adjusted_score):
#     if is_acceptable_align_len(align_len, align_len_threshold):
#         if is_acceptable_rmsd(rmsd, align_len, is_length_adjusted_score, merging_rmsd_threshold):
#             return True

#     return False

def align_to_ordering(component_ordering, connected_components_features, cluster_pairwise_alignment_details):

    ### Order loops according to the best component first, and align all of them on it
    merged_component_features = []
    merged_component_features.append([])    
    edges_in_merged_components = []
    edges_in_merged_components.append([])
    loop_ordering = [] 
    loop_ordering.append([])

    edges_among_all_components = []

    if len(component_ordering) > 0:
        (i, _, _, _, _) = component_ordering[0][0]
        central_node, _, _, _ = connected_components_features[i]
        node_to_start = central_node
        loop_ordering[-1].append((node_to_start, None))
        i, r1 = node_to_start
        for component_list in component_ordering:
            for (k, _, _, _, _) in component_list:
                # central_node, _, component = connected_components_features[i]
                _, _, component_nodes, _ = connected_components_features[k]
                _, _, component_nodes, _ = connected_components_features[k]
                first_item = True
                for j, r2 in component_nodes:
                    if (j, r2) != node_to_start:
                        align_rmsd, align_len = get_rmsd_align_len(i, r1, j, r2, cluster_pairwise_alignment_details)                    
                        loop_ordering[-1].append(((j, r2), node_to_start))
                        if first_item == True:
                            edges_in_merged_components[-1].append((node_to_start, (j, r2), (align_rmsd, align_len)))
                        first_item = False

                merged_component_features[-1].append(connected_components_features[k])
    
    return loop_ordering, merged_component_features, edges_among_all_components, edges_in_merged_components

def get_align_to_rmsd_info(cluster_id, cluster_pairwise_alignment_details, alignment_data, fp_align_len_threshold = None):
    global scanx_align_to_superimposition
    prev_status = scanx_align_to_superimposition
    scanx_align_to_superimposition = True
    ordered_dependency_list, merged_components_features, edges_among_all_components, edges_in_merged_components = generate_loop_print_dependency_v2(cluster_id, cluster_pairwise_alignment_details, alignment_data, fp_align_len_threshold)
    avg_rmsd, total_alignment_length = get_family_rmsd_and_alignment_summary(ordered_dependency_list, cluster_pairwise_alignment_details)

    scanx_align_to_superimposition = prev_status
    return avg_rmsd, total_alignment_length

def generate_loop_print_dependency_v2(cluster_id, cluster_pairwise_alignment_details, alignment_data, fp_align_len_threshold = None):
    align_len_threshold = generate_align_length_threshold(cluster_pairwise_alignment_details)

    if fp_align_len_threshold != None:
        fp_align_len_threshold.write(cluster_id + ',' + str(align_len_threshold) + '\n')

    adjacency_list, directed_adjacency_list, transpose_directed_adjacency_list = create_best_alignment_graph(cluster_pairwise_alignment_details, align_len_threshold)
    weighted_directed_adjacency_matrix = create_weighted_complete_directed_graph(cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)

    connected_components = get_connected_components(adjacency_list)
    # print('\nConnected Components')
    # print('Component count = ' + str(len(connected_components)))

    # Find features (central_node, cycle_nodes_of_the_component, component_nodes, component_directed_adjacency_list) of components
    connected_components_features = get_component_features(cluster_pairwise_alignment_details, directed_adjacency_list, transpose_directed_adjacency_list, connected_components)
 
    # Write some analysis on the length and rmsd of components to variation data log file
    # analyze_component_features(cluster_id, alignment_data, connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, align_dir, log_file_list, is_length_adjusted_score)
    
    component_ordering = merge_and_order_connected_components(connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)

    # component_ordering = get_best_component_ordering(connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)

    if scanx_align_to_superimposition:
        return align_to_ordering(component_ordering, connected_components_features, cluster_pairwise_alignment_details)

    ### Order loops according to the component first, and then the weighted BFS. Merge components accordingly.
    merged_component_features = []
    edges_among_all_components = []
    edges_in_merged_components = []
    loop_ordering = []
    first_component = True
    # rmsd_list = []
    # rmsd_list1 = []
    # print_a_list(component_ordering)
    for component_list in component_ordering:
        loop_ordering.append([])
        merged_component_features.append([])
        edges_in_merged_components.append([])
        first_item = True
        for (i, node_to_start, parent, align_rmsd, align_len) in component_list:
            
            if not first_item:
                edges_in_merged_components[-1].append((parent, node_to_start, (align_rmsd, align_len)))
            first_item = False

            # central_node, _, component = connected_components_features[i]
            central_node, _, component_nodes, component_directed_adjacency_list = connected_components_features[i]

            if first_component:
                node_to_start = central_node
            else:
                edges_among_all_components.append((parent, node_to_start, (align_rmsd, align_len)))
                # rmsd_list.append(align_rmsd)
                # rmsd_list1.append((align_rmsd, node_to_start, parent))

            # dfs(component_directed_adjacency_list, node_to_start, visited, parent, traverse_list_with_parent)
            # traverse_list_with_parent, traverse_list_with_parent_and_weight = bfs_weighted(component_directed_adjacency_list, node_to_start, parent, align_len_threshold, is_length_adjusted_score)
            # traverse_list_with_parent = Prims_MST(component_directed_adjacency_list, node_to_start, parent, align_len_threshold, is_length_adjusted_score)
            traverse_list_with_parent, traverse_list_with_parent_and_weight = dijkstra(component_directed_adjacency_list, weighted_directed_adjacency_matrix, node_to_start, parent, align_len_threshold, is_length_adjusted_score)

            # for node, parent, weight in traverse_list_with_parent_and_weight:
            #     # print node, parent, weight
            #     if parent != None:
            #         rmsd_list.append(weight[0])
            #         rmsd_list1.append((weight[0], node, parent))

            merged_component_features[-1].append(connected_components_features[i])
            loop_ordering[-1] += traverse_list_with_parent
            first_component = False 
    # generate_subfamily_representative(cluster_id, loop_ordering, alignment_data, cluster_pairwise_alignment_details, align_len_threshold)
    return loop_ordering, merged_component_features, edges_among_all_components, edges_in_merged_components

# def generate_loop_print_dependency(cluster_id, cluster_pairwise_alignment_details, alignment_data, fp_align_len_threshold = None):
#     print(cluster_id)
#     align_len_threshold = generate_align_length_threshold(cluster_pairwise_alignment_details)

#     # merging_rmsd_threshold = generate_merging_rmsd_threshold(cluster_pairwise_alignment_details)
#     # fp_rmsd_threshold = open("Familywise_RMSD_Threshold_" + loop_type + ".txt", 'a')
#     # fp_rmsd_threshold.write(cluster_id + ',' + str(merging_rmsd_threshold) + '\n')
#     # fp_rmsd_threshold.close()

#     if fp_align_len_threshold != None:
#         fp_align_len_threshold.write(cluster_id + ',' + str(align_len_threshold) + '\n')

#     adjacency_list, directed_adjacency_list, transpose_directed_adjacency_list = create_best_alignment_graph(cluster_pairwise_alignment_details, align_len_threshold)
#     weighted_directed_adjacency_matrix = create_weighted_complete_directed_graph(cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)

#     connected_components = get_connected_components(adjacency_list)
#     print('\nConnected Components')
#     print('Component count = ' + str(len(connected_components)))

#     # Find features (central_node, cycle_nodes_of_the_component, component_nodes, component_directed_adjacency_list) of components
#     connected_components_features = get_component_features(cluster_pairwise_alignment_details, directed_adjacency_list, transpose_directed_adjacency_list, connected_components)
 
#     # Write some analysis on the length and rmsd of components to variation data log file
#     # analyze_component_features(cluster_id, alignment_data, connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, align_dir, log_file_list, is_length_adjusted_score)

#     filtered_connected_components_features = {}
#     for i in connected_components_features:
#         _, cycle_nodes_of_the_component, component_nodes, _ = connected_components_features[i]
#         filtered_connected_components_features[i] = cycle_nodes_of_the_component, component_nodes

#     component_ordering = get_best_component_ordering(filtered_connected_components_features, cluster_pairwise_alignment_details, align_len_threshold, is_normalized_score)

#     # if scanx_align_to_superimposition:
#     #     return align_to_ordering(component_ordering, connected_components_features, cluster_pairwise_alignment_details)

#     ### Order loops according to the component first, and then the weighted BFS. Merge components accordingly.
#     merged_component_features = []
#     edges_among_all_components = []
#     edges_in_merged_components = []
#     loop_ordering = []
#     first_component = True
#     rmsd_list = []
#     rmsd_list1 = []
#     for (i, node_to_start, parent, align_rmsd, align_len) in component_ordering:
#         # central_node, _, component = connected_components_features[i]
#         central_node, _, component_nodes, component_directed_adjacency_list = connected_components_features[i]

#         if first_component:
#             node_to_start = central_node
#         else:
#             edges_among_all_components.append((parent, node_to_start, (align_rmsd, align_len)))
#             rmsd_list.append(align_rmsd)
#             rmsd_list1.append((align_rmsd, node_to_start, parent))

#         # dfs(component_directed_adjacency_list, node_to_start, visited, parent, traverse_list_with_parent)
#         # traverse_list_with_parent, traverse_list_with_parent_and_weight = bfs_weighted(component_directed_adjacency_list, node_to_start, parent, align_len_threshold, is_length_adjusted_score)
#         # traverse_list_with_parent = Prims_MST(component_directed_adjacency_list, node_to_start, parent, align_len_threshold, is_length_adjusted_score)
#         traverse_list_with_parent, traverse_list_with_parent_and_weight = dijkstra(component_directed_adjacency_list, weighted_directed_adjacency_matrix, node_to_start, parent, align_len_threshold, is_length_adjusted_score)

#         for node, parent, weight in traverse_list_with_parent_and_weight:
#             # print node, parent, weight
#             if parent != None:
#                 rmsd_list.append(weight[0])
#                 rmsd_list1.append((weight[0], node, parent))

#         if first_component or not merge_components or not acceptable_edge_for_merging(align_rmsd, align_len, align_len_threshold, is_length_adjusted_score):
#         # if first_component or not merge_components or not acceptable_edge_for_merging(align_rmsd, align_len, align_len_threshold, merging_rmsd_threshold, is_length_adjusted_score):
#             loop_ordering.append([])
#             merged_component_features.append([])
#             edges_in_merged_components.append([])
#             print('Component ' + str(i) + ' Not Merged')
#         else:
#             edges_in_merged_components[-1].append((parent, node_to_start, (align_rmsd, align_len)))
#             print('Component ' + str(i) + ' merged')
#             # parent = None
#             # node_to_start = central_node

#         # print(align_rmsd)
#         # print(align_len)

#         merged_component_features[-1].append(connected_components_features[i])
#         loop_ordering[-1] += traverse_list_with_parent
#         first_component = False 

#     # print('RMSD Z-Scores')
#     # print(cluster_id)
#     # z_scores = list(map(lambda x: (round(x[0],1),round(x[1],1)), get_z_scores(rmsd_list)))
#     # rmsd_with_loop = list(map(lambda x: (round(x[0],1),x[1][1],x[2][1]), rmsd_list1))
#     # if cluster_id == 'kink-turn':
#     # if cluster_id == 'sarcin-ricin':
#     # if cluster_id == 'reverse-kink-turn':
#         # sys.exit()

#     generate_subfamily_representative(cluster_id, loop_ordering, alignment_data, cluster_pairwise_alignment_details, align_len_threshold)

#     return loop_ordering, merged_component_features, edges_among_all_components, edges_in_merged_components

def get_component_organism_stat(ordered_dependency_list, pdb_organism_details):   #stat by org_type only
    component_organism_stat = {}
    for component_id, component in enumerate(ordered_dependency_list):
        component_organism_stat[component_id] = {}
        for component_loop_id, ((i, loop), parent) in enumerate(component):
            pdb_chain = loop.strip().split(':')[0]
            if pdb_chain in pdb_organism_details:
                RNA_Types, organism, org_class, org_type, pdb_source = pdb_organism_details[pdb_chain]
            else:
                org_type = 'Unknown'
            if org_type not in component_organism_stat[component_id]:
                component_organism_stat[component_id][org_type] = 0
            component_organism_stat[component_id][org_type] += 1

    return component_organism_stat

def write_alignments_to_file(i, r1, j, r2, aln1, aln2, rmsd, parent_load_id, load_id, f):
    ###########################################################
    f.write(str(parent_load_id).ljust(25) + '\t\t')
    if input_index_type == 'fasta':
        f.write(r2.rjust(30) + '(FASTA)\t\t')
    r2_pdb_ind = convert_a_loop_from_FASTA_to_PDB(r2)
    f.write(r2_pdb_ind.rjust(30) + '(PDB)\t\t' + aln2.rjust(40) + '\t\t' + str(round(rmsd, 2)).rjust(15) + '\n')

    f.write(str(load_id).ljust(25) + '\t\t')
    if input_index_type == 'fasta':
        f.write(r1.rjust(30) + '(FASTA)\t\t')
    r1_pdb_ind = convert_a_loop_from_FASTA_to_PDB(r1)
    f.write(r1_pdb_ind.rjust(30) + '(PDB)\t\t' + aln1.rjust(40) + '\n')

    f.write('\n')
    ###########################################################

def get_boundary_canonical_interactions(boundary_index_list, pdb_index_list, index_residue_dict):
    boundary_interaction_list = []
    temp_s = None
    prev_e = None
    for i, (s, e) in enumerate(boundary_index_list):
        if i == 0:
            temp_s = s
            prev_e = e
            continue

        a = prev_e
        b = s
        bp = a + "-" + b
        line = bp + "," + index_residue_dict[a] + "-" + index_residue_dict[b] + ",W/W" + ",cis" + "\n"
        boundary_interaction_list.append(line)

        prev_e = e

    a = temp_s
    b = prev_e
    bp = a + "-" + b
    line = bp + "," + index_residue_dict[a] + "-" + index_residue_dict[b] + ",W/W" + ",cis" + "\n"
    boundary_interaction_list.append(line)

    return boundary_interaction_list

def write_alignments_with_interactions_to_file(load_id, r, aln, rmsd, pdb_res_map_dict, f, categorize_annotation = True):

    r_pdb_ind = convert_a_loop_from_FASTA_to_PDB(r)

    if(len(load_id) > 0):
        f.write(str(load_id).ljust(25) + '\t\t')
        if input_index_type == 'fasta':
            f.write(r + '(FASTA)')
        
        f.write(r_pdb_ind.rjust(50) + '(PDB)')
    
    else:
        if input_index_type == 'fasta':
            f.write(r + ' (FASTA)\n')
        
        f.write(r_pdb_ind + ' (PDB)')

    if len(aln) > 0:
        f.write('\t\t' + aln.rjust(40))
        f.write('\t\t' + str(round(rmsd, 2)).rjust(8))
    f.write('\n')

    # f4.write(str(parent_load_id).ljust(25) + '\t\t' + r2.rjust(30) + '\t\t')
    # r2_pdb_ind = convert_a_loop_from_FASTA_to_PDB(r2)
    # f4.write(r2_pdb_ind.rjust(25) + '\t\t' + aln2.rjust(40) + '\t\t' + str(round(rmsd, 2)).rjust(15) + '\n')

    f_loop = open(os.path.join(loop_dir, r + ".smf"))
    loop_info_lines = f_loop.readlines()
    f_loop.close()

    pdb_chain, loop_regions = r.strip().split(":")

    if pdb_chain not in pdb_res_map_dict:
        pdb_res_map_dict[pdb_chain] = load_pdb_res_map(pdb_chain)

    pdb_res_map = pdb_res_map_dict[pdb_chain]
    loop_regions = loop_regions.strip().split("_")
    loop_regions_seq = "".join(loop_info_lines[1].strip().split("..."))

    loop_regions_pdb_ind = []
    # index_residue_dict = {}
    pdb_index_list = []
    boundary_index_list = []
    for part, region in enumerate(loop_regions):
        region_pdb_index = []
        
        s, e = region.strip().split("-")
        s = int(s)
        e = int(e)
        for fasta_indx in range(s, e+1):
            if len(pdb_res_map[fasta_indx][2]) > 0:
                pdb_index_list.append(str(pdb_res_map[fasta_indx][1]) + '.' + str(pdb_res_map[fasta_indx][2]))
            else:
                pdb_index_list.append(pdb_res_map[fasta_indx][1])

        start_ind = str(pdb_res_map[s][1])
        if len(pdb_res_map[s][2]) > 0:
            start_ind += '.' + str(pdb_res_map[s][2])

        end_ind = pdb_res_map[e][1]
        if len(pdb_res_map[e][2]) > 0:
            end_ind += '.' + str(pdb_res_map[e][2])
            
        loop_regions_pdb_ind.append(start_ind + "-" + end_ind)
        boundary_index_list.append((start_ind, end_ind))

    index_residue_dict = dict(zip(pdb_index_list, loop_regions_seq))

    cat1_lines = []
    cat2_lines = []
    other_lines = []
    stack_lines = []
    in_stack = False
    for line in loop_info_lines:
        if line.startswith(">"):
            continue
        elif line.startswith("#info=stacking") or in_stack == True:
            # break
            pieces = line.strip().split(",")
            if len(pieces) == 2:
                a, b = pieces[0].strip().split("-")
                a = pdb_index_list[int(a)]
                b = pdb_index_list[int(b)]
                bp = a + "-" + b
                line = bp + "," + index_residue_dict[a] + "-" + index_residue_dict[b] + "," + pieces[1] + '\n'

            stack_lines.append(line)
            in_stack = True
        else:
            pieces = line.strip().split(",")
            if len(pieces) == 3:
                a, b = pieces[0].strip().split("-")
                a = pdb_index_list[int(a)]
                b = pdb_index_list[int(b)]
                bp = a + "-" + b
                line = bp + "," + index_residue_dict[a] + "-" + index_residue_dict[b] + "," + pieces[1] + "," + pieces[2] + '\n'

                if categorize_annotation == True and (pieces[1].strip() == 'S/H' or pieces[1].strip() == 'H/S' or pieces[1].strip() == 'S/S'):
                    cat1_lines.append(line)
                else:
                    cat2_lines.append(line)
                #     f4.write(line)
            else:
                other_lines.append(line)
            #     f4.write(line)
            
            # f.write(line)
    if len(other_lines) > 0:
        f.write(''.join(other_lines))
    if len(cat1_lines) > 0:
        f.write('\n')
        f.write(''.join(cat1_lines))
        f.write('\n')
    boundary_interaction_list = get_boundary_canonical_interactions(boundary_index_list, pdb_index_list, index_residue_dict)
    f.write(boundary_interaction_list[-1])
    if len(cat2_lines) > 0:
        f.write(''.join(cat2_lines))
    f.write(''.join(boundary_interaction_list[:-1]))
    if len(stack_lines) > 0:
        # f.write('\n')
        f.write(''.join(stack_lines))

def extract_rmsd_from_dict(rmsd_data_list_dict, i, r1, j, r2):
    _, rmsd_list = rmsd_data_list_dict[(i, r1)]
    rmsd = -1.0
    for (x, rx, rmsdx, _) in rmsd_list:
        if j == x and rx == r2:
            rmsd = rmsdx
            break
    return rmsd

def get_align_len(aln1, aln2):
    len1 = len(aln1.replace('-',''))
    len2 = len(aln2.replace('-',''))
    return len1 if len1 < len2 else len2;

# def is_better_rmsd_len(rmsd_len1, rmsd_len2, align_loop_count):
#     rmsd1, len1 = rmsd_len1
#     rmsd2, len2 = rmsd_len2

#     len1 /= align_loop_count
#     len2 /= align_loop_count

#     if abs(rmsd1 - rmsd2) < 0.000000001:
#         return len1 > len2
#     return rmsd1 < rmsd2

def is_better_rmsd_len(rmsd_len1, rmsd_len2, align_loop_count):
    rmsd1, len1 = rmsd_len1
    rmsd2, len2 = rmsd_len2
    len1 /= align_loop_count
    len2 /= align_loop_count

    if abs(len1 - len2) < 0.000000001:
        return rmsd1 < rmsd2
    return len1 > len2

def generate_subfamily_representative(fp_representative, cluster_id, component_id, component, alignment_data, rmsd_data_list_dict, pdb_res_map_dict, align_len_threshold):

    # if is_generate_subfamily_representative == False:
        # return

    # output_dir = os.path.join(superimposition_output_dir, 'componentwise_analysis')
    # output_dir = representative_dir
    
    # create_directory(output_dir)

    # f = open(os.path.join(output_dir, str(cluster_id) + "_representatives.txt"), "w")
    # pdb_res_map_dict = {}
    # for component_id, component in enumerate(ordered_dependency_list):

    representative_loop = ''
    represetative_i = -1
    max_acceptable_align_count = 0
    best_weighted_avg_rmsd = (20000.0, 0)
    loop_count = len(component) - 1
    for component_loop_id1, ((i, r1), _) in enumerate(component):
        rmsd_align_len_list = []
        
        for component_loop_id2, ((j, r2), _) in enumerate(component):
            if component_loop_id1 != component_loop_id2:
                (t1, t2, zscore, cr1, cr2, aln2, aln1, score) = alignment_data[cluster_id][strToNode(r2)][strToNode(r1)]
                rmsd = extract_rmsd_from_dict(rmsd_data_list_dict, i, r1, j, r2)
                align_len = get_align_len(aln1, aln2);

                if is_acceptable_align_len(align_len, align_len_threshold):
                    rmsd_align_len_list.append((rmsd,align_len))

        weighted_avg_rmsd = get_weighted_avg_rmsd(rmsd_align_len_list)
        acceptable_align_count = len(rmsd_align_len_list)

        # Try to make sure the represetative has acceptable alignment with at least half of the members
        if ((max_acceptable_align_count < int((loop_count + 1) / 2) and 
            (acceptable_align_count > max_acceptable_align_count or 
            (acceptable_align_count == max_acceptable_align_count and is_better_rmsd_len(weighted_avg_rmsd , best_weighted_avg_rmsd, loop_count))))
            or
            (acceptable_align_count >= int((loop_count + 1) / 2) and is_better_rmsd_len(weighted_avg_rmsd , best_weighted_avg_rmsd, loop_count))):
            max_acceptable_align_count = acceptable_align_count
            best_weighted_avg_rmsd = weighted_avg_rmsd
            representative_loop = r1
            represetative_i = i

    avg_rmsd, total_align_len = best_weighted_avg_rmsd
    if cluster_id in known_motif_fullname:
        fp_representative.write(known_motif_fullname[cluster_id] + '-Sub' + str(component_id + 1) + ':\n')
    else:
        fp_representative.write(cluster_id + '-Sub' + str(component_id + 1) + ':\n')
    if output_env == 'local':
        fp_representative.write("Acceptable Align Count = " + str(max_acceptable_align_count) + "/" + str(loop_count) + ",\nWeighted Average RMSD = (" + str(round(avg_rmsd, 2)) + ',' + str(total_align_len) + ")\n")

    write_alignments_with_interactions_to_file("", representative_loop, "", best_weighted_avg_rmsd, pdb_res_map_dict, fp_representative, False)
    fp_representative.write("\n\n")

    return represetative_i, representative_loop

def generate_representative_loop_image(time_in_distance_calc, representative_dir, rotation_version, cluster_id, component_id, i, r, loop_display_info_dict, draw_figures, show_extended_loop, show_label):
    if draw_figures == False:
        return 0

    image_fname = os.path.join(representative_dir, add_rotation_version_prefix(rotation_version) + cluster_id + '-Sub' + str(component_id + 1) + '_repr.png')

    display_load_name, align_load_name, chain_load_name = loop_display_info_dict[(i,r)]
    # print(display_load_name, align_load_name, chain_load_name)
    # pymol.cmd._do('hide all')
    # pymol.cmd.sync()
    # wait_for_certain_time_according_to_wait_factor()
    # time.sleep(.200)
    
    display_color = 'gray'
    cano_atom_color = 'orange'
    other_atom_color = 'blue'
    bp_atom_color = 'green'
    bp_line_color = 'red'

    cano_seqnums = []
    other_seqnums = []
    bp_seqnums = []
    bp_list = []
    
    pdb_chain, loop_regions = r.strip().split(":")
    pdb_res_map = load_pdb_res_map(pdb_chain)
    loop_regions = loop_regions.strip().split('_')
    pdb_index_list = []
    for region in loop_regions:
        s, e = region.strip().split("-")
        s = int(s)
        e = int(e)
        cano_seqnums.append(pdb_res_map[s][1] + pdb_res_map[s][2])
        cano_seqnums.append(pdb_res_map[e][1] + pdb_res_map[e][2])
        for fasta_indx in range(s+1, e):
            other_seqnums.append(pdb_res_map[fasta_indx][1] + pdb_res_map[fasta_indx][2])
        for fasta_indx in range(s, e+1):
            if fasta_indx != s and fasta_indx != e:
                other_seqnums.append(pdb_res_map[fasta_indx][1] + pdb_res_map[fasta_indx][2])
            pdb_index_list.append((pdb_res_map[fasta_indx][1], pdb_res_map[fasta_indx][2]))

    f_loop = open(os.path.join(loop_dir, r + ".smf"))
    loop_info_lines = f_loop.readlines()
    f_loop.close()

    for line in loop_info_lines:
        if line.startswith(">"):
            continue
        elif line.startswith("#info=stacking"):
            break
        else:
            pieces = line.strip().split(",")
            if len(pieces) == 3:
                a, b = pieces[0].strip().split("-")
                seqnum1, icode1 = pdb_index_list[int(a)]
                seqnum2, icode2 = pdb_index_list[int(b)]
                a = seqnum1 + icode1
                b = seqnum2 + icode2
                bp_seqnums.append(a)
                bp_seqnums.append(b)
                bp_list.append((a, b))
    
    other_bp_load_name = cluster_id + '_sub' + str(component_id + 1) + '_other_bp'
    bp_atoms_load_name = cluster_id + '_sub' + str(component_id + 1) + '_bp_atoms'
    cano_bp_load_name = cluster_id + '_sub' + str(component_id + 1) + '_cano_bp'
    dist_load_name = ''
    pymol.cmd.select(other_bp_load_name, chain_load_name + ' and (%s)' % ' or '.join(list(map(lambda x: 'resi '+x, other_seqnums))))
    pymol.cmd.select(bp_atoms_load_name, chain_load_name + ' and (%s)' % ' or '.join(list(map(lambda x: 'resi '+x, bp_seqnums))))
    pymol.cmd.select(cano_bp_load_name, chain_load_name + ' and (%s)' % ' or '.join(list(map(lambda x: 'resi '+x, cano_seqnums))))

    config_pymol_cartoon(display_color, True)
    # pymol.cmd._do('set stick_transparency, 0.5')

    if show_extended_loop:
        pymol.cmd.show('cartoon', chain_load_name)
    else:
        pymol.cmd.show('cartoon', display_load_name)

    # pymol.cmd.show('stick', bp_atoms_load_name)

    pymol.cmd.color(display_color, chain_load_name)
    pymol.cmd.color(other_atom_color, other_bp_load_name)
    pymol.cmd.color(bp_atom_color, bp_atoms_load_name)
    pymol.cmd.color(cano_atom_color, cano_bp_load_name)

    a_load_name = 'a_atoms'
    b_load_name = 'b_atoms'
    # print(bp_list)
    dist_load_names = []
    atom_list_for_distance_measure = ["N1", "N3", "N7", "C2", "C4", "C5" "C6", "C8"]
    for a, b in bp_list:
        # print(a,b)
        a_name_list = []
        b_name_list = []
        pymol.cmd.select(a_load_name, chain_load_name + ' and resi ' + a)
        pymol.cmd.select(b_load_name, chain_load_name + ' and resi ' + b)

        name_dict = { 'a_name_list' : [], 'b_name_list' : [] }
        pymol.cmd.iterate(a_load_name, 'a_name_list.append(name)', space=name_dict)
        pymol.cmd.iterate(b_load_name, 'b_name_list.append(name)', space=name_dict)
        # print(name_dict)
        if len(name_dict['a_name_list']) == 0 or len(name_dict['b_name_list']) == 0:
            logger.warning('Skipping min_distance calculation. Check.')
            logger.warning(r)
            logger.warning(str(a) + ', ' + str(b))
            continue
        min_distance = 1000
        min_a_name = ''
        min_b_name = ''
        a_single_name = 'a_atom'
        b_single_name = 'b_atom'
        time_ss = time.time()
        for a_name in name_dict['a_name_list']:
            # if not (a_name.startswith('C') or a_name.startswith('N')):
            if a_name not in atom_list_for_distance_measure:
                continue
            for b_name in name_dict['b_name_list']:
                # if not (b_name.startswith('C') or b_name.startswith('N')):
                if b_name not in atom_list_for_distance_measure:
                    continue
                pymol.cmd.select(a_single_name, a_load_name + ' and name ' + a_name)
                pymol.cmd.select(b_single_name, b_load_name + ' and name ' + b_name)
                distance = pymol.cmd.distance('dist', a_single_name, b_single_name)
                if distance < min_distance:
                    min_distance = distance
                    min_a_name = a_name
                    min_b_name = b_name
                pymol.cmd.delete('dist')
        # print(min_distance, min_a_name, min_b_name)
        time_dd = time.time() - time_ss
        time_in_distance_calc += time_dd

        pymol.cmd.select(a_single_name, a_load_name + ' and name ' + min_a_name)
        pymol.cmd.select(b_single_name, b_load_name + ' and name ' + min_b_name)
        dist_load_name = cluster_id + '_sub' + str(component_id + 1) + 'dist_' + min_a_name.replace("'", "p") + '_' + min_b_name.replace("'", "p")
        pymol.cmd.distance(dist_load_name, a_single_name, b_single_name)
        # pymol.cmd._do('hide labels, ' + dist_load_name)
        pymol.cmd.hide('labels', dist_load_name)
        # print(bp_line_color, dist_load_name)
        pymol.cmd.color(bp_line_color, dist_load_name)

        pymol.cmd.delete(a_single_name)
        pymol.cmd.delete(b_single_name)
        dist_load_names.append(dist_load_name)
    
    # pymol.cmd._do('zoom')
    pymol.cmd.zoom()
    wait_for_certain_time_according_to_wait_factor(len(atom_list_for_distance_measure))     #remove after changing distance code
    pymol.cmd.sync()
    # time.sleep(.100)
    pymol.cmd.png(image_fname, 1200, 1200, dpi=300, ray=1, quiet=1)
    # time.sleep(.100)
    wait_for_certain_files_to_be_generated([image_fname], False)
    pymol.cmd.sync()

    if save_pymol_session == True:
        pymol.cmd.deselect()
        session_fname = os.path.join(representative_dir, cluster_id + '-Sub' + str(component_id + 1) + '_repr.pse')
        pymol.cmd._do('save ' + os.path.join(representative_dir, session_fname))
        # time.sleep(.100)
        wait_for_certain_files_to_be_generated([session_fname], False)
        pymol.cmd.sync()

    # sys.exit()

    pymol.cmd.delete(a_load_name)
    pymol.cmd.delete(b_load_name)

    pymol.cmd.delete(other_bp_load_name)
    pymol.cmd.delete(bp_atoms_load_name)
    pymol.cmd.delete(cano_bp_load_name)
    for item in dist_load_names:
        pymol.cmd.delete(item)

    config_pymol_cartoon(display_color, show_label)
    pymol.cmd.hide()
    # print('deleting ', a_load_name)
    # print('deleting ', b_load_name)
    # print('deleting ', other_bp_load_name)
    # print('deleting ', bp_atoms_load_name)
    # print('deleting ', cano_bp_load_name)
    # print('deleting ', dist_load_names)

    wait_for_certain_time_according_to_wait_factor(1)
    pymol.cmd.sync()
    # time.sleep(.100)
    return time_in_distance_calc

def generate_table(summary_dir, subfamily_details_table_data, loop_type, is_latex=False):

    if is_latex == True:
        subfamily_tbl_fn = os.path.join(summary_dir, 'Subfamily_summary_table.txt')
    else:
        subfamily_tbl_fn = os.path.join(summary_dir, 'Subfamily_summary_table.tsv')

    fw = open(subfamily_tbl_fn, 'w')

    #  Table top
    if is_latex == True:
        fw.write('\\begin{table*}[b]\n')
        fw.write('\\tableparts{%\n')
        fw.write('\\caption{Subfamilies of the known motif families (Merging Threshold: RMSD  = ' + str(rmsd_threshold_for_merging) + ', Align Len Z-Score = ' + str(align_len_zscore_threshold) + ', Connectivity = ' + str(connectivity_test_threshold) + ', Max_align_len_in_equation? = ' + str(use_max_align_len_in_equation) + ')}\n')
        fw.write('\\label{table:subfamily}%\n')
        fw.write('}{%\n')
        fw.write('\\begin{tabular*}{0.83\\textwidth}{@{}llclccc@{}}\n')
        fw.write('\\toprule\n')
        fw.write('% \\cline{1-7}\n')

        # Table Title
        fw.write('\\multirow{2}{*}{Loop Type} & \\multirow{2}{*}{Motif Family} & \n')
        fw.write('\\multirow{2}{*}{\\begin{tabular}[c]{@{}c@{}} No of \\\\ Motifs\\end{tabular}} & \n')
        fw.write('\\multirow{2}{*}{\\begin{tabular}[l]{@{}l@{}} Subfamily \\\\ Count (Sizes)\\end{tabular}} & \n')
        fw.write('\\multicolumn{3}{c}{Superimposition Avg. RMSD / Aligned Length} \\\\ \n')
        fw.write('& & & & No Ordering & Ordered & Subfamilies \\\\\n')
        fw.write('\\colrule\n')
    else:
        fw.write('Subfamilies of the known motif families (Merging Threshold: RMSD  = ' + str(rmsd_threshold_for_merging) + ', Align Len Z-Score = ' + str(align_len_zscore_threshold) + ', Connectivity = ' + str(connectivity_test_threshold) + ')\n')
        fw.write('Loop Type \tMotif Family \tNo of Motifs \tSubfamily Count (Sizes) \tSuperimposition Avg. RMSD / Aligned Length\n')
        fw.write('\t\t\t\tNo Ordering \tOrdered \tSubfamilies\n')

    # Table Data
    if loop_type == 'IL':
        fw.write('Internal Loop (IL)')
    elif loop_type == 'HL':
        fw.write('Hairpin Loop (HL)')
    elif len(loop_type.strip()) > 0:
        fw.write(loop_type)

    # if is_latex == False:
    #     fw.write('\t')

    priority_order = ["GNRA", "GNAA", "GNGA", "Kink-turn", "reverse-Kink-turn", "Sarcin-ricin", "C-loop", "E-loop", "Hook-turn", "Tandem-shear",
     "Tetraloop-receptor", "L1-complex", "T-loop", "Rope-sling"]

    # Print families in priority order that has 2 or more subfamilies
    for cluster_id in priority_order:
        if cluster_id in subfamily_details_table_data:
            row = subfamily_details_table_data[cluster_id]
            # If there are more than one subfamily
            if ',' in row[2]:
                # if len(loop_type.strip()) > 0:
                if is_latex == True:
                    fw.write(' & ')
                else:
                    fw.write('\t')
                if is_latex == True:
                    fw.write(' & '.join(row) + ' \\\\\n')
                else:
                    fw.write('\t'.join(row) + '\n')

    # Print families NOT in priority order that has 2 or more subfamilies
    for cluster_id in subfamily_details_table_data:
        if cluster_id not in priority_order:
            row = subfamily_details_table_data[cluster_id]
            # If there are more than one subfamily
            if ',' in row[2]:
                # if len(loop_type.strip()) > 0:
                if is_latex == True:
                    fw.write(' & ')
                else:
                    fw.write('\t')
                if is_latex == True:
                    fw.write(' & '.join(row) + ' \\\\\n')
                else:
                    fw.write('\t'.join(row) + '\n')


    # Print families in priority order that has only 1 subfamily
    for cluster_id in priority_order:
        if cluster_id in subfamily_details_table_data:
            row = subfamily_details_table_data[cluster_id]
            # If there is exactly one subfamily
            if ',' not in row[2]:
                # if len(loop_type.strip()) > 0:
                if is_latex == True:
                    fw.write(' & ')
                else:
                    fw.write('\t')
                if is_latex == True:
                    fw.write(' & '.join(row) + ' \\\\\n')
                else:
                    fw.write('\t'.join(row) + '\n')

    # Print families NOT in priority order that has only 1 subfamily
    for cluster_id in subfamily_details_table_data:
        if cluster_id not in priority_order:
            row = subfamily_details_table_data[cluster_id]
            # If there is exactly one subfamily
            if ',' not in row[2]:
                # if len(loop_type.strip()) > 0:
                if is_latex == True:
                    fw.write(' & ')
                else:
                    fw.write('\t')
                if is_latex == True:
                    fw.write(' & '.join(row) + ' \\\\\n')
                else:
                    fw.write('\t'.join(row) + '\n')


    if is_latex == True:
        fw.write('& & & & & & \\\\\n')
    else:
        fw.write('\n')

    # fw.write('Internal Loop (IL) & Kink-Turn (KT) & 61 & 3 (52, 4, 5) & 2.269 & 0.690 & 0.447 \\\\\n')
    # fw.write('& reverse Kink-Turn (rKT) & 14 & 3 (6, 2, 6) & 1.060 & 0.900 & 0.629 \\\\\n')
    # fw.write('& Sarcin-Ricin (SR) & 72 & 4 (7, 55, 7, 3) & 1.724 & 0.625 & 0.419 \\\\\n')
    # fw.write('& C-Loop (CL) & 43 & 6 (9, 8, 11, 4, 9, 2) & 1.795 & 1.194 & 0.527 \\\\\n')
    # fw.write('& E-Loop (EL) & 49 & 3 (4, 43, 2) & 1.709 & 0.721 & 0.498 \\\\\n')
    # fw.write('& Tetraloop Receptor (TR) & 19 & 2 (17, 2) & 0.992 & 0.539 & 0.536 \\\\\n')
    # fw.write('& L1-Complex (L1C) & 6 & 2 (4, 2) & 2.400 & 1.701 & 0.833 \\\\\n')
    # fw.write('& Hook Turn (HT) & 33 & 3 (29, 2, 2) & 1.453 & 0.883 & 0.775 \\\\\n')
    # fw.write('& Tandem Sheared (TS) & 46 & 1 (46) & 1.016 & 0.495 & 0.495 \\\\\n')
    # fw.write('& T-Loop (TL) & 2 & 1 (2) & 0.500 & 0.500 & 0.500 \\\\\n')
    # fw.write('& Rope Sling (RS) & 9 & 1 (9) & 0.632 & 0.409 & 0.408 \\\\\n')
    # fw.write('& & & & & \\\\\n')
    # fw.write('Hairpin Loop (HL) & GNAA & 230 & 3 (226, 2, 2) & 0.774 & 0.465 & 0.269 \\\\\n')
    # fw.write('& GNGA & 64 & 3 (13, 41, 10) & 0.867 & 0.371 & 0.329 \\\\\n')
    # fw.write('& T-Loop (TL) & 109 & 7 (3, 86, 2, 2, 7, 4, 5) & 0.865 & 0.533 & 0.428 \\\\\n')

    if is_latex == True:
        # Table Ending
        fw.write('\\botrule\n')
        fw.write('\\end{tabular*}%\n')
        fw.write('}\n')
        fw.write('{}\n')
        fw.write('% {This is a table footnote}\n')
        fw.write('\\end{table*}\n')

    fw.close()

def generate_subfamily_details_table_data(cluster_id, ordered_dependency_list, avg_rmsd_align_to, total_alignment_length_align_to, avg_rmsd, total_alignment_length, avg_rmsd_sufamily_only, total_alignment_length_subfamily_only, subfamily_details_table_data):

    subfamily_instance_count = []
    total_loop_count = 0
    
    for component_id, component in enumerate(ordered_dependency_list):
        subfamily_instance_count.append(str(len(component)))
        total_loop_count += len(component)
    
    instance_count_string = ','.join(subfamily_instance_count)

    cluster_full_name = cluster_id
    cluster_shortcode = cluster_id
    if cluster_id in known_motif_fullname:
        cluster_full_name = known_motif_fullname[cluster_id]
    if cluster_id in known_motif_shortcode:
        cluster_shortcode = known_motif_shortcode[cluster_id]

    row_items = []
    cluster_name = cluster_full_name
    if cluster_full_name != cluster_shortcode:
        cluster_name += ' (' + cluster_shortcode + ')'

    avg_alignment_length_align_to = total_alignment_length_align_to / (total_loop_count - 1)
    avg_alignment_length = total_alignment_length / (total_loop_count - 1)
    avg_alignment_length_subfamily_only = total_alignment_length_subfamily_only / (total_loop_count - len(ordered_dependency_list))

    row_items.append(cluster_name)
    row_items.append(str(total_loop_count))
    row_items.append(str(len(ordered_dependency_list)) + ' (' + instance_count_string + ')')
    row_items.append("{:.3f}".format(round(avg_rmsd_align_to, 3)) + " / " + "{:.0f}".format(round(avg_alignment_length_align_to)))
    row_items.append("{:.3f}".format(round(avg_rmsd, 3)) + " / " + "{:.0f}".format(round(avg_alignment_length)))
    row_items.append("{:.3f}".format(round(avg_rmsd_sufamily_only,3)) + " / " + "{:.0f}".format(round(avg_alignment_length_subfamily_only)))
    subfamily_details_table_data[cluster_id] = row_items

def generate_componentwise_analysis_files(superimposition_output_dir, cluster_id, ordered_dependency_list, load_id_dict, alignment_data, rmsd_data):
    if generate_bp_ann_files == False:
        return

    output_dir = os.path.join(superimposition_output_dir, 'componentwise_analysis')
    
    create_directory(output_dir)

    if generate_loop_source_info == True:
        f_source1 = open(os.path.join(output_dir, str(cluster_id) + "_loop_source_data_summary.txt"), "w")
        f_source2 = open(os.path.join(output_dir, str(cluster_id) + "_loop_source_data_details.txt"), "w")

    source_dict = {}   

    _, rmsd_data_list_dict = rmsd_data[cluster_id]

    pdb_res_map_dict = {}
    for component_id, component in enumerate(ordered_dependency_list):

        f1 = open(os.path.join(output_dir, str(cluster_id) + "_" + str(component_id+1) + "_alignments_only.txt"), "w")
        f2 = open(os.path.join(output_dir, str(cluster_id) + "_" + str(component_id+1) + "_alignments_with_grouped_interactions.txt"), "w")
        f3 = open(os.path.join(output_dir, str(cluster_id) + "_" + str(component_id+1) + "_inter_subfam_alignments_only.txt"), "w")
        f4 = open(os.path.join(output_dir, str(cluster_id) + "_" + str(component_id+1) + "_alignments_with_interactions.txt"), "w")

        if generate_loop_source_info == True:
            source_dict['DeNovo'] = []
            source_dict['R3D'] = []
            source_dict['Both'] = []
            source_dict['N/A'] = []

        is_first = True
        for component_loop_id, ((i, r1), parent) in enumerate(component):

            if generate_loop_source_info == True:
                source_dict[get_loop_cluster_source(r1)].append(r1)

            if is_first == True:
                is_first = False
                # inter-subfamily edge
                if parent != None:
                    (j, r2) = parent
                    parent_load_id, parent_loop_id = load_id_dict[(j, r2)]
                    load_id, loop_id = load_id_dict[(i, r1)]
                    (t1, t2, zscore, cr1, cr2, aln2, aln1, score) = alignment_data[cluster_id][strToNode(r2)][strToNode(r1)]
                    rmsd = extract_rmsd_from_dict(rmsd_data_list_dict, i, r1, j, r2)
                    write_alignments_to_file(i, r1, j, r2, aln1, aln2, rmsd, parent_load_id, load_id, f3)
            else:
            # if parent != None:
                (j, r2) = parent
                parent_load_id, parent_loop_id = load_id_dict[(j, r2)]
                load_id, loop_id = load_id_dict[(i, r1)]
                (t1, t2, zscore, cr1, cr2, aln2, aln1, score) = alignment_data[cluster_id][strToNode(r2)][strToNode(r1)]
                rmsd = extract_rmsd_from_dict(rmsd_data_list_dict, i, r1, j, r2)
                write_alignments_to_file(i, r1, j, r2, aln1, aln2, rmsd, parent_load_id, load_id, f1)

                ###########################################################
                #*************************r2******************************#
                write_alignments_with_interactions_to_file(parent_load_id, r2, aln2, rmsd, pdb_res_map_dict, f2)
                write_alignments_with_interactions_to_file(parent_load_id, r2, aln2, rmsd, pdb_res_map_dict, f4, False)

                #*************************r2******************************#
                f2.write('\n')
                f4.write('\n')
                write_alignments_with_interactions_to_file(load_id, r1, aln1, rmsd, pdb_res_map_dict, f2)
                write_alignments_with_interactions_to_file(load_id, r1, aln1, rmsd, pdb_res_map_dict, f4, False)

                #*************************r1******************************#
                f2.write('\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                f2.write('\n\n\n')
                f4.write('\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                f4.write('\n\n\n')
                ###########################################################
  
        f1.close()
        f2.close()
        f4.close()
        f3.close()

        if generate_loop_source_info == True:
            f_source1.write(str(cluster_id) + '_' + str(component_id+1) + ' loop source summary: \n')
            if len(source_dict['DeNovo']) > 0:
                f_source1.write('From DeNovo: ' + str(len(source_dict['DeNovo'])) + '\n')
            if len(source_dict['R3D']) > 0:
                f_source1.write('From R3D: ' + str(len(source_dict['R3D'])) + '\n')
            if len(source_dict['Both']) > 0:
                f_source1.write('From Both: ' + str(len(source_dict['Both'])) + '\n')
            if len(source_dict['N/A']) > 0:
                f_source1.write('N/A: ' + str(len(source_dict['N/A'])) + '\n')
            f_source1.write('Total: ' + str(len(component)) + '\n')
            f_source1.write('\n\n')

            f_source2.write(str(cluster_id) + '_' + str(component_id+1) + ' loop source details: (FASTA index)\n')
            if len(source_dict['DeNovo']) > 0:
                f_source2.write('From DeNovo: (' + str(len(source_dict['DeNovo'])) + ')\n' + ', '.join(source_dict['DeNovo']) + '\n')
            if len(source_dict['R3D']) > 0:
                f_source2.write('From R3D: (' + str(len(source_dict['R3D'])) + ')\n' + ', '.join(source_dict['R3D']) + '\n')
            if len(source_dict['Both']) > 0:
                f_source2.write('From Both: (' + str(len(source_dict['Both'])) + ')\n' + ', '.join(source_dict['Both']) + '\n')
            if len(source_dict['N/A']) > 0:
                f_source2.write('N/A: (' + str(len(source_dict['N/A'])) + ')\n' + ', '.join(source_dict['N/A']) + '\n')
            f_source2.write('Total: ' + str(len(component)) + '\n')
            f_source2.write('\n\n')

    if generate_loop_source_info == True:
        f_source1.close()
        f_source2.close()

# write_alignments_with_interactions_to_file("", representative_loop, "", best_weighted_avg_rmsd, pdb_res_map_dict, f, False)

def generate_componentwise_bp_annotation_files(subfamily_details_dir, cluster_id, ordered_dependency_list, load_id_dict):
    if generate_bp_ann_files == False:
        return
        # return ordered_dependency_list

    # mapping_dir = pdb_fasta_mapping_dir
    # loop_dir = loop_dir
    
    create_directory(subfamily_details_dir)

    ordered_dependency_list_with_bp_ann = []

    pdb_res_map_dict = {}
    fw_bp_details = open(os.path.join(subfamily_details_dir, str(cluster_id) + "_bp_details.txt"), "w")

    for component_id, component in enumerate(ordered_dependency_list):
        new_component = []
        fw = open(os.path.join(subfamily_details_dir, str(cluster_id) + "-Sub" + str(component_id+1) + "_bp_ann.txt"), "w")
        fw_bp_details.write(str(cluster_id) + "-Sub" + str(component_id+1) + ':\n')
        fw_bp_details.write('Total motifs: ' + str(len(component)) + '\n')
        bp_dict = {}

        for component_loop_id, ((i, r1), parent) in enumerate(component):
            load_id, loop_id = load_id_dict[(i, r1)]

            f_loop = open(os.path.join(loop_dir, r1 + ".smf"))
            loop_info_lines = f_loop.readlines()
            f_loop.close()

            t_bp_dict = {}
            in_section = False
            for line in loop_info_lines:
                if line.startswith('#info=basepair'):
                    in_section = True
                elif line.startswith('#info=stacking'):
                    in_section = False
                    break
                elif in_section == True:
                    _, bp, orientation = line.strip().split(',')

                    bp_edges = bp.strip().split('/')
                    bp_edges.sort()
                    bp = '/'.join(bp_edges)

                    if len(orientation) > 0:
                        bp = orientation[0] + bp

                    if bp not in t_bp_dict:
                        t_bp_dict[bp] = 0
                    t_bp_dict[bp] += 1

            for bp in t_bp_dict:
                if bp not in bp_dict:
                    bp_dict[bp] = []
                bp_dict[bp].append(t_bp_dict[bp])

            if output_env == 'local':
                fw.write(str(load_id) + "\n")

            write_alignments_with_interactions_to_file("", r1, "", 0.0, pdb_res_map_dict, fw, False)
            fw.write('\n\n')

        for bp in sorted(bp_dict):
            fw_bp_details.write(bp + ' (' + str(len(bp_dict[bp])) + ' motifs, ' + str(sum(bp_dict[bp])) + ' occurences): ' + ','.join(map(lambda x: str(x), bp_dict[bp])) + '\n')

        fw_bp_details.write('\n')

    fw_bp_details.close()

    # return ordered_dependency_list_with_bp_ann

# def generate_componentwise_bp_annotation_files(cluster_id, ordered_dependency_list, load_id_dict):
#     if generate_bp_ann_files == False:
#         return ordered_dependency_list

#     # mapping_dir = pdb_fasta_mapping_dir
#     # loop_dir = loop_dir
#     output_dir = os.path.join(superimposition_output_dir, 'subfamilywise_bp_ann')
    
#     create_directory(output_dir)

#     ordered_dependency_list_with_bp_ann = []

#     pdb_res_map_dict = {}
#     for component_id, component in enumerate(ordered_dependency_list):
#         new_component = []
#         fw = open(os.path.join(output_dir, str(cluster_id) + "_" + str(component_id+1) + ".txt"), "w")
#         for component_loop_id, ((i, r1), parent) in enumerate(component):
#             load_id, loop_id = load_id_dict[(i, r1)]

#             f_loop = open(os.path.join(loop_dir, r1 + ".smf"))
#             loop_info_lines = f_loop.readlines()
#             f_loop.close()

#             fw.write(str(load_id) + "\t" + r1 + "\t\t")

#             pdb_chain, loop_regions = r1.strip().split(":")

#             if pdb_chain not in pdb_res_map_dict:
#                 pdb_res_map_dict[pdb_chain] = load_pdb_res_map(pdb_chain)

#             pdb_res_map = pdb_res_map_dict[pdb_chain]
#             loop_regions = loop_regions.strip().split("_")
#             loop_regions_seq = "".join(loop_info_lines[1].strip().split("..."))

#             loop_regions_pdb_ind = []
#             # index_residue_dict = {}
#             pdb_index_list = []
#             for part, region in enumerate(loop_regions):
#                 region_pdb_index = []
                
#                 s, e = region.strip().split("-")
#                 s = int(s)
#                 e = int(e)
#                 for fasta_indx in range(s, e+1):
#                     pdb_index_list.append(pdb_res_map[fasta_indx][1])

#                 start_ind = str(pdb_res_map[s][1])
#                 if len(pdb_res_map[s][2]) > 0:
#                     start_ind += '.' + str(pdb_res_map[s][2])

#                 end_ind = pdb_res_map[e][1]
#                 if len(pdb_res_map[e][2]) > 0:
#                     end_ind += '.' + str(pdb_res_map[e][2])

#                 loop_regions_pdb_ind.append(start_ind + "-" + end_ind)

#             index_residue_dict = dict(zip(pdb_index_list, loop_regions_seq))

#             r1_pdb_ind = pdb_chain + ":" + "_".join(loop_regions_pdb_ind)
#             fw.write(r1_pdb_ind + "\n")
            
#             bp_ann_list = []
#             for line in loop_info_lines:
#                 if line.startswith(">"):
#                     continue
#                 elif line.startswith("#info=stacking"):
#                     break
#                 else:
#                     pieces = line.strip().split(",")
#                     if len(pieces) == 3:
#                         a, b = pieces[0].strip().split("-")
#                         a = pdb_index_list[int(a)]
#                         b = pdb_index_list[int(b)]
#                         bp = a + "-" + b
#                         line = bp + "," + index_residue_dict[a] + "-" + index_residue_dict[b] + "," + pieces[1] + "," + pieces[2] + "\n"
#                         bp_ann_list.append(((a, b), index_residue_dict[a] + "-" + index_residue_dict[b], pieces[1], pieces[2]))
                    
#                     fw.write(line)
#                     if "..." in line:
#                         print_a_dict_sorted(index_residue_dict, fw, separator=": ")
            
#             new_component.append(((i, r1), parent, index_residue_dict, bp_ann_list))
#             fw.write("\n")
#         fw.close()

#         ordered_dependency_list_with_bp_ann.append((component_id, new_component))

#     return ordered_dependency_list_with_bp_ann

def get_rmsd_align_len(i, r1, j, r2, cluster_pairwise_alignment_details):
    avg_rmsd, pairwise_align_details = cluster_pairwise_alignment_details[(i, r1)]
    for (j_c, r2_c, rmsd, align_length) in pairwise_align_details:
        if j_c == j and r2 == r2_c:
            return rmsd, align_length

    return 0, 0

def get_family_rmsd_and_alignment_summary(ordered_dependency_list, cluster_pairwise_alignment_details, subfamily_only = False):
    rmsd_align_len_list = []
    for component_id, component in enumerate(ordered_dependency_list):
        is_first = True
        for component_loop_id, ((i, r1), parent) in enumerate(component):
            if parent != None:
                j, r2 = parent
                if is_first == False or subfamily_only == False:
                    rmsd, align_len = get_rmsd_align_len(i, r1, j, r2, cluster_pairwise_alignment_details)
                    rmsd_align_len_list.append((rmsd, align_len))
            is_first = False
        
    return get_weighted_avg_rmsd(rmsd_align_len_list)


def similar_to_existing_color(subfamily_colors_rgb, new_color):
    for color in subfamily_colors_rgb:
        is_similar = True
        # Check if all values are same or not
        for i, val in enumerate(color):
            if abs(color[i] - new_color[i]) > 0.01:
                is_similar = False
                break

        # Maybe check if value ratio is same or not, if is_similar == False

    return is_similar

def get_sub_family_colors(max_component_count = 100):
    # PyMol and pyplot common colors
    # (pyplot order) 'red','green','blue','brown','firebrick','salmon','darksalmon','chocolate','orange','wheat','olive','yellow','limegreen','aquamarine','teal','cyan','lightblue','skyblue','violet','purple','magenta','pink','lightpink'
    # (sorted) 'aquamarine','blue','brown','chocolate','cyan','darksalmon','firebrick','green','lightblue','lightpink','limegreen','magenta','olive','orange','pink','purple','red','salmon','skyblue','teal','violet','wheat','yellow'
    subfamily_colors_by_name = ['red', 'green', 'blue', 'cyan', 'brown', 'lightpink', 'magenta', 'wheat', 'teal', 'orange', 'purple', 'lightblue', 'violet', 'darksalmon', 'olive', 'yellow']
    subfamily_colors_dict = {'red': [1.0, 0.0, 0.0], 
                            'green': [0.0, 1.0, 0.0], 
                            'blue': [0.0, 0.0, 1.0], 
                            'cyan': [0.0, 1.0, 1.0], 
                            'brown': [0.65, 0.32, 0.17], 
                            'lightpink': [1.00, 0.75, 0.87],                             
                            'purple': [0.75, 0.00, 0.75],
                            'wheat': [0.99, 0.82, 0.65],
                            'teal': [0.00, 0.75, 0.75],
                            'orange': [1.0, 0.5, 0.0],
                            'lightblue': [0.75, 0.75, 1.0],                    
                            'violet': [1.0, 0.5, 1.0],
                            'magenta': [1.0, 0.0, 1.0], 
                            'darksalmon': [0.73, 0.55, 0.52],
                            'olive': [0.77, 0.70, 0.00],
                            'yellow': [1.0, 1.0, 0.0]}

    subfamily_colors_rgb = []
    for color_name in subfamily_colors_by_name:
        subfamily_colors_rgb.append(subfamily_colors_dict[color_name])

    random.seed(3)

    for i in range (max_component_count * 10):
        if len(subfamily_colors_rgb) >= max_component_count:
            break

        rand_ind1 = random.randint(0,len(subfamily_colors_rgb) - 1)
        rand_ind2 = random.randint(0,len(subfamily_colors_rgb) - 1)
        rand_ind3 = random.randint(0,len(subfamily_colors_rgb) - 1)

        # mix 3 random existing color to generate a new color
        new_color = [(x + y + z)/3.0 for x, y, z in zip(subfamily_colors_rgb[rand_ind1], subfamily_colors_rgb[rand_ind2], subfamily_colors_rgb[rand_ind3])]

        if not similar_to_existing_color(subfamily_colors_rgb, new_color):
            subfamily_colors_rgb.append(new_color)
            subfamily_colors_by_name.append(new_color)       

    # If need more class, assign random colors
    for i in range (max_component_count * 100):
        if len(subfamily_colors_rgb) >= max_component_count:
            break

        new_color = [random.random(), random.random(), random.random()]
        if not similar_to_existing_color(subfamily_colors_rgb, new_color):
            subfamily_colors_rgb.append(new_color)
            subfamily_colors_by_name.append(new_color)       

    # Generate tuple to make it compatible with pyplot (as list is needed for pymol)
    subfamily_colors_tuple = []
    for item in subfamily_colors_by_name:
        if type(item) == 'list' and len(item) == 3:
            subfamily_colors_tuple.append((item[0], item[1], item[2]))
        else:
            subfamily_colors_tuple.append(item)

    return subfamily_colors_by_name, subfamily_colors_tuple

def get_component_graph(component_features, load_id_dict, color):
    cycle_node_list = []
    uncycled_node_list = []
    edge_list = []

    _, cycle_node_list_original, component_nodes, component_directed_adjacency_list = component_features

    for node in component_nodes:

        # family_name, loop_id = load_id_dict[node][0].split('_')
        family_name, loop_id = separate_family_name_and_loop_id(load_id_dict[node][0])

        loop_short_name = get_motif_family_short_code(family_name) + '_' + loop_id

        if node in cycle_node_list_original:    
            cycle_node_list.append(loop_short_name)
        else:
            uncycled_node_list.append(loop_short_name)

    for node1 in component_directed_adjacency_list:
        for node2 in component_directed_adjacency_list[node1]:
            
            # family_name1, loop_id1 = load_id_dict[node1][0].split('_')
            family_name1, loop_id1 = separate_family_name_and_loop_id(load_id_dict[node1][0])
            # family_name2, loop_id2 = load_id_dict[node2][0].split('_')
            family_name2, loop_id2 = separate_family_name_and_loop_id(load_id_dict[node2][0])

            loop_short_name1 = get_motif_family_short_code(family_name1) + '_' + loop_id1
            loop_short_name2 = get_motif_family_short_code(family_name2) + '_' + loop_id2
            edge_list.append((loop_short_name1, loop_short_name2, round(component_directed_adjacency_list[node1][node2][0],2)))

    return (cycle_node_list, uncycled_node_list, edge_list, color)

def generate_subfamily_image(image_file_list, pdb_organism_details, cluster_id, subfamily_dir, draw_figures, suffix = '', is_graph_image=False, show_pdb_info=False, show_image_caption=True):
    if draw_figures == False:
        return
    # if create_subfamily_images == False:
        # return

    # superimposition_output_dir = image_file_list[0][:-1*len(os.path.basename(image_file_list[0]))]
    # subfamily_dir = os.path.join(superimposition_output_dir, 'subfamily')

    create_directory(subfamily_dir)

    if len(suffix) > 0:
        create_collage(image_file_list, pdb_organism_details, os.path.join(subfamily_dir, str(cluster_id) + '_' + suffix + '.png'), show_pdb_info, is_graph_image, show_image_caption)
    else:
        create_collage(image_file_list, pdb_organism_details, os.path.join(subfamily_dir, str(cluster_id) + '.png'), show_pdb_info, is_graph_image, show_image_caption)
    # for image_fname in image_file_list:
        # copy_subfamily_image(image_fname, suffix)

def get_index_list(loop):
    ind_list = []
    segments = loop.strip().split(':')[1].strip().split('_')
    for segment in segments:
        a, b = segment.strip().split('-')
        ind_list.append((int(a), int(b)))

    return ind_list

def union_segments(r1_union_ind_list, r1_ind_list, cr1_ind_list):
    index = 0
    for (a, b) in r1_ind_list:
        if a > b:
            logger.error('ERROR: Range is reversed (' + str(a) + '-' + str(b) + ')')
            sys.exit()
        for (c, d) in cr1_ind_list:
            if c > d:
                logger.error('ERROR: Range is reversed (' + str(c) + '-' + str(d) + ')')
                sys.exit()
            # if a <= c and d <= b:
            if c <= b and d >= a: # has overlap
                overlap_start = max(a, c)
                overlap_end = min(b, d)                
                (x, y) = r1_union_ind_list[index]                
                if overlap_start < x:
                    r1_union_ind_list[index] = (overlap_start, y)
                if overlap_end > y:
                    r1_union_ind_list[index] = (r1_union_ind_list[index][0], overlap_end)
        index += 1

    return r1_union_ind_list

def generate_loop_boundary(items, cluster_alignment_data):
    # print cluster_alignment_data
    loop_boundary = {}
    loop_boundary_original = {}

    for (i, r1) in sorted(items, key=lambda x: items[x][0]):
        r1_ind_list = get_index_list(r1)
        loop_boundary_original[(i, r1)] = r1_ind_list
        if len(items) > 1:
            r1_union_ind_list = []    
            for (a, b) in r1_ind_list:
                r1_union_ind_list.append((b, a))
            target_loop = strToNode(r1)
            _, pairwise_align_details = items[(i, r1)]
            for (j, r2, _, _) in pairwise_align_details:
                mobile_loop = strToNode(r2)
                (_, _, _, cr1, cr2, _, _, _) = cluster_alignment_data[target_loop][mobile_loop]
                cr1_ind_list = get_index_list(cr1)
                r1_union_ind_list = union_segments(r1_union_ind_list, r1_ind_list, cr1_ind_list)
            r1_final_union_ind_list = []    
            for (a, b) in r1_union_ind_list:
                if a <= b:
                    r1_final_union_ind_list.append((a,b))

            loop_boundary[(i, r1)] = r1_final_union_ind_list
        else:
            loop_boundary[(i, r1)] = r1_ind_list

    return loop_boundary, loop_boundary_original

def get_loop_boundary_pdb_index(loop_boundary_fasta):
    loop_boundary_pdb = {}
    for (i, r1) in loop_boundary_fasta:
        loop_boundary_pdb[(i, r1)] = []
        pdb_chain = r1.strip().split(':')[0]
        pdb_res_map = load_pdb_res_map(pdb_chain)
        for (a, b) in loop_boundary_fasta[(i, r1)]:
            for ind in range(a, b+1):
                if ind in pdb_res_map:
                    loop_boundary_pdb[(i, r1)].append(pdb_res_map[ind])

    return loop_boundary_pdb

def reset_pymol():
    pymol.cmd.sync()
    # pymol.commanding.sync()
    pymol.cmd.deselect()
    pymol.cmd.delete('all')
    pymol.cmd.reinitialize()
    pymol.cmd.bg_color('white')

def load_pdb_in_pymol(partial_pdbx_dir, pdb_chain, pdb_align_ind_list, pdb_disp_ind_list, backbone_atom_list, sugar_atom_list, load_name_id, is_cif, loop_name):
    pdb_id, chain_id = pdb_chain.strip().split('_')
    pdb_load_name = 'pdb_' + str(load_name_id)
    chain_load_name = 'rna_' + str(load_name_id)
    target_load_name = 'target_' + str(load_name_id)
    display_load_name = 'display_target_' + str(load_name_id)

    if is_cif:
        # pymol.cmd.load(os.path.join(pdb_dir, pdb_id+'.cif'), pdb_load_name)
        loop_name = str(strToNode(loop_name))
        pymol.cmd.load(os.path.join(partial_pdbx_dir, loop_name+'.cif'), pdb_load_name)
    else:
        pymol.cmd.load(os.path.join(partial_pdbx_dir, pdb_id+'.pdb'), pdb_load_name)

    pymol.cmd.select(chain_load_name, pdb_load_name + ' and chain %s' % chain_id)
    pymol.cmd.select(target_load_name, chain_load_name + ' and (%s)' % ' or '.join(list(map(lambda x: 'resi '+x, pdb_align_ind_list))))
    pymol.cmd.select(display_load_name, chain_load_name + ' and (%s)' % ' or '.join(list(map(lambda x: 'resi '+x, pdb_disp_ind_list))))

    pymol.cmd.hide('everything', pdb_load_name) 
    # pymol.commanding.sync()

    centroid = get_centroid_of_loop(target_load_name, pdb_align_ind_list, backbone_atom_list, sugar_atom_list)

    return pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid
    # return test, pdb_load_name

def get_ordered_loop_coordinates(coordinates, atom_list, res, target_load_name):
    stored.sel = []
    for atom in atom_list[res]:
        pymol.cmd.select('t', target_load_name + ' and resi %s and name %s' % (res, atom))
        pymol.cmd.iterate_state(1, 't', 'stored.sel.append([x,y,z])')
    
    if len(stored.sel) > 0:
        coordinates.append(numpy.sum(stored.sel, axis=0)/float(len(stored.sel)))

def get_ordered_loop_backbone_coordinates(target_load_name, pdb_align_ind_list, backbone_atom_list, sugar_atom_list):
    coordinates = []
    for res in pdb_align_ind_list:
        get_ordered_loop_coordinates(coordinates, backbone_atom_list, res, target_load_name)
        get_ordered_loop_coordinates(coordinates, sugar_atom_list, res, target_load_name)

    return coordinates

def get_centroid_of_loop(target_load_name, pdb_align_ind_list, backbone_atom_list, sugar_atom_list):
    coord_list = get_ordered_loop_backbone_coordinates(target_load_name, pdb_align_ind_list, backbone_atom_list, sugar_atom_list)
    centroid = numpy.sum(coord_list, axis=0)/float(len(coord_list))
    return centroid

def get_seqnums_from_indices(indices):
    seqnums = []
    # print(indices)
    for chain, index, icode in indices:
        seqnums.append(index + icode)

    return seqnums

def load_one_pdb_in_pymol(partial_pdbx_dir, i, r1, loop_boundary_dict, loop_boundary_original_dict, load_name_id, is_cif, load_color):
    pdb_chain = r1.strip().split(':')[0]
    align_target = get_seqnums_from_indices(loop_boundary_dict[(i, r1)])
    align_mobile = get_seqnums_from_indices(loop_boundary_original_dict[(i, r1)])

    backbone_atom_list = {}
    sugar_atom_list = {}
    for index in align_target:
        # backbone_atom_list[index] = ['P']
        # sugar_atom_list[index] = ["C1'"]
        backbone_atom_list[index], sugar_atom_list[index] = get_backbone_and_sugar_atoms()

    pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid = load_pdb_in_pymol(partial_pdbx_dir, pdb_chain, align_target, align_mobile, backbone_atom_list, sugar_atom_list, load_name_id, is_cif, r1)
    # pymol.commanding.sync()
    pymol.cmd.sync()
    pymol.cmd.color(load_color, chain_load_name)

    return (pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, r1, load_name_id)

def translate_coords(pdb_data, target_load_name, pdb_align_ind_list, backbone_atom_list, sugar_atom_list, centroid):
    coord_list = get_ordered_loop_backbone_coordinates(target_load_name, pdb_align_ind_list, backbone_atom_list, sugar_atom_list)
    coord_list -= centroid
    pdb_translated = pdb_data - centroid

    return coord_list, pdb_translated

def alter_structure(pdb_translated, pdb_load_name):
    stored.res_list = pdb_translated.tolist()
    pymol.cmd.alter_state(1, pdb_load_name, '(x,y,z)=stored.res_list.pop(0)')

def translate_and_show_single_loop(pymol_load_info, align_target, disp_target, load_id, image_fname, show_extended_loop, show_label, display_color = 'gray', align_color = 'red'):
    pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, lp, load_name_id = pymol_load_info
    disp_target = get_seqnums_from_indices(disp_target)
    align_target = get_seqnums_from_indices(align_target)

    #load pdb file in pymol
    pdb_data = get_pdb_coordinates(pdb_load_name)
    backbone_atom_list = {}
    sugar_atom_list = {}
    for index in align_target:
        # backbone_atom_list[index] = ['P']
        # sugar_atom_list[index] = ["C1'"]
        backbone_atom_list[index], sugar_atom_list[index] = get_backbone_and_sugar_atoms()

    #Translate coordinates around centroid (to be used to find superimposition matrix)
    sel, pdb_translated = translate_coords(pdb_data, target_load_name, align_target, backbone_atom_list, sugar_atom_list, centroid)
    alter_structure(pdb_translated, pdb_load_name)
    # display_load_name = alter_and_display_structure(pdb_translated, load_id, pdb_load_name, chain_load_name, disp_target, color1, target_load_name, color2, show_label, "C2'")
    # display_load_name = alter_and_display_structure(pdb1_translated, load_id + '_1', pdb_load_name, chain_load_name, disp_target, 'gray', target_load_name, 'tv_yellow')

    if draw_input_images == True:
        show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, target_load_name, image_fname, show_extended_loop, show_label, "C2'", display_color, align_color)
    else:
        set_loop_color(display_color, align_color, display_load_name, target_load_name)

    # pymol.cmd._do('zoom')

    # pymol.commanding.sync()
    # pymol.cmd.png(image_fname, 1200, 1200, dpi=300, ray=1, quiet=0)
    # pymol.commanding.sync()

    return (display_load_name, target_load_name, chain_load_name), (pdb_load_name, chain_load_name, target_load_name, display_load_name, np.zeros(3), pdb_chain, lp, load_name_id)

def config_pymol_cartoon(display_color, show_label):
    # arrow,dash,loop,putty,skip,automatic,dumbbell,oval,rectangle,tube
    # pymol.cmd.cartoon('dumbbell')     
    # pymol.cmd._do('set cartoon_nucleic_acid_color, ' + display_color)
    # pymol.cmd._do('cartoon oval')

    # pymol.cmd._do('set cartoon_ring_color, '+ display_color)

    if show_label:
        # pymol.cmd._do('set cartoon_ring_mode, 1') # (or 2 or 3) 
        # pymol.cmd._do('set cartoon_ring_transparency, 0.5')
        pymol.cmd.set(name='cartoon_ring_mode',value=1,quiet=1)
        pymol.cmd.set(name='cartoon_ring_transparency',value=0.5,quiet=1)
    else:
        pymol.cmd.set(name='cartoon_ring_mode',value=0,quiet=1)

    # pymol.cmd._do('set cartoon_tube_radius,0.8')
    # pymol.cmd._do('set cartoon_ring_finder, 1') # (or 2 or 3 or 4)
    # pymol.cmd._do('set cartoon_nucleic_acid_mode, 4') # (or 1 or 2 or 3 or 4)
    # pymol.cmd._do('set cartoon_fancy_helices=1')
    # pymol.cmd._do('set cartoon_highlight_color, gray')

    # pymol.cmd._do('rotate y, 145, all')
    # pymol.cmd._do('rotate z, -45, all')
    pass

def set_loop_color(display_color, align_color, display_load_name, align_load_name):

    if type(display_color) is list and len(display_color) == 3:
        pymol.cmd.set_color(display_load_name + '_color', display_color)
        display_color = display_load_name + '_color'

    # print(display_color, align_color)

    if type(align_color) is list and len(align_color) == 3:
        pymol.cmd.set_color(align_load_name + '_color', align_color)
        align_color = align_load_name + '_color'

    pymol.cmd.color(display_color, display_load_name)
    pymol.cmd.color(align_color, align_load_name)
        
def show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, label_atom, display_color = 'gray', align_color = 'red'):

    set_loop_color(display_color, align_color, display_load_name, align_load_name)

    if show_label:
        pymol.cmd.label(display_load_name + " and name " + label_atom, "'%s-%s' %(resn, resi)")
        config_pymol_cartoon(display_color, show_label)
    # else:
        # pymol.cmd.label(display_load_name + " and name " + "dummy", "'%s-%s' %(resn, resi)")

    if show_extended_loop:
        pymol.cmd.show('cartoon', chain_load_name)
        # print('showing cartoon')
    else:
        pymol.cmd.show('cartoon', display_load_name)

    # pymol.cmd._do('zoom')
    pymol.cmd.sync()
    pymol.cmd.zoom()

    # pymol.commanding.sync()
    pymol.cmd.sync()
    # pymol.cmd._do('set ray_opaque_background, 0')
    # pymol.cmd.set(name='ray_opaque_background',value=0,quiet=1)
    pymol.cmd.png(image_fname, 1200, 1200, dpi=300, ray=1, quiet=1)
    # pymol.commanding.sync()
    # wait_for_certain_files_to_be_generated([image_fname], True)
    pymol.cmd.sync()

def get_rotatation_matrix(axis = 'x', angle = 0):
    cos_t = math.cos(math.radians(angle))
    sin_t = math.sin(math.radians(angle))

    mat = np.zeros((3,3))
    
    if axis.lower() == 'x':
        mat[0][0] = 1
        mat[1][1] = cos_t
        mat[2][2] = cos_t
        mat[1][2] = -1 * sin_t
        mat[2][1] = sin_t

    elif axis.lower() == 'y':
        mat[1][1] = 1
        mat[0][0] = cos_t
        mat[2][2] = cos_t
        mat[0][2] = sin_t
        mat[2][0] = -1 * sin_t

    elif axis.lower() == 'z':
        mat[2][2] = 1
        mat[0][0] = cos_t
        mat[1][1] = cos_t
        mat[0][1] = -1 * sin_t
        mat[1][0] = sin_t
    else:
        logger.error('Invalid axis. Please choose x, y, or z.')
        sys.exit()

    return mat

def get_multiple_orientation_rotation_matrices():
    rotation_matrices = []
    rotation_matrices.append(get_rotatation_matrix())
    # generate multiple orientation when any specific view file not found
    
    # if r1_view == None and generate_multiple_orientation == True:
    if generate_multiple_orientation:
        rotation_matrices.append(numpy.dot(get_rotatation_matrix('x',45), get_rotatation_matrix('y',45))) # x: 45, y: 45
        # rotation_matrices.append(get_rotatation_matrix('x',90)) # x: 90
        # rotation_matrices.append(get_rotatation_matrix('y',90)) # x: 90

    # rotation_matrices.append(numpy.dot(get_rotatation_matrix('x',90), get_rotatation_matrix('y',90))) # x: 135, y: 135

    return rotation_matrices

# def get_rotation_matrices():
#     rotation_matrices = []
#     # rotation_matrices.append(get_rotatation_matrix())
#     rotation_start, rotation_end, step = rotation_angle_start_end_step

#     for angle_x in range(rotation_start, rotation_end, step):
#         for angle_y in range(rotation_start, rotation_end, step):
#             for angle_z in range(rotation_start, rotation_end, step):
#                 rotation_matrices.append(numpy.dot(get_rotatation_matrix('x', angle_x), get_rotatation_matrix('y', angle_y), get_rotatation_matrix('z', angle_z)))

#     return rotation_matrices

# def generate_matrices_and_loop_rotations(ordered_dependency_list):
#     first_motif = None
#     if len(ordered_dependency_list) > 0 and len(ordered_dependency_list[0]) > 0:
#         first_motif = ordered_dependency_list[0][0][1]
#     else:
#         logger.warning('No motif found in component.')
#         sys.exit()

#     rotation_matrices = get_rotation_matrices()
#     for v, rotation_matrix in enumerate(rotation_matrices):
#         rotation_version = 'rotation_v' + str('{:05d}'.format(v+1))

#         ### Rotate, group and show the subfamilies 
        
#             pymol.cmd._do('hide all')
#             pymol.commanding.sync()
#             time.sleep(.100)
#             pdb_load_name1, chain_load_name1, target_load_name1, display_load_name1, centroid1, pdb_chain1, lp1, load_name_id = pymol_load_info_dict[(i, r1)]
#         ### superimposition subfamilies
#         # file_name = os.path.join(superimposition_output_dir,  str(cluster_id) + '_' + str(component_id + 1) + '_' + str(component_loop_id + 1) + '_' + str(family_loop_id))
#         image_fname = os.path.join(rotated_loop_image_dir, rotation_version + '_' + load_name_id + '_3.png')
#         text_fname = os.path.join(superimposition_output_dir, str(cluster_id) + '_' + str(component_id + 1) + '_' + str(component_loop_id + 1) + '_' + str(family_loop_id) + '.txt')

#         # Draw the first loop of the first  component  independent of any other loops
#         if component_id == 0 and component_loop_id == 0:
#             if draw_pymol_figure:
#                 display_load_name, align_load_name, chain_load_name = loop_display_info_dict[(i,r1)]
#                 rotate_first_loop(pymol_load_info_dict[(i, r1)], rotation_matrix)
#                 show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_label, "C2'", 'gray', subfamily_colors[component_id])

def get_pdb_coordinates(pdb_load_name):
    stored.pdb = []
    pymol.cmd.iterate_state(1, pdb_load_name, 'stored.pdb.append([x,y,z])')
    return stored.pdb

# Rotate the first loop of the cluster to define the orientation
def rotate_first_loop(pymol_load_info, rotation_matrix):
    pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, lp, load_name_id = pymol_load_info
    pdb_data = get_pdb_coordinates(pdb_load_name)
    pdb_translated = pdb_data - centroid
    pdb_rotated = numpy.dot(pdb_translated, rotation_matrix)
    alter_structure(pdb_rotated, pdb_load_name)

def write_pdb_organisgm_details(pdb_organism_details, loop, fp):
    pdb_chain = loop.strip().split(':')[0]
    if pdb_chain in pdb_organism_details:
        RNA_Types, organism, org_class, org_type, pdb_source = pdb_organism_details[pdb_chain]
        fp.write(pdb_chain + '\t' + RNA_Types + '\t' + organism + '\t' + org_class + '\t' + org_type + '\t' + pdb_source + '\t')
    
    else:
        fp.write(pdb_chain + '\t\t\t\t\t\t')

def write_dock_file_list(c_id, m_id, i, r1, j, r2, align_len, zscore, score, fp1, cluster_pairwise_alignment_details, pdb_organism_details, t1, t2):
    rmsd = 0.0

    # fp1 = open(text_fname, 'a')

    if(r2 != 'None'):
        rmsd, align_len2 = get_rmsd_align_len(i, r1, j, r2, cluster_pairwise_alignment_details)
        if align_len != align_len2 :
            logger.error('ERROR!!! Align length inconsistency!!!')
            # print('Error message generated from script rmsd_based_benchmarking.py')
            sys.exit()

    if(r2 != 'None'):
        fp1.write(str(c_id) + '\t')
        write_pdb_organisgm_details(pdb_organism_details, r1, fp1)
        fp1.write(str(m_id) + ',' + '(' + str(i) + ',' + r1 + ') <- (' + str(j) + ',' + r2 + '), rmsd: ' + str(round(rmsd, 3)) + ', align len: ' + str(align_len) + ', score: ' + str(score)+ ', zscore: ' + str(zscore))
        fp1.write('\t' + t1 + ',' + t2 + '\n')

    else:
        if c_id == 1:
            # Write Title
            fp1.write('c_id \t PDB_chain \t RNA_Types \t Organism \t Org_class \t Org_type \t PDB_source \t Loop \t Align Order \n')

        fp1.write(str(c_id) + '\t')
        write_pdb_organisgm_details(pdb_organism_details, r1, fp1)
        fp1.write(str(m_id) + ',' + '(' + str(i) + ',' + r1 + ')' + '\n')

    # print(text_fname)

    # fp1.close()

def get_pdb_residues(pdb_dir, pdb_chain, residue_list, structures, is_cif, loop):
    # pdb_dir = '../../DataCollection/nrpdbs'
    pdb_id, chain_id = pdb_chain.strip().split('_')
    # print(pdb_id, chain_id)
    # print(loop)

    residue_dict = None
    # if (pdb_id, chain_id) in structures:
        # residue_dict = structures[(pdb_id, chain_id)]

    if loop in structures:
        residue_dict = structures[loop]
    else:
        pdb_fn = None
        if is_cif:
            pdb_fn = os.path.join(pdb_dir, str(strToNode(loop)) + '.cif')
            # pdb_fn = os.path.join(pdb_dir, pdb_id + '.cif')
        else:
            pdb_fn = os.path.join(pdb_dir, pdb_id + '.pdb')

        parser = None
        if is_cif:
            parser = FastMMCIFParser()
        else:
            parser = PDBParser()
        
        # print(pdb_fn)
        structure = parser.get_structure('struct', pdb_fn)      
        
        chain = structure[0][chain_id]
        residues = chain.get_residues()

        residue_dict = {}
        
        for res in residues:
            hetflag, ind, icode = res.get_id()
            residue_dict[(ind, icode)] = res

        structures[loop] = residue_dict
        # structures[(pdb_id, chain_id)] = residue_dict

    return residue_dict

def get_atom_list(atom_names, residue_dict, residue_list):
    # residues = chain.get_residues()
    atom_list_dict = {}

    for index in residue_list:
        # if chain_id == 'n' and index == 'a':
        #     continue
        ind, icode = get_separated_index_icode(index)
        
        atom_list = []

        if (ind, icode) in residue_dict:
            for atom in atom_names:
                if atom in residue_dict[(ind, icode)]:
                    atom_list.append(atom)

            atom_list_dict[index] = atom_list

    return atom_list_dict

def get_backbone_atom_list(pdb_dir, pdb_chain, residue_list, structures, is_cif, loop):

    backbone_atoms, sugar_atoms = get_backbone_and_sugar_atoms()

    residue_dict = get_pdb_residues(pdb_dir, pdb_chain, residue_list, structures, is_cif, loop)

    backbone_atom_list = get_atom_list(backbone_atoms, residue_dict, residue_list)
    sugar_atom_list = get_atom_list(sugar_atoms, residue_dict, residue_list)

    return backbone_atom_list, sugar_atom_list

def get_common_atom_list(atom_list1, atom_list2, pdb1_align_ind_list, pdb2_align_ind_list, lp1, lp2):

    common_atom_list1 = {}
    common_atom_list2 = {}

    for i, index1 in enumerate(pdb1_align_ind_list):
        index2 = pdb2_align_ind_list[i]
        common_atom_list1[index1] = []
        common_atom_list2[index2] = []      
        if index1 not in atom_list1 or index2 not in atom_list2:
            # print 'residue index ' + str(index2) + ' not found in ' + lp2
            continue
            # sys.exit()
        else:
            for atom in atom_list1[index1]:
                if atom in atom_list2[index2]:
                    common_atom_list1[index1].append(atom)
                    common_atom_list2[index2].append(atom)

    return common_atom_list1, common_atom_list2

def collect_common_backbone_atom_list(pdb_dir, pdb_chain1, pdb_chain2, pdb1_align_ind_list, pdb2_align_ind_list, is_cif, structures, lp1, lp2):
    # pdb_id, chain_id = pdb_chain.strip().split('_')
    backbone_atom_list1, sugar_atom_list1 = get_backbone_atom_list(pdb_dir, pdb_chain1, pdb1_align_ind_list, structures, is_cif, lp1)
    backbone_atom_list2, sugar_atom_list2 = get_backbone_atom_list(pdb_dir, pdb_chain2, pdb2_align_ind_list, structures, is_cif, lp2)
    
    common_backbone_atom_list1, common_backbone_atom_list2 = get_common_atom_list(backbone_atom_list1, backbone_atom_list2, pdb1_align_ind_list, pdb2_align_ind_list, lp1, lp2)
    common_sugar_atom_list1, common_sugar_atom_list2 = get_common_atom_list(sugar_atom_list1, sugar_atom_list2, pdb1_align_ind_list, pdb2_align_ind_list, lp1, lp2)

    return common_backbone_atom_list1, common_backbone_atom_list2, common_sugar_atom_list1, common_sugar_atom_list2

def rotate_pdb(sel1, sel2, pdb_translated):
    e0 = numpy.sum( numpy.sum(sel1 * sel1,axis=0),axis=0) + numpy.sum( numpy.sum(sel2 * sel2,axis=0),axis=0)
    v, s, wt = numpy.linalg.svd( numpy.dot( numpy.transpose(sel2), sel1))
    reflect = float(str(float(numpy.linalg.det(v) * numpy.linalg.det(wt))))

    if reflect == -1.0:
        s[-1] = -s[-1]
        v[:,-1] = -v[:,-1]

    u = numpy.dot(v, wt)

    return numpy.dot((pdb_translated), u)

# Rotates the second loop
def rotate_loop(partial_pdbx_dir, pymol_load_info1, pymol_load_info2, align_target, align_mobile, disp_mobile, is_cif=True):

    # Stores the atoms of loops
    structures = {}

    align_target = get_seqnums_from_indices(align_target)
    align_mobile = get_seqnums_from_indices(align_mobile)
    disp_mobile = get_seqnums_from_indices(disp_mobile)

    pdb_load_name1, chain_load_name1, target_load_name1, display_load_name1, centroid1, pdb_chain1, lp1, load_name_id = pymol_load_info1

    pdb_load_name2, chain_load_name2, target_load_name2, display_load_name2, centroid2, pdb_chain2, lp2, load_name_id = pymol_load_info2

    #Make a list  of common atoms for all the aligned residues
    backbone_atom_list1, backbone_atom_list2, sugar_atom_list1, sugar_atom_list2 = collect_common_backbone_atom_list(partial_pdbx_dir, pdb_chain1, pdb_chain2, align_target, align_mobile, is_cif, structures, lp1, lp2)

    #load 1st pdb data from pymol
    pdb1_data = get_pdb_coordinates(pdb_load_name1)
    #Translate coordinates around centroid (to be used to find superimposition matrix)
    # sel1 = get_ordered_loop_backbone_coordinates(chain_load_name1, align_target, backbone_atom_list1, sugar_atom_list1)
    sel1, pdb1_translated = translate_coords(pdb1_data, chain_load_name1, align_target, backbone_atom_list1, sugar_atom_list1, centroid1)

    # load 2nd pdb data from pymol and translate
    pdb2_data = get_pdb_coordinates(pdb_load_name2)
    # sel2 = get_ordered_loop_backbone_coordinates(chain_load_name2, align_mobile, backbone_atom_list2, sugar_atom_list2)
    centroid2 = get_centroid_of_loop(target_load_name2, align_mobile, backbone_atom_list2, sugar_atom_list2)
    sel2, pdb2_translated = translate_coords(pdb2_data, chain_load_name2, align_mobile, backbone_atom_list2, sugar_atom_list2, centroid2)

    #rotate the 2nd pdb to align with the 1st one
    pdb2_rotated = rotate_pdb(sel1, sel2, pdb2_translated)

    #Adjust coordinate to fit to the centroid shift due to local alignment
    centroid_resolve = get_centroid_of_loop(target_load_name1, align_target, backbone_atom_list1, sugar_atom_list1)
    pdb2_rotated += centroid_resolve

    # display_load_name_cur = alter_and_display_structure(pdb2_rotated, load_id + "_2", pdb_load_name2, chain_load_name2, disp_mobile, 'gray', target_load_name2, color, show_label, "O4'")
    alter_structure(pdb2_rotated, pdb_load_name2)

def generate_and_add_family_image(superimposition_output_dir, cluster_id, image_file_list, component_list, display_load_name_list, rotation_version, draw_figures, show_extended_loop):

    # Show each components

    # Show all loops in one image
    image_fname = os.path.join(superimposition_output_dir, add_rotation_version_prefix(rotation_version) + cluster_id + '__all.png')
    total_loop_count = 0
    for component in component_list:
        total_loop_count += len(component)
    if draw_figures:
        # pymol.cmd._do('hide all')
        # pymol.cmd.hide()
        # pymol.cmd.sync()
        if show_extended_loop:
            pymol.cmd.show('cartoon', 'all')
        else:
            # for component in display_load_name_list:
            for component in component_list:
                for (i, r1), parent in component:
                    display_load_name, _, _ = display_load_name_list[(i, r1)]
                    pymol.cmd.show('cartoon', display_load_name)
        
        wait_for_certain_time_according_to_wait_factor(total_loop_count)
        # pymol.cmd._do('zoom')
        pymol.cmd.zoom()
        # pymol.cmd.zoom('everything')
        pymol.cmd.sync()
        pymol.cmd.png(image_fname, 1200, 1200, dpi=300, ray=1, quiet=1)
        wait_for_certain_files_to_be_generated([image_fname], False)
        pymol.cmd.sync()

    image_file_list.append((None, image_fname, total_loop_count))
    return image_file_list

    # new_image_file_list = [image_fname]
    # new_image_file_list += image_file_list
    # return new_image_file_list

def get_loop_count_for_rmsd_data_dict(current_rmsd_data_dict):
    loop_count = 0
    for cluster_id in sorted(current_rmsd_data_dict):
        _, cluster_pairwise_alignment_details = current_rmsd_data_dict[cluster_id]
        loop_count += len(cluster_pairwise_alignment_details)
    return loop_count

def write_rmsd_and_alignment_summary(rmsd_and_alignment_summary_dict, current_rmsd_data_dict):
    
    rmsd_align_len_list = []
    # max_cid_len = max([len(x) for x in rmsd_and_alignment_summary_dict])
    # print('Family Name'.ljust(max_cid_len) + '\tAvg. RMSD\tTotal Aln Len')
    for cluster_id in rmsd_and_alignment_summary_dict:
        # rmsd, aln_len = rmsd_and_alignment_summary_dict[cluster_id]
        # print(str(cluster_id).ljust(max_cid_len) + '\t' + "{:.3f}".format(round(rmsd, 3)).rjust(9) + '\t' + str(aln_len).rjust(13))
        rmsd_align_len_list.append(rmsd_and_alignment_summary_dict[cluster_id])

    avg_rmsd, total_alignment_length = get_weighted_avg_rmsd(rmsd_align_len_list)
    print('Total Loops: ' + str(get_loop_count_for_rmsd_data_dict(current_rmsd_data_dict)))
    print('Total superimposition average RMSD: ' + str(round(avg_rmsd, 3)) + '\n')
    # print('Total Superimposition Alignment Length: '.ljust(24) + str(total_alignment_length))

def load_view_file(ordered_dependency_list):
    if len(ordered_dependency_list) < 1 or len(ordered_dependency_list[0]) < 1:
        logger.info('No motif found for superimposition in ' + cluster_id + '.')
        return None

    (i, r1), _ = ordered_dependency_list[0][0]

    view_fn = os.path.join(views_dir, str(r1) + '.view')
    if os.path.isfile(view_fn):
        
        if input_index_type == 'pdb':
            logger.info('View file found for ' + convert_a_loop_from_FASTA_to_PDB(r1) + '. Loading view of this motif to set orientation for all motifs in this cluster.')
        else:
            logger.info('View file found for ' + r1 + '. Loading view of this motif to set orientation for all motifs in this cluster.')

        fv = open(view_fn)
        view_lines = fv.readlines()
        fv.close()
        return view_lines[0]

    return None

# def get_view_from_user(partial_pdbx_dir, cluster_id, i, r, loop_boundary_dict, loop_boundary_original_dict, load_name_id, show_extended_loop, is_cif, load_color):
def get_view_from_user(partial_pdbx_dir, cluster_id, ordered_dependency_list, loop_boundary_dict, loop_boundary_original_dict, show_extended_loop, is_cif):    
    if len(ordered_dependency_list) < 1 or len(ordered_dependency_list[0]) < 1:
        logger.info('No loop found for superimposition in ' + cluster_id + '.')
        return

    (i, r), _ = ordered_dependency_list[0][0]
    load_name_id = 'first_loop'

    display_color = 'gray'
    align_color = 'red'
    

    # pymol.finish_launching()
    # pymol.cmd.hide('all')
    # pymol.cmd.hide()

    (pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, r, load_name_id) = load_one_pdb_in_pymol(partial_pdbx_dir, i, r, loop_boundary_dict, loop_boundary_original_dict, load_name_id, is_cif, display_color)
    pymol.cmd.color(align_color, target_load_name)
    pymol.cmd.deselect()

    view_fname = os.path.join(views_dir, r + '.view')
    if os.path.isfile(view_fname):

        if input_index_type == 'pdb':
            logger.info('View file found for ' + convert_a_loop_from_FASTA_to_PDB(r))
        else:
            logger.info('View file found for ' + r)

        fp = open(view_fname)
        lines = fp.readlines()
        fp.close()
        if len(lines) > 0:
            pymol.cmd.set_view(lines[0])
    
    if show_extended_loop:
        pymol.cmd.show('cartoon', chain_load_name)
        # pymol.cmd.show('cartoon', pdb_load_name)
    else:
        pymol.cmd.show('cartoon', display_load_name)
    
    pymol.cmd.zoom(chain_load_name)

    set_adj_res_view(r, (pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, r, load_name_id))
    
    inp = 'N'
    # while(True):
    print('Please provide desired orientation for ' + cluster_id + '.')
    print('(Yes: To save current orientation / No: To keep previous orientation)')
    inp = input('Continue? (Yes/No): ')
    inp = inp.lower()
    if inp == 'y' or inp == 'yes':
        view = pymol.cmd.get_view()
        # print('got view')
        # print(view)
        fp = open(view_fname, 'w')
        fp.write(str(view))
        fp.close()

    # pymol.cmd.hide('all')
    pymol.cmd.hide()
    pymol.cmd.delete(pdb_load_name)
    pymol.cmd.delete(chain_load_name)
    pymol.cmd.delete(target_load_name)
    pymol.cmd.delete(display_load_name)

    wait_for_certain_time_according_to_wait_factor(1)
    pymol.cmd.sync()
    # sys.exit()

def generate_formatted_superimposition_details(superimposition_details_dir, cluster_id, prev_text_name, pdb_organism_details):
    fp_t = open(prev_text_name)
    fp_superimposition_details = open(os.path.join(superimposition_details_dir, str(cluster_id) + '_superimposition_details.tsv'), 'w')
    if len(pdb_organism_details) == 0:
        fp_superimposition_details.write('Subfamily Id\tLoop (Child)\tSuperimposition Reference (Parent)\tAlignment/Superimposition Details')
    else:
        fp_superimposition_details.write('Subfamily Id\tPDB Chain\tRNA Types\tOrganism Name\tLoop (Child)\tSuperimposition Reference (Parent)\tAlignment/Superimposition Details')

    if input_index_type == 'fasta':
        fp_superimposition_details.write('\tLoop (Child) (FASTA)\tSuperimposition Reference (Parent) (FASTA)')

    fp_superimposition_details.write('\n')
    lines = fp_t.readlines()
    fp_t.close()

    pieces = lines[1].split('\t')
    child_loop = pieces[-1].split(',')[-1].strip().split(')')[0].strip()
    # if input_index_type == 'pdb':
    child_loop_pdb = convert_a_loop_from_FASTA_to_PDB(child_loop)

    if len(pdb_organism_details) == 0:
        fp_superimposition_details.write('\t'.join(pieces[:1]) + '\t' + child_loop_pdb + '\t\t')
    else:
        fp_superimposition_details.write('\t'.join(pieces[:4]) + '\t' + child_loop_pdb + '\t\t')

    if input_index_type == 'fasta':
        fp_superimposition_details.write('\t' + child_loop)

    fp_superimposition_details.write('\n')

    for line in lines[2:]:
        pieces = line.split('\t')
        child_loop, parent_loop = pieces[-1].split(',')
        child_loop = child_loop.strip()
        parent_loop = parent_loop.strip()

        # if input_index_type == 'pdb':
        child_loop_pdb = convert_a_loop_from_FASTA_to_PDB(child_loop)
        parent_loop_pdb = convert_a_loop_from_FASTA_to_PDB(parent_loop)

        # print(pieces[7])
        rmsd, align_len, score, zscore = pieces[7].split(',')[4:8]
        rmsd = str(round(float(rmsd.strip().split(':')[1].strip()), 2))
        align_len = align_len.strip().split(':')[1].strip()
        score = str(round(float(score.strip().split(':')[1].strip()), 2))
        zscore = str(round(float(zscore.strip().split(':')[1].strip()), 2))

        if len(pdb_organism_details) == 0:
            fp_superimposition_details.write('\t'.join(pieces[:1]))
        else:
            fp_superimposition_details.write('\t'.join(pieces[:4]))

        fp_superimposition_details.write('\t' + child_loop_pdb + '\t' + parent_loop_pdb + '\t' + 'rmsd: ' + rmsd + ', align_len: ' + align_len + ', score: ' + score + ', zscore: ' + zscore)

        if input_index_type == 'fasta':
            fp_superimposition_details.write('\t' + child_loop + '\t' + parent_loop)

        fp_superimposition_details.write('\n')
    
    fp_superimposition_details.close()

def load_pdb_fasta_mapping_and_fasta_seq_dict(cluster_id, alignment_data):
    pdb_chain_dict = {}
    pdb_res_mapping_dict = {}
    fasta_seq_dict = {}
    for l1 in alignment_data[cluster_id]:
        pdb_chain, _ = str(l1).strip().split(':')
        
        pdb_id, chain_id = pdb_chain.strip().split('_')
        if pdb_id not in pdb_chain_dict:
            pdb_chain_dict[pdb_id] = []
        pdb_chain_dict[pdb_id].append(chain_id)

        if pdb_chain not in pdb_res_mapping_dict:
            pdb_res_mapping_dict[pdb_chain] = load_pdb_res_map(pdb_chain)

    for pdb_id in pdb_chain_dict:
        fasta_seq_dict.update(load_fasta_seq(pdb_id, pdb_chain_dict[pdb_id]))

    return pdb_res_mapping_dict, fasta_seq_dict

def get_inter_subfamily_parent_info(ordered_dependency_list, edges_among_all_components):
    parent_usage = {}

    for component in ordered_dependency_list:
        loop_list = []
        for (i, r1), parent in component:
            loop_list.append((i, r1))

        for (i, r1), (j, r2), _ in edges_among_all_components:
            if (i, r1) in loop_list and (j, r2) not in loop_list:
                if (i, r1) not in parent_usage:
                    parent_usage[(i, r1)] = 0
                parent_usage[(i, r1)] += 1

    return parent_usage

def delete_motif_from_pymol(pymol_load_info_motif):
    pdb_load_name, chain_load_name, target_load_name, display_load_name, _, _, _, _ = pymol_load_info_motif
    pymol.cmd.delete(display_load_name)
    pymol.cmd.delete(target_load_name)
    pymol.cmd.delete(chain_load_name)
    pymol.cmd.delete(pdb_load_name)

def add_rotation_version_prefix(rotation_version):
    if len(rotation_version) > 0:
        return rotation_version + '__'
    return ''

def add_rotation_version_suffix(rotation_version):
    if len(rotation_version) > 0:
        return '__' + rotation_version
    return ''

def split_subfamily_cumulative_count(ordered_dependency_list):
    splitted_subfamily_cumulative_count = []
    
    for component in ordered_dependency_list:
        current_list = []
        remaining_motifs = len(component)
        # start_ind = 0
        previous_val = 0
        while remaining_motifs > max_no_of_motifs_in_superimposition:
            next_split_size = max_no_of_motifs_in_superimposition

            if remaining_motifs < 2 * max_no_of_motifs_in_superimposition:
                next_split_size = (remaining_motifs + 1) / 2
            
            # end_ind = start_ind + next_split_size
            remaining_motifs -= next_split_size

            # current_list.append(component[start_ind:end_ind])
            current_list.append(previous_val + next_split_size)
            previous_val += next_split_size
            # start_ind = end_ind

        # current_list.append(component[start_ind:])
        current_list.append(previous_val + remaining_motifs)

        splitted_subfamily_cumulative_count.append(current_list)

    # for i, component in enumerate(ordered_dependency_list):
    #     print(str(len(component)) + ': '),
    #     print(','.join(map(lambda x: str(x), splitted_subfamily_cumulative_count[i])))

    return splitted_subfamily_cumulative_count

def set_adj_res_view(r, pymol_load_info):

    fname = os.path.join(neares_protein_data_dir, r + '.adj_info')
    
    if os.path.isfile(fname):

        pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, lp, load_name_id = pymol_load_info

        fp = open(fname)
        lines = fp.readlines()
        fp.close()

        ring_mode, ring_transparency = lines[3].strip().split(',')
        ring_mode, ring_transparency = int(ring_mode), float(ring_transparency)

        items = lines[2].strip().split(',')
        adj_chain_str = ''
        for item in items:
            if len(item) > 0:
                chain_id = item.strip().split(':')[0]
                adj_chain_str += ' or chain ' + chain_id

        if len(adj_chain_str) == 0:
            adj_chain_str = target_load_name
        else:
            adj_chain_str = "(" + target_load_name + adj_chain_str + ")"
 
        base_atoms = ["C1'", "C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"]
        Giovanni_Ciriello_2010 = ["P", "OP1", "OP2", "O5'", "C5'", "C4'"] 

        select_str = ''
        for i, atom in enumerate(base_atoms):
            if i > 0:
                select_str += " or" 
            select_str += " name " + atom
        
        pymol.cmd.select(display_load_name + '_temp', target_load_name + " and (" + select_str + ")")
        pymol.cmd.select(display_load_name + '_temp2', 'byres ' + pdb_load_name + ' within 5 of ' + display_load_name + '_temp')
        pymol.cmd.select(display_load_name + '_temp3', 'byres ' + pdb_load_name + ' within 10 of ' + display_load_name + '_temp')
        pymol.cmd.select(display_load_name + '_label', display_load_name + '_temp2' + ' and ' + adj_chain_str)
        pymol.cmd.select(display_load_name + '_zoom', display_load_name + '_temp3' + ' and ' + adj_chain_str)


        pymol.cmd.set(name='cartoon_ring_mode',value=ring_mode,quiet=1)
        pymol.cmd.set(name='cartoon_ring_transparency',value=ring_transparency,quiet=1)
        pymol.cmd.label(display_load_name + '_label' + " and (name C2' or name CA)", "'%s-%s' %(resn, resi)")
        pymol.cmd.zoom(display_load_name + '_zoom')

        pymol.cmd.deselect()
        pymol.cmd.show('cartoon', pdb_load_name)

        pymol.cmd.delete(display_load_name + '_zoom')
        pymol.cmd.delete(display_load_name + '_label')
        pymol.cmd.delete(display_load_name + '_temp3')
        pymol.cmd.delete(display_load_name + '_temp2')
        pymol.cmd.delete(display_load_name + '_temp')

        return True
    return False

def check_and_save_pymol_figure_of_a_loop_with_protein(r, pymol_load_info, image_fname, show_extended_loop, show_label, label_atom, display_color, align_color, superimposition_output_dir):
    # fname = os.path.join(neares_protein_data_dir, r + '.adj_info')
    
    # if os.path.isfile(fname):
    #     pdb_load_name, chain_load_name, target_load_name, display_load_name, centroid, pdb_chain, lp, load_name_id = pymol_load_info

    #     fp = open(fname)
    #     lines = fp.readlines()
    #     fp.close()

    #     # items = lines[1].strip().split(',')
    #     ring_mode, ring_transparency = lines[3].strip().split(',')
    #     ring_mode, ring_transparency = int(ring_mode), float(ring_transparency)

    #     items = lines[2].strip().split(',')
    #     adj_chain_str = ''
    #     for i, item in enumerate(items):
    #         chain_id = item.strip().split(':')[0]
            
    #         if i > 0:
    #             adj_chain_str += ' or'

    #         adj_chain_str += ' chain ' + chain_id

    #     # select_str = ''
    #     # for item in items:
    #     #     chain_id, regions = item.strip().split(':')
    #     #     regions = regions.strip().replace('_', ',')
    #     #     select_str += ' or (resi ' + regions + ' and chain ' + chain_id + ')'

    #     # pymol.cmd.select(display_load_name + '_zoom', chain_load_name + select_str)

    #     base_atoms = ["C1'", "C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"]

    #     select_str = ''
    #     for i, atom in enumerate(base_atoms):
    #         if i > 0:
    #             select_str += " or" 
    #         select_str += " name " + atom
        
    #     pymol.cmd.select(display_load_name + '_temp', target_load_name + " and (" + select_str + ")")
    #     pymol.cmd.select(display_load_name + '_temp2', 'byres all within 5 of ' + display_load_name + '_temp')
    #     pymol.cmd.select(display_load_name + '_temp3', 'byres all within 10 of ' + display_load_name + '_temp')
    #     pymol.cmd.select(display_load_name + '_zoom', display_load_name + '_temp3' + ' and (' + target_load_name + ' or ' + adj_chain_str + ')')
    #     pymol.cmd.select(display_load_name + '_label', display_load_name + '_temp2' + ' and (' + target_load_name + ' or ' + adj_chain_str + ')')


    #     pymol.cmd.set(name='cartoon_ring_mode',value=ring_mode,quiet=1)
    #     pymol.cmd.set(name='cartoon_ring_transparency',value=ring_transparency,quiet=1)
    #     pymol.cmd.label(display_load_name + '_label' + " and (name C2' or name CA)", "'%s-%s' %(resn, resi)")
    #     pymol.cmd.zoom(display_load_name + '_zoom')
        
    #     pymol.cmd.deselect()
    #     pymol.cmd.show('cartoon', pdb_load_name)

    #     pymol.cmd.delete(display_load_name + '_zoom')
    #     pymol.cmd.delete(display_load_name + '_label')
    #     pymol.cmd.delete(display_load_name + '_temp3')
    #     pymol.cmd.delete(display_load_name + '_temp2')
    #     pymol.cmd.delete(display_load_name + '_temp')
    if set_adj_res_view(r, pymol_load_info):
        image_dir_with_adj_info = os.path.join(superimposition_output_dir, 'rotated_loop_images_with_adj_info/')
        create_directory(image_dir_with_adj_info)

        image_fname = os.path.join(image_dir_with_adj_info, os.path.basename(image_fname))
        image_fname = '.'.join(image_fname.split('.')[:-1]) + '_' + r +'.png'
        session_fname = '.'.join(image_fname.split('.')[:-1]) + '.pse'
        pymol.cmd.png(image_fname, 1200, 1200, dpi=300, ray=1, quiet=1)
        pymol.cmd.sync()
        if save_pymol_session == True:
            pymol.cmd.save(session_fname)
            pymol.cmd.sync()
        pymol.cmd.hide()
        config_pymol_cartoon('red', show_label)

def generate_pymol_images(time_in_distance_calc, removable_text_file_list, partial_pdbx_dir, summary_dir, superimposition_output_dir, subfamily_details_dir, superimposition_details_dir, representative_dir, pymol_session_dir, current_rmsd_data_dict, alignment_data, pdb_organism_details, loop_type, set_view_manually, draw_figures, show_extended_loop, is_cif=True):

    if generate_similarity_graph_image:
        plt_fig = plt.figure(frameon = False)
        plt_fig = plt.figure()
        plt_fig.set_size_inches(9, 9)

    fp_align_len_threshold = None
    subfamily_cluster_fp = None
    show_label = False

    if set_view_manually == True:
        pymol.finish_launching()

    else:
        if draw_figures == True: # or generate_rotation_matrices:
            # pymol.finish_launching()
            pymol.finish_launching(['pymol', '-cqQ'])

            # If all cluster length is <= 2, show labels (this is supposed to happen for testing purpose clusters)
            show_label = True
            for cluster_id in current_rmsd_data_dict:
                _, cluster_pairwise_alignment_details = current_rmsd_data_dict[cluster_id]
                if len(cluster_pairwise_alignment_details) > 2:
                    show_label = False
                    break

            if save_pymol_session == True:
                create_directory(pymol_session_dir)

        # To save alignment length threshold
        familywise_align_len_threshold_fn = os.path.join(summary_dir, 'Familywise_Align_Length_Threshold.txt')
        fp_align_len_threshold = open(familywise_align_len_threshold_fn, 'w')

        # File to write components(subfamilies) as cluster
        subfamily_cluster_fname = os.path.join(superimposition_output_dir, 'subfamily_cluster.csv')
        subfamily_cluster_fp = open(subfamily_cluster_fname, 'w')

    initial_loop_image_dir = os.path.join(superimposition_output_dir, 'initial_loop_images/')
    rotated_loop_image_dir = os.path.join(superimposition_output_dir, 'rotated_loop_images/')
    
    if draw_input_images == True:
        create_directory(initial_loop_image_dir)
    
    create_directory(rotated_loop_image_dir)
    create_directory(representative_dir)

    cluster_image_file_list = {}
    rmsd_and_alignment_summary_dict = {}
    subfamily_details_table_data = {}
    pdb_res_map_dict = {}

    representative_pymol_load_info = {}
    subfamily_pymol_load_info = {}
    all_pymol_load_info = {}

    for cluster_id in sorted(current_rmsd_data_dict):

        # logger.info('Started processing ' + cluster_id)
        time_start_for_cluster = time.time()

        representative_pymol_load_info[cluster_id] = {}
        subfamily_pymol_load_info[cluster_id] = {}
        cluster_image_file_list[cluster_id] = {}

        logger.info('Loading pdb-fasta seq dict')
        pdb_res_mapping_dict, fasta_seq_dict = load_pdb_fasta_mapping_and_fasta_seq_dict(cluster_id, alignment_data)
        _, cluster_pairwise_alignment_details = current_rmsd_data_dict[cluster_id]

        load_id_dict = {} # Store the id assigned for each loop
        loop_list = cluster_pairwise_alignment_details.keys()
        for loop_id, (i, r1) in enumerate(loop_list):
            load_id = str(cluster_id) + '__' + str(loop_id + 1)
            load_id_dict[(i, r1)] = (load_id, loop_id + 1)

        logger.info('Generating loop boundaries')
        ### Generate boundary of each loop in the context of the cluster (family)
        # loop_boundary_dict, loop_boundary_original_dict = generate_loop_boundary(cluster_pairwise_alignment_details, alignment_data[cluster_id.strip().split('_')[0]])
        loop_boundary_dict, loop_boundary_original_dict = generate_loop_boundary(cluster_pairwise_alignment_details, alignment_data[cluster_id])
        loop_boundary_dict = get_loop_boundary_pdb_index(loop_boundary_dict)
        loop_boundary_original_dict = get_loop_boundary_pdb_index(loop_boundary_original_dict)
        # loop_coord_dict = get_translated_coordinates(loop_boundary_dict, pdb_dir, is_cif)

        logger.info('Generating subfamilies for ' + cluster_id)
        ### Generate the optimal ordering of loops to show the best possible superimposition
        # To Do: Return the components of the graph
        # merged_components_features -> (central_node, cycle_nodes_of_the_component, component_nodes, component_directed_adjacency_list) of components
        time_start_for_dependency_list = time.time()
        ordered_dependency_list, merged_components_features, edges_among_all_components, edges_in_merged_components = generate_loop_print_dependency_v2(cluster_id, cluster_pairwise_alignment_details, alignment_data, fp_align_len_threshold)
        logger.info('Completed generating subfamilies. Time taken: ' + str(round(time.time() - time_start_for_dependency_list, 3)) + ' seconds.' )
        print('')
        
        splitted_subfamily_cumulative_count = split_subfamily_cumulative_count(ordered_dependency_list)

        parent_usage = get_inter_subfamily_parent_info(ordered_dependency_list, edges_among_all_components)
        
        # set view for the first loop of traversal
        if set_view_manually == True:
            get_view_from_user(partial_pdbx_dir, cluster_id, ordered_dependency_list, loop_boundary_dict, loop_boundary_original_dict, show_extended_loop, is_cif)
            continue
        
        r1_view = None
        if draw_figures == True:
            r1_view = load_view_file(ordered_dependency_list)

        generate_componentwise_bp_annotation_files(subfamily_details_dir, cluster_id, ordered_dependency_list, load_id_dict)
        if output_env == 'local':
            generate_componentwise_analysis_files(superimposition_output_dir, cluster_id, ordered_dependency_list, load_id_dict, alignment_data, current_rmsd_data_dict) 

        avg_rmsd, total_alignment_length = get_family_rmsd_and_alignment_summary(ordered_dependency_list, cluster_pairwise_alignment_details)
        rmsd_and_alignment_summary_dict[cluster_id] = (avg_rmsd, total_alignment_length)

        subfamily_colors, subfamily_colors_tuple = get_sub_family_colors(len(ordered_dependency_list))

        # Generate subfamily details table data
        avg_rmsd_align_to, total_alignment_length_align_to = get_align_to_rmsd_info(cluster_id, cluster_pairwise_alignment_details, alignment_data)
        avg_rmsd_sufamily_only, total_alignment_length_subfamily_only = get_family_rmsd_and_alignment_summary(ordered_dependency_list, cluster_pairwise_alignment_details, True)
        generate_subfamily_details_table_data(cluster_id, ordered_dependency_list, avg_rmsd_align_to, total_alignment_length_align_to, avg_rmsd, total_alignment_length, avg_rmsd_sufamily_only, total_alignment_length_subfamily_only, subfamily_details_table_data)
        
        if generate_similarity_graph_image:

            # Generate component graph image
            graph_image_list = []
            subfamily_dir = os.path.join(superimposition_output_dir, 'subfamily')
            graph_dir = os.path.join(subfamily_dir, 'graph')
            create_directory(graph_dir)

            # all_components_info_list = []
            for merged_component_id, a_merged_component_features in enumerate(merged_components_features):
                component_id = 1
                # merged_components_info_list = []
                total_loop = 0
                for component_features in a_merged_component_features:
                    image_fname = os.path.join(graph_dir,  str(cluster_id) + '__' + str(merged_component_id + 1) + '_' + str(component_id)+'.png')
                    graph_image_list.append((None, image_fname, len(component_features[2])))
                    component_info = get_component_graph(component_features, load_id_dict, subfamily_colors_tuple[merged_component_id]) 
                    # merged_components_info_list.append(component_info)
                    # all_components_info_list.append(component_info)

                    draw_graph([component_info], [], None, image_fname, plt_fig)
                    component_id += 1
                    total_loop += len(component_features[2])

            generate_subfamily_image(graph_image_list, pdb_organism_details, cluster_id, os.path.join(superimposition_output_dir, 'subfamily'), draw_figures, 'step2_subfamilies_identified', is_graph_image=True)

        ### Load all the loops in PDB and save corresponding images
        pymol_load_info_dict = {} # Store all pymol data load related info
        loop_display_info_dict = {} # Store only the info to show the loops in pymol

        if draw_figures:
            reset_pymol()        

        fp_representative = open(os.path.join(representative_dir, str(cluster_id) + "_representatives.txt"), "w")

        text_fname = os.path.join(superimposition_output_dir, str(cluster_id) + '.txt')
        fp_text_fname = open(text_fname, 'w')
        removable_text_file_list.append(text_fname)
        
        prev_component_count = 0
        initial_loop_image_dict = {}
        rotated_loop_image_list_dict = {}
        component_image_file_list_dict = {}

        displayed_motifs = {}

        for component_id, component in enumerate(ordered_dependency_list):
            # write subcluster_id to file
            if cluster_id in known_motif_fullname:
                subfamily_cluster_fp.write(known_motif_fullname[cluster_id] + '-Sub' + str(component_id + 1))
            else:
                subfamily_cluster_fp.write(str(cluster_id) + '-Sub' + str(component_id + 1))

            ###### To Do: add code for deleting loops from pymol
            if draw_figures == True and whole_family_superimpose == False:
                for (i, r1) in displayed_motifs:
                    if (i, r1) not in parent_usage or parent_usage[(i, r1)] < 1:
                        delete_motif_from_pymol(displayed_motifs[(i, r1)])
            # if draw_figures == True:
            #     # pymol.cmd._do('hide all')
            #     pymol.cmd.hide()
            #     wait_for_certain_time_according_to_wait_factor(prev_component_count)
            #     pymol.cmd.sync()

            if len(component) > 0:
                _, parent = component[0]
                if parent != None:
                    if parent not in parent_usage:
                        logger.error('Parent not found.')
                        sys.exit()
                    parent_usage[parent] -= 1

            subfamily_pymol_load_info[cluster_id][component_id] = []
            # Show the list in the original loop_id order
            for (i, r1), parent in component:
                
                if input_index_type == 'pdb':
                    r1_pdb_ind = convert_a_loop_from_FASTA_to_PDB(r1)
                    subfamily_cluster_fp.write(',' + r1_pdb_ind)
                
                else:    
                    subfamily_cluster_fp.write(',' + r1)

                load_id, _ = load_id_dict[(i, r1)]
                image_fname = os.path.join(initial_loop_image_dir, load_id + '.png')
                initial_loop_image_dict[(i, r1)] = (r1, image_fname, 1)

                if draw_figures:
                    pymol_load_info_dict[(i, r1)] = load_one_pdb_in_pymol(partial_pdbx_dir, i, r1, loop_boundary_dict, loop_boundary_original_dict, load_id, is_cif, 'gray')
                    loop_display_info_dict[(i, r1)], pymol_load_info_dict[(i, r1)] = translate_and_show_single_loop(pymol_load_info_dict[(i, r1)], loop_boundary_dict[(i,r1)], loop_boundary_original_dict[(i, r1)], load_id, image_fname, show_extended_loop, show_label, 'gray', 'gray')

                    displayed_motifs[(i, r1)] = pymol_load_info_dict[(i, r1)]
                    subfamily_pymol_load_info[cluster_id][component_id].append((i, r1))

                    if draw_input_images == True:
                        # pymol.cmd._do('hide all')
                        pymol.cmd.hide()
                        pymol.cmd.sync()
            
            subfamily_cluster_fp.write('\n')

            represetative_i, representative_loop = generate_subfamily_representative(fp_representative, cluster_id, component_id, component, alignment_data, cluster_pairwise_alignment_details, pdb_res_map_dict, generate_align_length_threshold(cluster_pairwise_alignment_details))
            representative_pymol_load_info[cluster_id][component_id] = (represetative_i, representative_loop)

            rotation_matrices = get_multiple_orientation_rotation_matrices()
            ## Generate images for different orientation
            for v, rotation_matrix in enumerate(rotation_matrices):
                rotation_version = ''
                if len(rotation_matrices) > 1:
                    rotation_version = 'v' + str(v + 1)
                
                # if draw_figures == True:
                #     pymol.cmd.hide()
                #     pymol.cmd.sync()

                ### Rotate, group and show the subfamilies 
                # family_loop_id = 0
                if rotation_version not in rotated_loop_image_list_dict:
                    rotated_loop_image_list_dict[rotation_version] = {}
                if component_id not in rotated_loop_image_list_dict[rotation_version]:
                    rotated_loop_image_list_dict[rotation_version][component_id] = []

                if rotation_version not in cluster_image_file_list[cluster_id]:
                    cluster_image_file_list[cluster_id][rotation_version] = []
                if rotation_version not in component_image_file_list_dict:
                    component_image_file_list_dict[rotation_version] = []


                ## Save rotated motifs to generate side-by-side image
                # current_cumulative_index = 0
                for component_loop_id, ((i, r1), parent) in enumerate(component):
            
                    # family_loop_id += 1
                    load_name_id, _ = load_id_dict[(i, r1)]
                    # if draw_figures:
                        # pymol.cmd._do('hide all')
                        # pymol.cmd.hide()
                        # pymol.cmd.sync()
                        # pdb_load_name1, chain_load_name1, target_load_name1, display_load_name1, centroid1, pdb_chain1, lp1, load_name_id = pymol_load_info_dict[(i, r1)]
                    
                    ### superimposition subfamilies

                    
                    image_fname = os.path.join(rotated_loop_image_dir, add_rotation_version_prefix(rotation_version) + load_name_id + '__3.png')
                    
                    component_id_str = str(component_id + 1)
                    # if len(splitted_subfamily_cumulative_count[component_id]) > 1:
                    #     component_id_str += get_string_equivalent_index(current_cumulative_index)

                    file_name = os.path.join(superimposition_output_dir,  add_rotation_version_prefix(rotation_version) + str(cluster_id) + '__' + component_id_str + '_' + str(component_loop_id + 1))
                    adj_image_fname = file_name + '.png'
                    
                    # Draw the first loop of the first  component  independent of any other loops
                    if component_id == 0 and component_loop_id == 0:
                        
                        if draw_figures:
                            display_load_name, align_load_name, chain_load_name = loop_display_info_dict[(i,r1)]
                            if r1_view != None and v == 0:
                                pymol.cmd.set_view(r1_view.strip())
                            else:
                                rotate_first_loop(pymol_load_info_dict[(i, r1)], rotation_matrix)

                            if scanx_align_to_superimposition == False:
                                show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', subfamily_colors[component_id])
                                check_and_save_pymol_figure_of_a_loop_with_protein(r1, pymol_load_info_dict[(i, r1)], adj_image_fname, show_extended_loop, show_label, "C2'", 'gray', subfamily_colors[component_id], superimposition_output_dir)

                            else:
                                show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', 'red')

                        rotated_loop_image_list_dict[rotation_version][component_id].append((r1, image_fname, 1))
                        if v == 0:
                            write_dock_file_list(component_id + 1, component_loop_id + 1, i, r1, 0, 'None', 0, 0, 0, fp_text_fname, cluster_pairwise_alignment_details, pdb_organism_details, '', '')
                    
                    else:
                        mobile_loop = strToNode(r1)
                        # mobile loop has the best superimposition feature (e.g. rmsd, alignment length) with target loop
                        j, r2 = parent
                        target_loop = strToNode(r2)
                
                        (t1, t2, zscore, cr1, cr2, aln1, aln2, score) = alignment_data[cluster_id][mobile_loop][target_loop]
                        
                        # if output_env == 'local':
                        pdb1_pm, pdb2_pm, i1_pm, i2_pm = aln_residue_temp(pdb_res_mapping_dict, fasta_seq_dict, r1, r2, cr1, cr2, aln1, aln2, 0, len(aln1)-1, 0)
                        # else:
                        #     pdb1_pm, pdb2_pm, i1_pm, i2_pm = aln_residue(r1, r2, cr1, cr2, aln1, aln2, 0, len(aln1)-1, 0)

                        rmsd, align_len = get_rmsd_align_len(i, r1, j, r2, cluster_pairwise_alignment_details)

                        # rmsd_rank = get_rmsd_rank(rmsd, align_len, is_length_adjusted_score)
                        if draw_figures:
                            display_load_name, align_load_name, chain_load_name = loop_display_info_dict[(i,r1)]
                            # pymol_load_info_dict[(i, r1)], display_load_name, display_load_name_cur = superimposition_pair_of_loops(pymol_load_info_dict[(j, r2)], pymol_load_info_dict[(i, r1)], load_id, pdb2_pm, pdb1_pm, loop_boundary_original_dict[(i, r1)], pdb_dir, display_load_name, image_fname, is_cif, rmsd_rank, show_label, component_id, component_loop_id)     
                            rotate_loop(partial_pdbx_dir, pymol_load_info_dict[(j, r2)], pymol_load_info_dict[(i, r1)], pdb2_pm, pdb1_pm, loop_boundary_original_dict[(i, r1)])

                            if scanx_align_to_superimposition == False:
                                show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', subfamily_colors[component_id])
                                check_and_save_pymol_figure_of_a_loop_with_protein(r1, pymol_load_info_dict[(i, r1)], adj_image_fname, show_extended_loop, show_label, "C2'", 'gray', subfamily_colors[component_id], superimposition_output_dir)
                            else:
                                show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', 'green')

                        rotated_loop_image_list_dict[rotation_version][component_id].append((r1, image_fname, 1))
                        # subfamily_colors[component_id]

                        if v == 0:
                            write_dock_file_list(component_id + 1, component_loop_id + 1, i, r1, j, r2, align_len, zscore, score, fp_text_fname, cluster_pairwise_alignment_details, pdb_organism_details, t1, t2)

                    if draw_figures:
                        # pymol.cmd._do('hide all')
                        pymol.cmd.hide()
                        pymol.cmd.sync()

                time_in_distance_calc += generate_representative_loop_image(time_in_distance_calc, representative_dir, rotation_version, cluster_id, component_id, represetative_i, representative_loop, loop_display_info_dict, draw_figures, show_extended_loop, show_label)
                # time.sleep(.100)

                #### Draw subfamilies superimposed
                # component_image_file_list = []
                # if draw_figures == True:
                #     # pymol.cmd._do('hide all')
                #     pymol.cmd.hide()
                #     wait_for_certain_time_according_to_wait_factor(prev_component_count)
                #     pymol.cmd.sync()

                image_fname = None
                r1 = ''
                current_cumulative_index = 0
                previous_cumulative_value = 0
                ## Generate progressive images and subfamily superimposition
                for component_loop_id, ((i, r1), parent) in enumerate(component):

                    ### superimposition subfamilies
                    component_id_str = str(component_id + 1)
                    if len(splitted_subfamily_cumulative_count[component_id]) > 1:
                        component_id_str += get_string_equivalent_index(current_cumulative_index)

                    file_name = os.path.join(superimposition_output_dir,  add_rotation_version_prefix(rotation_version) + str(cluster_id) + '__' + component_id_str + '_' + str(component_loop_id + 1))
                    image_fname = file_name + '.png'

                    if draw_figures:
                        display_load_name, align_load_name, chain_load_name = loop_display_info_dict[(i,r1)]

                        if scanx_align_to_superimposition == True:
                            if component_id == 0 and component_loop_id == 0:
                                show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', 'red')
                            else:
                                show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', 'green')
                        else:
                            # print('showing and saving ' + r1)
                            show_and_save_pymol_fig_of_a_loop(chain_load_name, display_load_name, align_load_name, image_fname, show_extended_loop, show_label, "C2'", 'gray', subfamily_colors[component_id])

                    if component_loop_id + 1 == splitted_subfamily_cumulative_count[component_id][current_cumulative_index]:
                        number_of_motifs = splitted_subfamily_cumulative_count[component_id][current_cumulative_index] - previous_cumulative_value
                        component_image_file_list_dict[rotation_version].append((r1, image_fname, number_of_motifs))

                        if draw_figures:
                            wait_for_certain_files_to_be_generated([image_fname], True)
                            pymol.cmd.sync()

                            pymol.cmd.hide()
                            wait_for_certain_time_according_to_wait_factor(number_of_motifs)
                            pymol.cmd.sync()

                        previous_cumulative_value = splitted_subfamily_cumulative_count[component_id][current_cumulative_index]
                        current_cumulative_index += 1

                # if image_fname != None:
                #     number_of_motifs = splitted_subfamily_cumulative_count[component_id][current_cumulative_index] - previous_cumulative_value
                #     component_image_file_list_dict[rotation_version].append((r1, image_fname, number_of_motifs))

            if draw_figures == True:
                if save_pymol_session == True:
                    pymol.cmd.deselect()
                    # pymol.cmd._do('save ' + os.path.join(pymol_session_dir, str(cluster_id) + '-Sub' + str(component_id+1) + '.pse'))
                    pymol.cmd.save(os.path.join(pymol_session_dir, str(cluster_id) + '-Sub' + str(component_id+1) + '.pse'))
                    pymol.cmd.sync()

                # pymol.cmd._do('hide all')
                pymol.cmd.hide()
                wait_for_certain_time_according_to_wait_factor(len(component))
                pymol.cmd.sync()

            # prev_component_count = len(component)
        fp_text_fname.close()
        generate_formatted_superimposition_details(superimposition_details_dir, cluster_id, text_fname, pdb_organism_details)
        all_pymol_load_info[cluster_id] = pymol_load_info_dict

        for rotation_version in component_image_file_list_dict:
            if len(component_image_file_list_dict[rotation_version]) > 0:
                if whole_family_superimpose == True and len(component_image_file_list_dict[rotation_version]) > 1:
                    component_image_file_list_dict[rotation_version] = generate_and_add_family_image(superimposition_output_dir, cluster_id, component_image_file_list_dict[rotation_version], ordered_dependency_list, loop_display_info_dict, rotation_version, draw_figures, show_extended_loop)

                cluster_image_file_list[cluster_id][rotation_version] = component_image_file_list_dict[rotation_version]

        fp_representative.close()

        if draw_input_images == True:
            initial_loop_image_list = [] # Store the list of ungrouped loop
            for i, r1 in loop_list:
                initial_loop_image_list.append(initial_loop_image_dict[(i, r1)])

            # Collage of input motifs
            generate_subfamily_image(initial_loop_image_list, pdb_organism_details, cluster_id, os.path.join(superimposition_output_dir, 'subfamily'), draw_figures, 'step1_initial', show_pdb_info=True, show_image_caption=False)

        for rotation_version in rotated_loop_image_list_dict:
            combined_rotated_image_list = []
            for component_id in rotated_loop_image_list_dict[rotation_version]:
                component_collage_image_fname = '-Sub' + str(component_id + 1) + add_rotation_version_suffix(rotation_version)
                generate_subfamily_image(rotated_loop_image_list_dict[rotation_version][component_id], pdb_organism_details, cluster_id, subfamily_details_dir, draw_figures, component_collage_image_fname, show_image_caption=subfamily_side_by_side_image_caption)
                # generate_subfamily_image(rotated_loop_image_list_dict[rotation_version][component_id], pdb_organism_details, cluster_id, subfamily_details_dir, draw_figures, component_collage_image_fname, show_pdb_info=True, show_image_caption=True)
                combined_rotated_image_list += rotated_loop_image_list_dict[rotation_version][component_id]

            generate_subfamily_image(combined_rotated_image_list, pdb_organism_details, cluster_id, os.path.join(superimposition_output_dir, 'subfamily'), draw_figures, 'step3_loops_rotate_and_grouped' + add_rotation_version_suffix(rotation_version), show_image_caption=family_side_by_side_image_caption)

        time_diff_for_cluster = time.time() - time_start_for_cluster
        logger.info('Processed cluster: ' + cluster_id)
        logger.info('Total loops: ' + str(len(loop_list)))
        # logger.info('Time elapsed: ' + str(round(time_diff_for_cluster/60)) + ' minutes.')
        logger.info('Time elapsed: ' + str(round(time_diff_for_cluster, 3)) + ' seconds.')
    
    if set_view_manually == True:
        pymol.cmd.quit()
        return

    if output_env == 'local':
        generate_table(summary_dir, subfamily_details_table_data, loop_type, True)
    else:
        generate_table(summary_dir, subfamily_details_table_data, loop_type)

    subfamily_cluster_fp.close()
    fp_align_len_threshold.close()

    write_rmsd_and_alignment_summary(rmsd_and_alignment_summary_dict, current_rmsd_data_dict)

    for cluster_id in sorted(cluster_image_file_list):
        for rotation_version in sorted(cluster_image_file_list[cluster_id]):
            generate_subfamily_image(cluster_image_file_list[cluster_id][rotation_version], pdb_organism_details, cluster_id, os.path.join(superimposition_output_dir, 'subfamily'), draw_figures, 'step4_subfamily_superimposed' + add_rotation_version_suffix(rotation_version))
    
    subfamily_out_dir = os.path.join(superimposition_output_dir, 'subfamily')
    for cluster_id in sorted(cluster_image_file_list):
        for rotation_version in sorted(cluster_image_file_list[cluster_id]):
            image_file_list = []
            if draw_input_images == True:
                image_file_list.append(os.path.join(subfamily_out_dir, cluster_id + '_step1_initial.png'))
            image_file_list.append(os.path.join(subfamily_out_dir, cluster_id + '_step3_loops_rotate_and_grouped' + add_rotation_version_suffix(rotation_version) + '.png'))
            image_file_list.append(os.path.join(subfamily_out_dir, cluster_id + '_step4_subfamily_superimposed' + add_rotation_version_suffix(rotation_version) + '.png'))
            create_combined_subfamily_collage(os.path.join(subfamily_out_dir, 'Subfamily_Combined'), image_file_list, cluster_id, draw_figures, 'step5_combined_subfamily_output' + add_rotation_version_suffix(rotation_version))

    if draw_figures == True and save_pymol_session == True:
        optimize_pymol_session(pymol_session_dir, representative_dir, all_pymol_load_info, subfamily_pymol_load_info, representative_pymol_load_info)

    return time_in_distance_calc

def optimize_pymol_session(pymol_session_dir, representative_dir, all_pymol_load_info, subfamily_pymol_load_info, representative_pymol_load_info):
    logger.info('Optimizing PyMol sessions.')
    for cluster_id in subfamily_pymol_load_info:
        for component_id in subfamily_pymol_load_info[cluster_id]:
            pymol_session_name = str(cluster_id) + '-Sub' + str(component_id+1) + '.pse'
            pymol.cmd.load(os.path.join(pymol_session_dir, pymol_session_name))
            pymol.cmd.sync()

            for (i, r) in all_pymol_load_info[cluster_id]:
                if (i, r) not in subfamily_pymol_load_info[cluster_id][component_id]:
                    delete_motif_from_pymol(all_pymol_load_info[cluster_id][(i, r)])

            pymol.cmd.delete('t')
            pymol.cmd.sync()

            pymol.cmd.show('cartoon')
            pymol.cmd.sync()

            pymol.cmd.save(os.path.join(pymol_session_dir, pymol_session_name))
            pymol.cmd.sync()

    for cluster_id in representative_pymol_load_info:
        for component_id in representative_pymol_load_info[cluster_id]:
            pymol_session_name = str(cluster_id) + '-Sub' + str(component_id+1) + '_repr.pse'
            pymol.cmd.load(os.path.join(representative_dir, pymol_session_name))
            pymol.cmd.sync()

            for (i, r) in all_pymol_load_info[cluster_id]:
                if (i, r) != representative_pymol_load_info[cluster_id][component_id]:
                    delete_motif_from_pymol(all_pymol_load_info[cluster_id][(i, r)])

            pymol.cmd.delete('t')

            pymol.cmd.sync()
            pymol.cmd.save(os.path.join(representative_dir, pymol_session_name))
            pymol.cmd.sync()
            
