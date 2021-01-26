import argparse
import glob
import os
import sys
import pickle
sys.path.append('../../')
from util import *
from util_cluster import *
from residue import *

total_bp = 0
total_stk = 0



def generate_annotation_list(annotation, detailed):#, has_t_op):
    bp_ann_list = []
    stk_ann_list = []

    if detailed:
        bp_item_len = 5
        stk_item_len = 4
    else:
        bp_item_len = 4
        stk_item_len = 3

    for item in annotation:
        bp = "---"
        if len(item) == bp_item_len:

            if detailed:
                index1, index2, edges, orientation, bp = item
            else:
                index1, index2, edges, orientation = item

            i1 = min(index1, index2)
            i2 = max(index1, index2)

            if i1 != index1:
                edges, bp = get_inverse_bp_info(edges, bp)

            if detailed:
                item_to_append = (i1, i2, edges, orientation, bp)
            else:
                item_to_append = (i1, i2, edges, orientation)

            bp_ann_list.append(item_to_append)

        if len(item) == stk_item_len:
            if detailed:
                (index1, index2, direction, bp) = item
            else:
                (index1, index2, direction) = item

            i1 = min(index1, index2)
            i2 = max(index1, index2)

            if i1 != index1:
                direction, bp = get_inverse_stk_info(direction, bp)

            if detailed:
                item_to_append = (i1, i2, direction, bp)
            else:
                item_to_append = (i1, i2, direction)

            stk_ann_list.append(item_to_append)

    return bp_ann_list, stk_ann_list

def generate_stat(pdb_id, bp_ann_list, prefix, detailed=False):
    fp = open(prefix + "_residue_stat.txt", "a")
    fq = open(prefix + "_residue_edge_stat.txt", "a")
    fr = open(prefix + "_pairwise_edge_stat.txt", "a")

    residue_dict = {}
    residue_edge_dict = {}
    pairwise_edge_dict = {}

    for item in bp_ann_list:
        bp = "---"

        if detailed:
            (index1, index2, edges, orientation, bp) = item
        else:
            (index1, index2, edges, orientation) = item

        edge1 = edges[0]
        edge2 = edges[2]

        if (index1, bp[0]) not in residue_dict:
            residue_dict[(index1, bp[0])] = []
        residue_dict[(index1, bp[0])].append((bp, orientation[0] + edge1 + edge2)) # = residue_dict[index1] + 1
        if (index2, bp[2]) not in residue_dict:
            residue_dict[(index2, bp[2])] = []
        residue_dict[(index2, bp[2])].append((bp, orientation[0] + edge2 + edge1))    # = residue_dict[index2] + 1


        if (index1,bp[0],edge1) not in residue_edge_dict:
            residue_edge_dict[(index1,bp[0],edge1)] = []
        residue_edge_dict[(index1,bp[0],edge1)].append((bp[2], edge2))   # = residue_edge_dict[(index1,edge1)] + 1

        if (index2,bp[2],edge2) not in residue_edge_dict:
            residue_edge_dict[(index2,bp[2],edge2)] = []
        residue_edge_dict[(index2,bp[2],edge2)].append((bp[0], edge1))   # = residue_edge_dict[(index2,edge2)] + 1

        if (index1, index2, bp) not in pairwise_edge_dict:    
            pairwise_edge_dict[(index1, index2, bp)] = []
        pairwise_edge_dict[(index1, index2, bp)].append((orientation[0] + edge1 + edge2))

    for res_info in sorted(residue_dict):
        if len(residue_dict[res_info]) > 3:
            fp.write(pdb_id + "\t" + str(res_info) + "\t" + str(len(residue_dict[res_info])) + "\t" + str(residue_dict[res_info]) + "\n")
    
    for edge_info in sorted(residue_edge_dict):
        if len(residue_edge_dict[edge_info]) > 1:
            fq.write(pdb_id + "\t" + str(edge_info) + "\t" + str(len(residue_edge_dict[edge_info])) + "\t" + str(residue_edge_dict[edge_info]) + "\n")
    
    global edge_dict
    for pair_info in pairwise_edge_dict:
        if len(pairwise_edge_dict[pair_info]) > 1:
            # for edge in pairwise_edge_dict[pair_info]:
            #     if edge not in edge_dict:
            #         edge_dict[edge] = 0
            #     edge_dict[edge] = edge_dict[edge] + 1

            for ii, edge1 in enumerate(pairwise_edge_dict[pair_info]):
                if "." not in edge1:
                    if edge1 not in edge_dict:
                        edge_dict[edge1] = {}
                    for jj, edge2 in enumerate(pairwise_edge_dict[pair_info]):
                        if ii != jj:
                            if "." not in edge2:
                                if edge2 not in edge_dict[edge1]:
                                    edge_dict[edge1][edge2] = 0
                                edge_dict[edge1][edge2] = edge_dict[edge1][edge2] + 1

            fr.write(pdb_id + "\t" + str(pair_info) + "\t" + str(len(pairwise_edge_dict[pair_info])) + "\t" + str(sorted(pairwise_edge_dict[pair_info])) + "\n")

    fp.close()
    fq.close()
    fr.close()



def write_interaction_dict_in_file(category_dict, fn):
    # print_a_dict_sorted(category_dict)
    fp = open(fn, "w")
    for bp in sorted(category_dict):
        print bp + ": ",
        fp.write(bp + ": ")
        for interaction in sorted(category_dict[bp], key=lambda x: category_dict[bp][x], reverse=True):
            if category_dict[bp][interaction] > 0:
                print interaction + "-" + str(category_dict[bp][interaction]),
                fp.write(interaction + "-" + str(category_dict[bp][interaction]) + " ")
        print ""
        fp.write("\n")
    fp.close()

def write_bp_matrix(bp_dict, fn):
    # matrix_fn = fn[:-4] + "_matrix.txt"
    fx = open(fn, "w")
    for interaction in interaction_list:
        fx.write("\t" + interaction)
    fx.write("\n")
    for bp in bp_list:
        fx.write(bp)
        for interaction in interaction_list:
            fx.write("\t" + str(bp_dict[bp][interaction]))
        fx.write("\n")
    fx.close()

def resolve_pairwise_annotation_conflict_helper(pdb_id, ann_list, bp_interaction_category_rank, stk_interaction_category_rank, merged_annotation_dir, detailed):
    if detailed:
        bp_item_len = 5
        stk_item_len = 4
    else:
        print "Detailed BP info not available."
        sys.exit()

    pairwise_bp_dict = {}
    pairwise_stk_dict = {}

    print "Resolving " + pdb_id
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
            print "ERROR: invalid interact_info length in " + pdb_id
            sys.exit()

    pairwise_bp_dict = get_resolved_dict(pdb_id, pairwise_bp_dict, bp_interaction_category_rank, fp)
    pairwise_stk_dict = get_resolved_dict(pdb_id, pairwise_stk_dict, stk_interaction_category_rank, fp)

    resolved_annotation_list = []
    for (index1, index2) in pairwise_bp_dict:
        (bp, interaction) = pairwise_bp_dict[(index1, index2)][0]
        orientation = "cis" if interaction[0] == 'c' else "trans"
        edges = interaction[1] + "/" + interaction[2]
        resolved_annotation_list.append((index1, index2, edges, orientation, bp))

    for index1, index2 in pairwise_stk_dict:
        bp, direction = pairwise_stk_dict[(index1, index2)][0]
        resolved_annotation_list.append((index1, index2, direction, bp))

    resolved_annotation_list = sorted(resolved_annotation_list)
    # ann_dict[pdb_id] = resolved_annotation_list
    write_merged_annotation_to_file(pdb_id, resolved_annotation_list, merged_annotation_dir, detailed)

def resolve_pairwise_annotation_conflict(ann_dict, bp_interaction_category_rank, stk_interaction_category_rank, merged_annotation_dir, detailed):
    for pdb_id in ann_dict:
        resolve_pairwise_annotation_conflict_helper(pdb_id, ann_dict[pdb_id], bp_interaction_category_rank, stk_interaction_category_rank, merged_annotation_dir, detailed)
        
def get_resolved_dict(pdb_id, pairwise_dict, interaction_category_rank, fp):
    resolved_dict = {}
    for ind_pair in pairwise_dict:
        interactions = get_highest_freq_interactions(pairwise_dict[ind_pair])
        if len(interactions) > 1:
            best_rank = 1000
            best_rank_interaction = ("", "")
            for bp, interact in interactions:
                if bp not in interaction_category_rank:
                    print "ERROR: " + pdb_id + " te Ghapla ache."
                    print bp
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
            # print pairwise_dict[ind_pair]
            # print interactions
            # sys.exit()
        resolved_dict[ind_pair] = interactions

    return resolved_dict

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

def write_merged_annotation_to_file(pdb_id, ann_list, merged_annotation_dir, detailed):
    if detailed:
        bp_item_len = 5
        stk_item_len = 4
    else:
        print "Detailed BP info not available."
        sys.exit()

    # for pdb_id in ann_dict:
    fp = open(os.path.join(merged_annotation_dir, pdb_id + ".merged"), "w")
    for annotation in ann_list:
        if len(annotation) == bp_item_len:
            index1, index2, edges, orientation, bp = annotation
            fp.write(str(index1) + "\t" + str(index2) + "\t" + bp + "\t" + edges + "\t" + orientation + "\n")
        elif len(annotation) == stk_item_len:
            index1, index2, direction, bp = annotation
            fp.write(str(index1) + "\t" + str(index2) + "\t" + bp + "\t" + direction + "\n")
        else:
            print "ERROR: invalid interact_info length in " + pdb_id
            sys.exit()
    fp.close()

def initialize_interaction_dict(bp_category_dict, stk_category_dict):
    global bp_list
    global interaction_list

    res_list = ['A', 'C', 'G', 'U']
    edge_list = ['W', 'H', 'S']
    orientation_list = ['c', 't']
    stack_direction_list = ["upward", "downward", "inward", "outward"]

    for res1 in res_list:
        for res2 in res_list:
            bp_list.append(res1 + "-" + res2)

    for orientation in orientation_list:
        for edge1 in edge_list:
            for edge2 in edge_list:
                interaction_list.append(orientation + edge1 + edge2)

    for bp in bp_list:
        bp_category_dict[bp] = {}
        stk_category_dict[bp] = {}
        for interaction in interaction_list:
            bp_category_dict[bp][interaction] = 0
        for direction in stack_direction_list:
            stk_category_dict[bp][direction] = 0

    # print_a_dict_sorted(bp_category_dict)
    # print_a_dict_sorted(stk_category_dict)
    # sys.exit()

def initialize_log_files(prefix):#, loaded_t_op_index_dict):
    # if prefix == "fr3d" and len(loaded_t_op_index_dict) == 0:
    #     fr = open("fr3d_t_op_stat.txt", "w")
    #     fr.close()
    # fp = open(prefix + "_residue_stat.txt", "w")
    # fp.close()
    # fp = open(prefix + "_residue_edge_stat.txt", "w")
    # fp.close()
    fp = open(prefix + "_pairwise_edge_stat.txt", "w")
    fp.close()

    # if prefix == "merged":
    #     fp = open(prefix + "_residue_name_inconsistency.txt", "w")
    #     fp.close()

# ann_dict = {3:{}, 4:{}, 5:{}, 6:{}, 7:{}}
ann_dict = {}
# ann_category = {}
# ann_special_category = {}

# edge_dict = {}
bp_category_dict = {}
stk_category_dict = {}
total_bp = 0
total_stk = 0
bp_list = []
interaction_list = []

def main():

    # global ann_dict
    # global special_ann_dict
    # global ann_category
    # global ann_special_category

    # for i in range(3, 8):
    #     print_a_dict(ann_dict[i])

    # for item in special_ann_dict:
    #     print item
    #     item_cnt = min(len(special_ann_dict[item]), 10)
    #     print_a_list(special_ann_dict[item][:item_cnt])

    # print_a_dict(ann_category, "\t")

    # print_a_dict(ann_special_category, "\t")

    # prefix = "dssr"
    # prefix = "fr3d"
    prefix = "merged"

    detailed_info = True

    initialize_log_files(prefix)#, loaded_t_op_index_dict)

    dssr_dir = "../../Annotation/annDSSRnonPair/"
    fr3d_dir = "../../Annotation/fr3d_annotation/"
    merged_annotation_dir = "../../Annotation/dssr_fr3d_merged/"

    global ann_dict

    total_merged_interact = 0
    total_dssr_interact = 0
    total_fr3d_interact = 0

    initialize_interaction_dict(bp_category_dict, stk_category_dict)
    # test_list = ["4V9F.cif.out", "4WF9.cif.out"]
    # for fn in test_list:
    for fn in glob.glob(os.path.join(dssr_dir, "*.cif.out")):
        # fn = os.path.join(dssr_dir, "4V9F.cif.out")
        # fn = os.path.join(dssr_dir, "5K8H.cif.out")
        # fn = os.path.join(dssr_dir, "5AFI.cif.out")
        # fn = os.path.join(dssr_dir, "4R4V.cif.out")
        # fn = os.path.join(dssr_dir, "4WF9.cif.out")
        # fn = os.path.join(dssr_dir, fn)
        
        print "parsing %s" % fn
        pdb_id = os.path.basename(fn)[:-8]
        ann_dict[pdb_id] = parseDSSR(dssr_dir, fn, detailed_info)
        total_dssr_interact += len(ann_dict[pdb_id])
        # bp_ann_list, stk_ann_list = generate_annotation_list(ann_dict[pdb_id], detailed_info)

        if prefix != "dssr":
            fn = os.path.join(fr3d_dir, pdb_id + ".fr3d")
            if os.path.isfile(fn):

                fr3d_ann_list = parseFR3D(fr3d_dir, fn, detailed_info)
                total_fr3d_interact += len(fr3d_ann_list)
                # bp_fr3d_ann_list, stk_fr3d_ann_list = generate_annotation_list(fr3d_ann_dict, detailed_info)

                if prefix == "fr3d":
                    # bp_ann_list = bp_fr3d_ann_list
                    # stk_ann_list = stk_fr3d_ann_list
                    ann_dict[pdb_id] = fr3d_ann_list
                else:
                    # bp_ann_list.extend(bp_fr3d_ann_list)
                    # stk_ann_list.extend(stk_fr3d_ann_list)
                    ann_dict[pdb_id].extend(fr3d_ann_list)

            else:
                print "FR3D annotation not found for " + pdb_id

        
        # ann_dict[pdb_id] = list(set(ann_dict[pdb_id]))
        # ann_dict[pdb_id] = sorted(ann_dict[pdb_id])

        generate_category_stat(pdb_id, ann_dict[pdb_id], prefix, detailed_info)
        total_merged_interact += len(ann_dict[pdb_id])
        
        # break


    # if prefix == "fr3d":
    #     # print "Total PDB\t" + str(len(glob.glob(os.path.join(fr3d_dir, "*.fr3d"))))
    #     total_bp /= 2
    #     total_stk /= 2
    # else:
    print "Total PDB\t" + str(len(glob.glob(os.path.join(dssr_dir, "*.cif.out"))))
    # print "Total bp\t" + str(total_bp)
    # print "Total stack\t" + str(total_stk)
    print "Total DSSR interaction\t" + str(total_dssr_interact)
    print "Total FR3D interaction\t" + str(total_fr3d_interact)
    print "Total merged interaction\t" + str(total_merged_interact)
    print "Total merged (bp, stack)\t" + str(total_bp) + "\t" + str(total_stk)
    print "\n"
    # for e in edge_dict:
    #     print e
    #     print_a_dict(edge_dict[e])
    #     print ""
    write_interaction_dict_in_file(bp_category_dict, "bp_interaction_stat.txt")
    print ""
    write_interaction_dict_in_file(stk_category_dict, "stk_interaction_stat.txt")
    write_bp_matrix(bp_category_dict, "bp_interaction_stat_matrix.txt")
    bp_interaction_category_rank = generate_ranking_from_interaction_dict(bp_category_dict)
    stk_interaction_category_rank = generate_ranking_from_interaction_dict(stk_category_dict)

    resolve_pairwise_annotation_conflict(ann_dict, bp_interaction_category_rank, stk_interaction_category_rank, merged_annotation_dir, detailed_info)

if __name__ == "__main__":
    main()