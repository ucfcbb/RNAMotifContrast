import argparse
import sys
import os
import glob
import logging
import time

from config import *
sys.path.append(lib_dir)
sys.path.append(scripts_dir)
from my_log import *
from validators import *
from utils import *
from prepare_loops import *
from alignment_generator import *
from superimposition_generator import *
from partial_pdb_generator import *
from image_helper import *
from validators import *


from cif import *

import collections

def get_residue_count(d):
    return sum(map(lambda x: len(d[x]), d))

def main():

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Prepare input for RNAMotifContrast')
    parser.add_argument('-i', nargs='?', default='sample1.in', const='sample1.in', help='Input file containing motifs')
    parser.add_argument('-o', nargs='?', default='', const='', help='Subdirectory inside Output directory to save the results')

    parser.add_argument('-p', nargs='?', default=False, const=True, help='Generate PyMOL images (prerequisite: installation of PyMOL)')
    parser.add_argument('-v', nargs='?', default=False, const=True, help='To set orientation in the PyMOL image of the first loop in superimposition')
    parser.add_argument('-e', nargs='?', default='5', const='5', help='How many residues to extend beyond loop boundary to generate the loop.cif file')
    parser.add_argument('-n', nargs='?', default=False, const=True, help='No extended residues in the PyMOL image')
    
    parser.add_argument('-a', nargs='?', default=False, const=True, help='Keep all cluster loops without filtering based on overall zscore and rmsd')

    parser.add_argument('-d', nargs='?', default='', const='', help='To provide alignment directory if alignments from external sources is used [Alignment format is described in README file]')

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()

    user_input_fname = args.i
    superimposition_output_subdir = args.o

    draw_figures = args.p
    set_view_manually = args.v
    loop_cif_extension = int(args.e)
    show_extended_loop = not args.n

    filter_cluster = not args.a

    alignment_dir_user = args.d

    if set_view_manually == True:
        draw_figures = False

    superimposition_output_dir = os.path.join(root_dir, 'output')
    partial_pdbx_dir = os.path.join(data_dir, 'pdbx_extracted_ext' + str(loop_cif_extension))
    graphs_and_pickles_dir = os.path.join(data_dir, 'graphs_and_pickles')

    global alignment_dir_auto

    if len(superimposition_output_subdir) > 0:
        superimposition_output_dir = os.path.join(superimposition_output_dir, superimposition_output_subdir)
        # alignment_dir_auto = os.path.join(alignment_dir_auto, superimposition_output_subdir)
        # pickles_dir = os.path.join(pickles_dir, superimposition_output_subdir)

    # summary_dir = os.path.join(superimposition_output_dir, 'summary')
    # representative_dir = os.path.join(superimposition_output_dir, 'representatives')
    # progressive_dir = os.path.join(superimposition_output_dir, 'progressive')
    # subfamilies_dir = os.path.join(superimposition_output_dir, 'subfamilies')
    # superimposition_details_dir = os.path.join(superimposition_output_dir, 'superimposition_details')
    # pymol_session_dir = os.path.join(superimposition_output_dir, 'pymol_sessions')
    # os.path.join(superimposition_output_dir, 'subfamilywise_bp_ann')

    summary_dir = os.path.join(superimposition_output_dir, 'summary')
    subfamilies_dir = os.path.join(superimposition_output_dir, 'subfamilies')
    progressive_dir = os.path.join(superimposition_output_dir, 'progressive')
    superimposition_details_dir = os.path.join(superimposition_output_dir, 'superimposition_details')
    representative_dir = os.path.join(superimposition_output_dir, 'subfamily_representatives')
    subfamily_details_dir = os.path.join(superimposition_output_dir, 'subfamily_details')

    pymol_session_dir = subfamily_details_dir
    subfamilywise_bp_ann_dir = subfamily_details_dir

    alignment_dir = alignment_dir_auto
    is_alignment_from_user = False
    if alignment_dir_user != '':
        alignment_dir = os.path.join(data_dir, alignment_dir_user)
        is_alignment_from_user = True

    if generate_loop_source_info == True and len(superimposition_output_subdir) > 0 and (superimposition_output_subdir == 'HL' or superimposition_output_subdir == 'IL'):    
        assign_cluster_source(os.path.join(cluster_source_dir, superimposition_output_subdir + "_r3d_diff_loops_cluster.csv"), "R3D")
        assign_cluster_source(os.path.join(cluster_source_dir, superimposition_output_subdir + "_denovo_diff_loops_cluster.csv"), "DeNovo")
        assign_cluster_source(os.path.join(cluster_source_dir, superimposition_output_subdir + "_common_loops_cluster.csv"), "Both")
    
    input_fname = os.path.join(data_dir, user_input_fname)
    input_fname_base = os.path.basename(input_fname)

    prepare_executables()
    validate_all(input_fname, draw_figures)
    create_required_directories(partial_pdbx_dir, alignment_dir, superimposition_output_dir, subfamily_details_dir, summary_dir, superimposition_details_dir, representative_dir, pymol_session_dir, graphs_and_pickles_dir, set_view_manually)

    print('')
    logger.info('Reading input from ' + input_fname[base_path_len:])
    print('')

    families = {}
    fp_input = open(input_fname)
    loop_list = csv_to_list(fp_input.readlines())
    fp_input.close()

    loop_count = 0
    for item in loop_list:
        if len(item) > 2:
            # families[item[0]] = map(lambda x: str(strToNode(x)), item[1:]) # item[1:]
            families[item[0]] = item[1:]

    prepare_data(families)
    if input_index_type == 'pdb':
        families = convert_a_cluster_from_PDB_to_FASTA(families)
        for family_id in families:
            families[family_id] = map(lambda x: str(strToNode(x)), families[family_id])

    # for family_id in families:
    #     summ = 0.0
    #     print('\n'+family_id)
    #     loops = families[family_id]
    #     for loop_fasta in loops:
    #         loop_len = get_fasta_loop_length(loop_fasta)
    #         summ += loop_len
    #         print(loop_fasta + '\t' + str(loop_len))

    #     print('Average: ' + str(summ/len(loops)))

    # sys.exit()
            
    if remove_homolog_subfamily == True:
        # check if a sunfamily contains all homologs; keep one of them and make supercluster
        families = get_homolog_filtered_families(families)
    
    loop_count = 0
    loop_node_list_str = []
    for family_id in families:
        loops = families[family_id]
        loop_count += len(loops)
        for loop in loops:
            loop = str(strToNode(loop))
            loop_node_list_str.append(loop)
    
    duplicates = [item for item, count in collections.Counter(loop_node_list_str).items() if count > 1]
    if len(duplicates) > 0:
        print('duplicates:')
        print(duplicates)

    loop_node_list_str = sorted(list(set(loop_node_list_str)))

    logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in ' + str(len(families)) + ' famil' + ('ies' if len(families) > 1 else 'y') + '.')
    print('')

    # get pdb and fasta files 
    # and
    # generate loop.cif files
    # and annotation files
    
    previous_graph_file_reused = True
    all_alignment_files_found = False
    
    # prepare_data(families)
    prepare_loop_files(loop_node_list_str)    #chkd
    prepare_partial_pdbs(partial_pdbx_dir, loop_node_list_str, loop_cif_extension)
    # get_alignment_files(alignment_dir, families, loop_node_list_str, is_alignment_from_user)
    alignment_data = None
    rmsd_data_dict = None
    attempt = 1
    while attempt <= 5:
        attempt += 1
        previous_graph_file_reused, all_alignment_files_found = generate_best_alignment_data(input_fname_base, graphs_and_pickles_dir, alignment_dir, families)
        if all_alignment_files_found == False:
            get_alignment_files(alignment_dir, families, loop_node_list_str, is_alignment_from_user)
            continue
        alignment_data, rmsd_data_dict = load_alignment_and_rmsd_data(families, loop_node_list_str, input_fname_base, partial_pdbx_dir, alignment_dir, graphs_and_pickles_dir, previous_graph_file_reused)
        if alignment_data == None:
            delete_graph_file(input_fname_base, graphs_and_pickles_dir) # as alignments need to be generated, it seems existing graph file is invalid
            get_alignment_files(alignment_dir, families, loop_node_list_str, is_alignment_from_user)
            continue
        break

    directories = (partial_pdbx_dir, summary_dir, subfamilies_dir, subfamily_details_dir, representative_dir, superimposition_output_dir, superimposition_details_dir, progressive_dir, pymol_session_dir)

    removable_text_file_list = []
    # generate_superimposition_images(input_fname_base, removable_text_file_list, partial_pdbx_dir, alignment_dir, superimposition_output_dir, subfamily_details_dir, superimposition_output_subdir, summary_dir, subfamilies_dir, superimposition_details_dir, representative_dir, progressive_dir, graphs_and_pickles_dir, pymol_session_dir, user_input_fname, families, loop_node_list_str, previous_graph_file_reused, is_alignment_from_user, draw_figures, filter_cluster, set_view_manually, show_extended_loop)
    generate_superimposition_images(families, loop_node_list_str, alignment_data, rmsd_data_dict, draw_figures, filter_cluster, set_view_manually, show_extended_loop, is_alignment_from_user, user_input_fname, removable_text_file_list, directories, superimposition_output_subdir)
    if output_env == 'global':
        # time.sleep(5)
        # wait_for_certain_time_according_to_wait_factor(len(families))
        if draw_figures == True:
            cleanup_images(superimposition_output_dir, progressive_dir, representative_dir, subfamilies_dir)
        cleanup_output_directories(removable_text_file_list, superimposition_output_dir, representative_dir, progressive_dir, subfamilies_dir, superimposition_details_dir, pymol_session_dir, draw_figures)
        time_consuming_info_filelist = glob.glob(os.path.join(root_dir, 'log_time_consuming_aligns_job*.txt'))
        for file in time_consuming_info_filelist:
            os.remove(file)
        
    # generate_superimposition_images_using_alignto(superimposition_output_dir, partial_pdbx_dir, families, draw_figures)

    # fp = open(os.path.join(superimposition_output_dir, 'summary.txt'), 'w')
    # threshold_percentage = 75.0
    # bacteria_minor_subfamilies = []
    # global global_component_organism_stat
    # for cluster_id in global_component_organism_stat:
    #     print(cluster_id, len(global_component_organism_stat[cluster_id]))
    #     for component_id in global_component_organism_stat[cluster_id]:
    #         # print component_id
    #         total_count = 0
    #         bacteria_count = 0
    #         for org_type in global_component_organism_stat[cluster_id][component_id]:
    #             total_count += global_component_organism_stat[cluster_id][component_id][org_type]
    #             if org_type == 'Bacteria':
    #                 bacteria_count += global_component_organism_stat[cluster_id][component_id][org_type]
    #         bacteria_percentage = round(bacteria_count * 100.0 / total_count, 2)
    #         # print bacteria_percentage
    #         nonbacteria_percentage = 100.0 - bacteria_percentage
    #         if nonbacteria_percentage >= threshold_percentage:
    #             bacteria_minor_subfamilies.append((cluster_id, component_id, global_component_organism_stat[cluster_id][component_id]))

    # print(str(len(bacteria_minor_subfamilies)) + ' bacteria minor subfamillies found.')
    # fp.write(str(len(bacteria_minor_subfamilies)) + ' bacteria minor subfamillies found.\n')
    # for cluster_id, component_id, org_stat in bacteria_minor_subfamilies:
    #     print(cluster_id + '_' + str(component_id+1) + ':')
    #     print_a_dict_sorted(org_stat, separator=': ')
    #     fp.write(cluster_id + '_' + str(component_id+1) + ':\n')
    #     print_a_dict_sorted(org_stat, fp=fp, separator=': ')

    # fp.close()
    logger.info('\nTotal time taken: ' + str(round((time.time() - process_start_time), 3)) + ' seconds.\n')
        
if __name__ == '__main__':
    main()
