import sys
import os
import platform
import logging
import shutil
import time
import glob
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.PDB import *

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from my_log import *
from classes import *

def prepare_executables():
    if output_env == 'local':
        os.chdir(dssr_dir)
        os.system('chmod +x x3dna-dssr')
        os.chdir(root_dir)

    scanx_aln_executable_fn = 'align_ga'
    if platform.system() == 'Darwin':
        scanx_aln_executable_fn = 'align_ga.mac'

    os.chdir(os.path.join(motifscanx_dir, 'bin'))
    os.system('chmod +x ' + scanx_aln_executable_fn)
    os.chdir(root_dir)

def wait_for_certain_time_according_to_wait_factor(n):
    wait_time = n * wait_factor
    wait_time = min(wait_time, max_wait_time)
    time.sleep(wait_time)
    # logger.info('waiting')

def is_all_files_generated(file_list):
    for file in file_list:
        if not os.path.isfile(file):
            return False
    return True

def wait_for_certain_files_to_be_generated(file_list, early_terminate=True):
    cur_wait_time = 0.0
    phase_wait = 10 * wait_factor
    while not is_all_files_generated(file_list):

        if early_terminate == True and cur_wait_time > wait_time:
            return False
        elif early_terminate == False and cur_wait_time > max_wait_time:
            return False

        logger.info('Waiting for files to be generated.')
        time.sleep(phase_wait)
        cur_wait_time += phase_wait
        # logger.info('waiting')

    return True

def create_directory(dir_to_create):
    if not os.path.exists(dir_to_create):
        os.makedirs(dir_to_create)

def delete_directory(dir_to_delete):
    if os.path.exists(dir_to_delete):
        shutil.rmtree(dir_to_delete)

def remove_all_from_dir(mypath, exception_dirs=[]):
    for root, dirs, files in os.walk(mypath):
        # for file in files:  #files in all dirs
        #     if file not in exception_files:
        #         os.remove(os.path.join(root, file))
        for dir in dirs:
            if dir not in exception_dirs:
                shutil.rmtree(os.path.join(root, dir))

def zscore(x, m, std):
    if std == 0:
        return 0
    return (x-m)/std

def isClose(a, b, precision):
    if abs(a-b) <= precision:
        return True
    return False

def csv_to_list(lines):
    list_of_lists = []
    for line in lines:
        pieces = line.strip().split(',')
        list_of_lists.append(list(map(lambda x: x.strip(), pieces)))

    return list_of_lists

def base_abbreviation(fn):
    """get the abbreviation of the modified residues from the 3DNA baselist"""
    ret = {}
    # fn = os.path.join(os.path.dirname(os.path.abspath( __file__ )), fn)
    fp = open(fn)
    for line in fp.readlines():
        line = line.rstrip()
        if line.startswith("#") or len(line) == 0:
            continue
        else:
            three_letter = line[:3].strip()
            one_letter = line[8]
            if one_letter == "T" or one_letter == "t":
                one_letter = "U"
            ret[three_letter] = one_letter.upper()
    fp.close()
    return ret

def amino_acid_collection(fn):
    """get the list of the amino acids from the file"""
    ret = {}
    # fn = os.path.join(os.path.dirname(os.path.abspath( __file__ )), fn)
    fp = open(fn)
    for line in fp.readlines():
        line = line.rstrip()
        if line.startswith("#") or len(line) == 0:
            continue
        else:
            three_letter = line[:3].strip().upper()
            one_letter = line[4]
            ret[three_letter] = one_letter.upper()
    fp.close()
    return ret

def get_aln_mapping(aln_seq1, aln_seq2):
    """
    :return ret1: dict[i]->j  i in seq1; j in seq2
    :return ret2: dict[j]->i  j in seq2; i in seq1
    """
    if len(aln_seq1) != len(aln_seq2):
        return None

    i = j = 0

    ret1 = {}
    ret2 = {}
    for k in range(len(aln_seq1)):
        if aln_seq1[k] == "-" and aln_seq2[k] != "-":
            j += 1
        elif aln_seq2[k] == "-" and aln_seq1[k] != "-":
            i += 1
        elif aln_seq1[k] != "-" and aln_seq2[k] != "-":
            ret1[i] = j
            ret2[j] = i
            i += 1
            j += 1

    return ret1, ret2

# def download_single_pdbx_or_fasta_file(fname, file_ext):
#     if file_ext == 'cif':
#         urlretrieve(pdbx_url + fname, os.path.join(pdbx_dir, fname))
#         fp = open(os.path.join(pdbx_dir, fname))
#         lines = fp.readlines()
#         fp.close()
#         if lines[0].startswith('data_' + fname.strip().split('.')[0]):
#             status = True

#     elif file_ext == 'fasta':
#         urlretrieve(fasta_url + fname.strip().split('.')[0], os.path.join(fasta_dir, fname))
#         fp = open(os.path.join(fasta_dir, fname))
#         lines = fp.readlines()
#         fp.close()
#         if lines[0].startswith('>' + fname.strip().split('.')[0]):
#             status = True
#     else:
#         logger.error('Invalid file_type provided to download.')

#     if status == False:
#         logger.error(fname + ' download unsuccessful.')

def getDSSRseqnum(input, with_res_name=False):
    seqnum = ''
    icode = ''
    resname = ''
    if '/' in input:
        temp = input.strip().split('/')
        resname = temp[0]
        seqnum = temp[1]
    elif len(input) > 2 and input[:2] in ['DA', 'DT', 'DU', 'DG', 'DC']:
        resname = input[:2]
        seqnum = input[2:]
    elif len(input) > 3 and input[2].isalpha():
        resname = input[:3]
        seqnum = input[3:]
    else:
        resname = input[0]
        seqnum = input[1:]
    if '^' in seqnum:
        temp = seqnum.strip().split('^')
        seqnum = temp[0]
        icode = temp[1]
        # fp = open("modified_residues.dat", "rb")
        # modified_residues = pickle.load(fp)
        # fp.close()
        # if pdb_id in m.keys() and m[pdb_id]:
        #     print "got it"
    if with_res_name:
        return resname, seqnum, icode
    return seqnum, icode

def load_fasta_seq(pdb_id, chains):
    fasta_seq_dict = {}
    fasta_fn = os.path.join(fasta_dir, pdb_id + '.fasta')
    for record in SeqIO.parse(fasta_fn, 'fasta'):
        # fasta_seq_dict[record.id.strip().split('|')[0].strip().split(':')[1]] = str(record.seq)
        # chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1].strip().split(',')
        # for chain_id in chain_ids:
        #     fasta_seq_dict[chain_id] = str(record.seq)
        chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1:]
        for chain_id in chain_ids:
            chain_id = chain_id.strip().strip(',')
            if '[' in chain_id:
                continue
                # chain_id = chain_id.split('[')[0].strip()
            elif ']' in chain_id:
                chain_id = chain_id.split(']')[0].strip()
            
            fasta_seq_dict[chain_id] = str(record.seq)

    return fasta_seq_dict
    
def get_loop_type(loop):
    _, regions = loop.strip().split(':')
    regions = regions.strip().split('_')
    region_cnt = len(regions)
    if region_cnt == 1:
        return 'HL'
    elif region_cnt == 2:
        return 'IL'
    elif region_cnt > 2:
        return 'ML'
    logger.error('Invalid loop')
    return ''

def strToNode(loop_str):

    (chain, region) = loop_str.split(':')
    segments = region.split('_')
    node = Node(chain, segments)

    return node

def get_local_alignment_index(loop_index, alignment_region):
    if len(alignment_region.strip()) == 0:
        return loop_index
    parts = loop_index.split(':')[1].strip().split('_')
    indices = []
    for part in parts:
        s, e = part.split('-')
        for i in range(int(s), int(e)+1):
            indices.append(i)
    converted_index = loop_index.split(':')[0] + ':'
    parts = alignment_region.split(',')
    for part in parts:
        s,e = part.split('-')
        if converted_index[-1] != ':':
            converted_index += '_'
        converted_index += str(indices[int(s)]) + '-' + str(indices[int(e)])
    return converted_index

def parse_scanx_alignment_block(lines, line_index):
    r1 = lines[line_index].split('::')[1].split(' and ')[0].strip().strip(':')
    r2 = lines[line_index].split('::')[1].split(' and ')[1].strip().strip(':')

    score_text = lines[line_index+1].split(':')[1].strip()
    if score_text == '':
        # score = -50.
        logger.error('ERROR: No alignment score found for: ' + r1 + ' and ' + r2)
        sys.exit()
    else:
        score = float(score_text)

    cr1 = get_local_alignment_index(r1, lines[line_index+3].split(':')[1].strip())
    cr2 = get_local_alignment_index(r2, lines[line_index+4].split(':')[1].strip())

    aln1 = lines[line_index+6].strip()
    aln2 = lines[line_index+7].strip()

    return r1, r2, cr1, cr2, aln1, aln2, score

def parse_scanx_alignment_block_raw(lines, line_index):
    r1 = lines[line_index].split('::')[1].split(' and ')[0].strip().strip(':')
    r2 = lines[line_index].split('::')[1].split(' and ')[1].strip().strip(':')

    score_text = lines[line_index+1].split(':')[1].strip()
    if score_text == '':
        # score = -50.
        logger.error('ERROR: No alignment score found for: ' + r1 + ' and ' + r2)
        sys.exit()
    else:
        score = float(score_text)

    cr1 = lines[line_index+3].split(':')[1].strip()
    cr2 = lines[line_index+4].split(':')[1].strip()

    aln1 = lines[line_index+6].strip()
    aln2 = lines[line_index+7].strip()

    matching_bp_info = []
    matching_stk_info = []
    i = line_index+10
    while not lines[i].startswith('#  Matched base-stacking interactions: '):
        matching_bp_info.append(list(map(lambda x: x.strip(), lines[i].strip().split('MATCHES'))))
        i += 1

    i += 1
    while not lines[i].startswith('Total Elapsed Time :'):
        matching_stk_info.append(list(map(lambda x: x.strip(), lines[i].strip().split('MATCHES'))))
        i += 1

    is_copied = False
    if lines[i].strip().endswith('(copied)'):
        is_copied = True
    elapsed_time = lines[i].strip().split(':')[1].strip().split(' ')[0].strip()

    return r1, r2, cr1, cr2, aln1, aln2, score, matching_bp_info, matching_stk_info, elapsed_time, is_copied, i

def get_loops_in_cluster(clusters):
    loops_in_cluster = []

    for c_id in clusters:
        for i in range(len(clusters[c_id])):
            node1 = strToNode(clusters[c_id][i])
            if node1 not in loops_in_cluster:
                loops_in_cluster.append(node1)

    return loops_in_cluster

def find_nodes_in_cluster(node, cluster_alignment_data):
    cid_nodelist_pair = []
    for c_id in cluster_alignment_data:
        if node in cluster_alignment_data[c_id]:
            cid_nodelist_pair.append((c_id, cluster_alignment_data[c_id][node]))
    
    return cid_nodelist_pair

def get_backbone_and_sugar_atoms():
    
    backbone_atoms = ["C3'", "C4'", "C5'", "O3'", "O5'", "P"]
    sugar_atoms = ["C1'", "C2'", "C3'", "C4'", "O4'"]
    # base_atoms = ["C1'", "C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"]

    return backbone_atoms, sugar_atoms

def get_z_scores(a_list, is_median=False):
    list_with_zscore = []
    mean = get_mean(a_list, is_median)
    std = np.std(a_list)
    for value in a_list:
        z_value = zscore(float(value), float(mean), float(std))
        list_with_zscore.append((value, z_value))
    return list_with_zscore

def get_fasta_loop_length(loop_fasta):
    pdb_chain, regions = loop_fasta.strip().split(':')
    return get_loop_length(regions.strip().split('_'))
    
def get_loop_length(segments):
    # will not work properly for pdb index
    loop_length = 0
    for segment in segments:
        pcs = segment.split("-")
        if pcs[0][-1].isalpha():
            s = int(pcs[0].strip().split(".")[0].strip())
        else:
            s = int(pcs[0].strip())
        if pcs[1][-1].isalpha():
            e = int(pcs[1].strip().split(".")[0].strip())
        else:
            e = int(pcs[1].strip())
        # s = int(pieces[0])
        # e = int(pieces[1])
        loop_length += (e-s+1)
    return loop_length

def get_zscore_rank(zscore):
    if zscore > 3.0:
        return 1
    elif zscore > 1.8:
        return 2
    elif zscore > 1.0:
        return 3
    elif zscore > 0.5:
        return 4
    elif zscore > 0.0:
        return 5
    else:
        return 100

def get_rmsd_rank(rmsd, align_length, is_length_adjusted_score):
    if is_length_adjusted_score:
        rmsd = rmsd * math.sqrt(align_length)

    if rmsd < 0.5:
        return 1
    elif rmsd < 1.0:
        return 2
    elif rmsd < 2.0:
        return 3
    elif rmsd < 4.0:
        return 4
    elif rmsd < 19.0:
        return 5
    else:
        return 100

def get_mean(a_list, is_median=False):
    if is_median:
        return round(np.median(a_list), 1)
    return round(np.mean(a_list), 1)

def print_a_dict_sorted(a_dict, fp=None, separator=": "):
    if fp == None:
        for key in sorted(a_dict):
            print(str(key) + separator + str(a_dict[key]))
        print('')
    else:
        for key in sorted(a_dict):
            fp.write(str(key) + separator + str(a_dict[key]) + "\n")
        fp.write("\n")

def print_a_list(a_list, fp=None):
    for item in a_list:
        if fp == None:
            print(item)
        else:
            fp.write(str(item) + "\n")
    if fp == None:
        print('')
    else:
        fp.write("\n")

def get_motif_family_short_code(family_name):
    if family_name.lower() in known_motif_shortcode:
        return known_motif_shortcode[family_name.lower()]
    else:
        return family_name

def rotate(l, x):
    return l[-x:] + l[:-x]

def get_all_loop_combination(loop):
    loop_combinations = []
    pdb_chain, regions = loop.strip().split(':')
    regions = regions.strip().split('_')
    loop = []
    for region in regions:
        s, e = region.strip().split('-')
        loop.append((s, e))
    
    for i in range(len(loop)):
        loop_i = rotate(loop, i)
        loop_combinations.append(pdb_chain + ':' + '_'.join(list(map(lambda x: '-'.join(x), loop_i))))

    # print(loop_combinations)
    return loop_combinations

def get_separated_index_icode(index):
    ind = icode = None
    if index[-1].isalpha():
        ind = int(index[:-1])
        icode = index[-1]
    else:
        ind = int(index)
        icode = ' '
    return ind, icode

def convert_a_cluster_from_FASTA_to_PDB(families):
    families_pdb = {}
    for family_id in families:
        families_pdb[family_id] = []
        loops = families[family_id]
        for loop in loops:
            loop_pdb = convert_a_loop_from_FASTA_to_PDB(loop)
            families_pdb[family_id].append(loop_pdb)
    return families_pdb

def convert_a_cluster_from_PDB_to_FASTA(families):
    families_pdb = {}
    for family_id in families:
        families_pdb[family_id] = []
        loops = families[family_id]
        for loop in loops:
            loop_pdb = convert_a_loop_from_PDB_to_FASTA(loop)
            families_pdb[family_id].append(loop_pdb)
    return families_pdb

def convert_a_loop_from_PDB_to_FASTA(loop):
    pdb_chain, segments = loop.strip().split(':')
    pdb_id, chain_id = pdb_chain.strip().split('_')

    mapping_file_name = pdb_chain + '.rmsx.nch'
    converter = PDB_FASTA_Index_Converter(pdb_fasta_mapping_dir, mapping_file_name)

    segments = segments.strip().split('_')
    converted_segments = []
    for segment in segments:
        a, b = segment.strip().split('-')
        icode_a = ''
        icode_b = ''
        if '.' in a:
            a, icode_a = a.strip().split('.')
        if '.' in b:
            b, icode_b = b.strip().split('.')

        a_pdb = Chainindex(chain_id, int(a), icode_a)
        b_pdb = Chainindex(chain_id, int(b), icode_b)

        a_fasta = converter.convert_PDBindx_To_FASTAindx(a_pdb)
        b_fasta = converter.convert_PDBindx_To_FASTAindx(b_pdb)

        converted_segments.append(str(a_fasta) + '-' + str(b_fasta))

    converted_segments = '_'.join(converted_segments)
    converted_loop = pdb_chain + ':' + converted_segments

    return converted_loop

def convert_a_loop_from_FASTA_to_PDB(loop):
    pdb_chain, segments = loop.strip().split(':')
    pdb_id, chain_id = pdb_chain.strip().split('_')

    mapping_file_name = pdb_chain + '.rmsx.nch'
    converter = PDB_FASTA_Index_Converter(pdb_fasta_mapping_dir, mapping_file_name)

    segments = segments.strip().split('_')
    converted_segments = []
    for segment in segments:
        a, b = segment.strip().split('-')

        a_fasta = a
        b_fasta = b

        a_pdb = converter.convert_FASTAindx_To_PDBindx(a_fasta)
        b_pdb = converter.convert_FASTAindx_To_PDBindx(b_fasta)

        if len(a_pdb.icode) == 0:
            a_pdb = str(a_pdb.seqnum)
        else:
            a_pdb = str(a_pdb.seqnum) + '.' + str(a_pdb.icode)

        if len(b_pdb.icode) == 0:
            b_pdb = str(b_pdb.seqnum)
        else:
            b_pdb = str(b_pdb.seqnum) + '.' + str(b_pdb.icode)

        converted_segments.append(a_pdb + '-' + b_pdb)
        
    converted_segments = '_'.join(converted_segments)
    converted_loop = pdb_chain + ':' + converted_segments

    return converted_loop

def get_loop_cluster_source(loop):
    if loop != None:
        pdb_chain = loop.strip().split(":")[0]
        if show_cluster_source == True:
            loop = strToNode(loop)
            if loop in loop_cluster_source:
                return loop_cluster_source[loop]
    return "N/A"

def assign_cluster_source(filename, source_str):
    global loop_cluster_source
    fp = open(filename)
    for line in fp.readlines():
        pieces = line.strip().split(",")
        for piece in pieces:
            if ":" in piece:
                loop_cluster_source[strToNode(piece.strip())] = source_str

def cleanup_output_directories(removable_text_file_list, superimposition_output_dir, representative_dir, progressive_dir, subfamilies_dir, superimposition_details_dir, pymol_session_dir, draw_figures):
    # remove_all_from_dir(superimposition_output_dir, ['best_alignment_graph', 'subfamilywise_bp_ann', 'pymol_sessions', 'representatives', 'subfamily'])
    # remove_all_from_dir(os.path.join(superimposition_output_dir, 'initial_loop_images'))
    # remove_all_from_dir(os.path.join(superimposition_output_dir, 'rotated_loop_images'))
    # print(os.walk(superimposition_output_dir))
    # text_filelist = glob.glob(os.path.join(superimposition_output_dir, '*.txt'))
    # for file in text_filelist:
    for file in removable_text_file_list:
        # if 'representatives' not in file:
        if os.path.isfile(file):
            os.remove(file)

    # time.sleep(5)
    wait_for_certain_time_according_to_wait_factor(1)
    # remove_all_from_dir(superimposition_output_dir, ['subfamilywise_bp_ann', os.path.basename(superimposition_details_dir), os.path.basename(pymol_session_dir), os.path.basename(representative_dir), os.path.basename(progressive_dir), os.path.basename(subfamilies_dir)])
    delete_directory(os.path.join(superimposition_output_dir, 'best_alignment_graph'))
    delete_directory(os.path.join(superimposition_output_dir, 'initial_loop_images'))
    delete_directory(os.path.join(superimposition_output_dir, 'rotated_loop_images'))
    delete_directory(os.path.join(superimposition_output_dir, 'subfamily'))
    delete_directory(temp_dir)
    
            # os.remove(file)
    # for root, dirs, files in os.walk(superimposition_output_dir):
    #     for file in files:
    #         if file.endswith('.txt') and 'representatives' not in file:
    #             os.remove(os.path.join(root, file))

def create_required_directories(partial_pdbx_dir, alignment_dir, superimposition_output_dir, subfamily_details_dir, summary_dir, superimposition_details_dir, representative_dir, pymol_session_dir, pickles_dir, set_view_manually):
    # create_directory(data_dir)
    # create_directory(views_dir)

    create_directory(pdbx_dir)
    create_directory(partial_pdbx_dir)
    create_directory(fasta_dir)
    create_directory(loop_dir)
    create_directory(pdb_fasta_mapping_dir)
    create_directory(alignment_dir)
    create_directory(annotation_dir)
    create_directory(pickles_dir)
    create_directory(subfamily_details_dir)

    if set_view_manually == False:
        create_directory(superimposition_output_dir)
        create_directory(summary_dir)
        create_directory(superimposition_details_dir)
        create_directory(representative_dir)
        if save_pymol_session == True:
            create_directory(pymol_session_dir)


    # if len(loop_type) > 0 or use_pickle_file == True:
    #     create_directory(os.path.join(data_dir, 'pickles', loop_type + '_Pickles'))

def get_string_equivalent_index(current_cumulative_index):
    id_str = ''
    current_cumulative_index += 1
    while current_cumulative_index:
        current_cumulative_index, mod_val = divmod(current_cumulative_index - 1, 26)
        id_str = chr(ord('a') + mod_val) + id_str

    return id_str

def get_loop_region_identifier_line(aligned_seq, regions, ext_len_list, delim):

    line = ''

    ch_ind_list = []
    for i in range(len(aligned_seq)):
        if aligned_seq[i] != '-':
            ch_ind_list.append(i)

    loop_indices = reduce(lambda y, z: y+z, map(lambda x: list(range(x[0], x[1]+1)), regions))
    # loop_indices_extended = reduce(lambda x, y: range(x[0]-extension_length, x[1]+1+extension_length)+range(y[0]-extension_length, y[1]+1+extension_length), regions)
    loop_indices_extended = []
    prev_ext_e = -1
    for i, (s, e) in enumerate(regions):
        ext_s = s - ext_len_list[i][0]
        ext_e = e + ext_len_list[i][1]

        if prev_ext_e != -1 and ext_s <= prev_ext_e:
            ext_s = prev_ext_e + 1

        loop_indices_extended += range(ext_s, ext_e+1)
        prev_ext_e = ext_e

    loop_indices = sorted(list(set(loop_indices)))
    loop_indices_extended = sorted(list(set(loop_indices_extended)))

    last_ind = 0
    in_loop_indices = False
    last_i = -1
    loop_character_indices = []
    for i, index in enumerate(loop_indices_extended):
        if in_loop_indices == False and index in loop_indices:
            in_loop_indices = True
            line += ' ' * (ch_ind_list[i] - last_ind)
            last_ind = ch_ind_list[i]
            last_i = i

        if in_loop_indices == True and index not in loop_indices:
            in_loop_indices = False
            line += delim * (ch_ind_list[i-1] - last_ind + 1)
            last_ind = ch_ind_list[i-1] + 1
            loop_character_indices += ch_ind_list[last_i:i]
            last_i = i

    if in_loop_indices == True:
        line += delim * (len(aligned_seq) - 1 - last_ind + 1)
        loop_character_indices += ch_ind_list[last_i:len(aligned_seq)]

    # print(loop_character_indices)
    return line, loop_character_indices

def update_homolog_data(homolog_set_list, l1, l2):
    set_found = False
    set_ind_for_l1 = -1
    set_ind_for_l2 = -1
    for i, homolog_set in enumerate(homolog_set_list):
        if l1 in homolog_set:
            set_ind_for_l1 = i

        if l2 in homolog_set:
            set_ind_for_l2 = i

        if set_ind_for_l1 > -1 and set_ind_for_l2 > -1:
            break

    if set_ind_for_l1 > -1 and set_ind_for_l2 > -1:
        if set_ind_for_l1 != set_ind_for_l2:
            set_for_l1 = homolog_set_list[set_ind_for_l1]
            set_for_l2 = homolog_set_list[set_ind_for_l2]

            homolog_set_list = [homolog_set for i, homolog_set in enumerate(homolog_set_list) if i != set_ind_for_l1 and i != set_ind_for_l2]
            homolog_set_list.append(set_for_l1 | set_for_l2)

    elif set_ind_for_l1 == -1 and set_ind_for_l2 == -1:
        new_homolog_set = set()
        new_homolog_set.add(l1)
        new_homolog_set.add(l2)
        homolog_set_list.append(new_homolog_set)

    else:
        if set_ind_for_l1 == -1:
            homolog_set_list[set_ind_for_l2].add(l1)
        else:
            homolog_set_list[set_ind_for_l1].add(l2)

def generate_sequence_alignment_for_all_pairs(family_id, loop_list, extension_length=50):

    homolog_set_list = []

    loop_node_list = map(lambda x: strToNode(x), loop_list)
    pdb_list = map(lambda x: x.strip().split(':')[0].strip().split('_')[0], loop_list)

    pdb_list = sorted(list(set(pdb_list)))
    loop_node_list_str = sorted(map(lambda x: str(x), list(set(loop_node_list))))

    ref_seq_dict = {}
    for pdb_id in pdb_list:
        fasta_fn = os.path.join(fasta_dir, pdb_id + '.fasta')
        if os.path.isfile(fasta_fn):
            ref_seq_dict[pdb_id] = {}
            for record in SeqIO.parse(fasta_fn, 'fasta'):
                # ref_seq_dict[pdb_id][record.id.strip().split('|')[0].strip().split(':')[1]] = str(record.seq)
                chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1].strip().split(',')
                for chain_id in chain_ids:
                    ref_seq_dict[chain_id] = str(record.seq)

    for i in range(len(loop_node_list_str)):
        l1 = loop_node_list_str[i]

        pdb_chain1, regions1 = l1.strip().split(':')
        pdb_id1, chain_id1 = pdb_chain1.strip().split('_')
        fasta_seq1 = ref_seq_dict[pdb_id1][chain_id1]

        regions1 = regions1.strip().split('_')
        regions1 = sorted(map(lambda x: tuple(map(lambda y: int(y), x.strip().split('-'))), regions1))
        
        ext_len_list1 = []
        for ii, (s, e) in enumerate(regions1):
            left_range = 0
            if ii > 0:
                left_range = s - (s - (regions1[ii-1][1] + 1) + 1) / 2

            right_range = len(fasta_seq1) - 1
            if ii < len(regions1) - 1:
                right_range = e + ((regions1[ii+1][0] - 1) - e + 1) / 2

            left_ext_length = min(extension_length, s - left_range)
            right_ext_length = min(extension_length, right_range - e)
            ext_len_list1.append((left_ext_length, right_ext_length))

        for j in range(i+1, len(loop_node_list_str)):
            l2 = loop_node_list_str[j]

            pdb_chain2, regions2 = l2.strip().split(':')
            pdb_id2, chain_id2 = pdb_chain2.strip().split('_')
            fasta_seq2 = ref_seq_dict[pdb_id2][chain_id2]

            regions2 = regions2.strip().split('_')
            regions2 = sorted(map(lambda x: tuple(map(lambda y: int(y), x.strip().split('-'))), regions2))
            
            ext_len_list2 = []
            for jj, (s, e) in enumerate(regions2):
                left_range = 0
                if jj > 0:
                    left_range = s - (s - (regions2[jj-1][1] + 1) + 1) / 2

                right_range = len(fasta_seq2) - 1
                if jj < len(regions2) - 1:
                    right_range = e + ((regions2[jj+1][0] - 1) - e + 1) / 2

                left_ext_length = min(extension_length, s - left_range)
                right_ext_length = min(extension_length, right_range - e)
                ext_len_list2.append((left_ext_length, right_ext_length))

            if len(ext_len_list1) != len(ext_len_list2):
                print('Comparing different type of loops is not feasible. Exitting.')
                sys.exit()

            ext_len_list = []
            for i in range(len(ext_len_list1)):
                left1, right1 = ext_len_list1[i]
                left2, right2 = ext_len_list2[i]
                left = min(left1, left2)
                right = min(right1, right2)
                ext_len_list.append((left, right))

            prev_e1 = -1
            prev_e2 = -1
            seq1 = ''
            seq2 = ''
            for i in range(len(ext_len_list)):
                left, right = ext_len_list[i]

                s1 = regions1[i][0] - left
                e1 = regions1[i][1] + right
                s2 = regions2[i][0] - left
                e2 = regions2[i][1] + right

                if prev_e1 != -1 and s1 <= prev_e1:
                    s1 = prev_e1 + 1
                if prev_e2 != -1 and s2 <= prev_e2:
                    s2 = prev_e2 + 1
                    
                seq1 += fasta_seq1[s1 : e1 + 1]
                seq2 += fasta_seq2[s2 : e2 + 1]

                prev_e1 = e1
                prev_e2 = e2

            # aln = pairwise2.align.globalms(seq1, seq2, 5, -3, -10, -1)
            aln = pairwise2.align.globalxx(seq1, seq2)
            # (aln_seq1, aln_seq2, _, _, _) = aln[0]
            aln = sorted(aln, key=lambda x: x[4])

            lines = pairwise2.format_alignment(*aln[0]).strip().split('\n')
            identifier_line1, loop1_character_indices = get_loop_region_identifier_line(lines[0], regions1, ext_len_list, 'v')
            identifier_line2, loop2_character_indices = get_loop_region_identifier_line(lines[2], regions2, ext_len_list, '^')
            # new_lines = identifier_line1 + '\n' + '\n'.join(lines[:3]) + '\n' + identifier_line2 + '\n' + lines[3] + '\n'
            
            # if len(loop_list) <= 50:
            #     fp_output = open(output_fname1, 'a')
            #     fp_output.write(str(family_id) + '\n')
            #     fp_output.write(l1 + ' and ' + l2 + '\n')
            #     fp_output.write('l1 seq length: ' + str(len(seq1)) + '\nl2 seq length: ' + str(len(seq2)) + '\n')
            #     # fp_output.write(pairwise2.format_alignment(*aln[0]))
            #     fp_output.write(new_lines)
            #     fp_output.write('\n')
            #     fp_output.close()

            score = aln[0][-3]
            min_seq_len = min(len(seq1), len(seq2))
            percentage1 = ((min_seq_len - score) * 100.0) / min_seq_len

            loop1_character_indices = set(loop1_character_indices)
            loop2_character_indices = set(loop2_character_indices)
            common_character_indices = loop1_character_indices & loop2_character_indices
            min_character_indices_len = min(len(loop1_character_indices), len(loop2_character_indices))
            percentage2 = (len(common_character_indices) * 100.0) / min_character_indices_len

            percentage2 = 100.0 - percentage2

            if percentage2 == 0:
                # fp_output = open(output_fname2, 'a')
                # fp_output.write(str(family_id) + '\n')
                # fp_output.write(l1 + ' and ' + l2 + '\n')
                # fp_output.write('l1 seq length: ' + str(len(seq1)) + '\nl2 seq length: ' + str(len(seq2)) + '\n')
                # # fp_output.write(pairwise2.format_alignment(*aln[0]))
                # fp_output.write(new_lines)
                # fp_output.write('\n\n')
                # fp_output.close()

                if percentage1 <= 5:
                    update_homolog_data(homolog_set_list, l1, l2)

            # elif percentage2 <= percentage_threshold_for_nearly_homologs:
            #     fp_output = open(output_fname3, 'a')
            #     fp_output.write(str(family_id) + '\n')
            #     fp_output.write(l1 + ' and ' + l2 + '\n')
            #     fp_output.write('l1 seq length: ' + str(len(seq1)) + '\nl2 seq length: ' + str(len(seq2)) + '\n')
            #     # fp_output.write(pairwise2.format_alignment(*aln[0]))
            #     fp_output.write(new_lines)
            #     fp_output.write('\n\n')
            #     fp_output.close()

                # if filter_nearly_homologs == True:
                #     update_homolog_data(homolog_set_list, l1, l2)

    return homolog_set_list

def get_homolog_filtered_families(families):
    logger.info('Filtering out the subfamilies with all homologs ...')
    
    new_families = {}
    for family_id in sorted(families):
        new_family_id = family_id.strip().split('-Sub')[0].strip()
        
        if new_family_id not in new_families:
            new_families[new_family_id] = []

        loop_list = families[family_id]
        print('Checking ' + family_id + '(' + str(len(loop_list)) + ' loops)')
        
        # full sequence match in smaller loop region and at least 95% alignment score for extended sequence
        homolog_set_list = generate_sequence_alignment_for_all_pairs(family_id, loop_list)

        if len(homolog_set_list) == 1:
            logger.info('Filtering ' + family_id)
            new_families[new_family_id].append(homolog_set_list[0].pop())

        else:
            new_families[new_family_id] += families[family_id]

    # print(families)
    # print(new_families)
    # sys.exit()
    return new_families
