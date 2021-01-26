import sys
import os
import glob
import time
import logging
import operator
import multiprocessing as mp
from functools import reduce

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from utils import *
from classes import *

def is_pdb_ind_within_ext_limit(a_pdb, b_pdb, i_pdb, extend):
    a_ind = int(a_pdb.seqnum)
    b_ind = int(b_pdb.seqnum)
    i_ind = int(i_pdb.seqnum)

    if i_pdb.icode != "":       # For the selected dataset, only 1HMH.cif has such case that seems somehow closed to the other residues.
        return True

    if i_ind >= a_ind - extend and i_ind <= b_ind + extend:
        return True
    return False

def add_adjacent_residue_info(pdb_id, loop, pdb_chain_adj_dict):
    fname = os.path.join(neares_protein_data_dir, loop + '.adj_info')
    if os.path.isfile(fname):
        fp = open(fname)
        lines = fp.readlines()
        fp.close()
        load_items = lines[2].strip().split(',')
        # print(load_items)
        if pdb_id not in pdb_chain_adj_dict:
            pdb_chain_adj_dict[pdb_id] = {}
            
        for item in load_items:
            chain_id = item.strip()
            regions = None
            if ':' in item:
                chain_id, regions = item.strip().split(':')

            if chain_id not in pdb_chain_adj_dict[pdb_id]:
                pdb_chain_adj_dict[pdb_id][chain_id] = {}
            # if loop not in pdb_chain_adj_dict[pdb_id][chain_id]:
            #     pdb_chain_adj_dict[pdb_id][chain_id][loop] = []

            if regions != None:
                regions = regions.strip().split('_')
                regions = map(lambda x: tuple(x.strip().split('-')), regions)
                regions = map(lambda x: (int(x[0].strip().split('.')[0]), int(x[1].strip().split('.')[0])), regions)

            pdb_chain_adj_dict[pdb_id][chain_id][loop] = regions

def generate_pdb_chain_dict(loop_list, extend):
    converterDict = {}
    pdb_chain_dict = {}
    pdb_chain_adj_dict = {}

    for loop in loop_list:
        pdb_chain, segments = loop.strip().split(':')
        pdb_id, chain_id = pdb_chain.strip().split('_')
        
        if pdb_id not in pdb_chain_dict:
            pdb_chain_dict[pdb_id] = {}
        if chain_id not in pdb_chain_dict[pdb_id]:
            pdb_chain_dict[pdb_id][chain_id] = {}
        if loop not in pdb_chain_dict[pdb_id][chain_id]:
            pdb_chain_dict[pdb_id][chain_id][loop] = []

        if pdb_chain not in converterDict.keys():
            mapping_file_name = pdb_chain + '.rmsx.nch'
            converter = PDB_FASTA_Index_Converter(pdb_fasta_mapping_dir, mapping_file_name)
            converterDict[pdb_chain] = converter

        segments = segments.strip().split('_')

        for segment in segments:
            a, b = segment.strip().split('-')
            a = int(a)
            b = int(b)
            if a > b:
                # print mapping_file_name
                # print a
                # print b
                logger.error('ERROR: sequence flipped for ' + str(loop))
                sys.exit()
            
            a_pdb = converterDict[pdb_chain].convert_FASTAindx_To_PDBindx(str(a))
            b_pdb = converterDict[pdb_chain].convert_FASTAindx_To_PDBindx(str(b))

            for i_fasta in range(a-extend, b+1+extend):
                i_pdb = converterDict[pdb_chain].convert_FASTAindx_To_PDBindx(str(i_fasta))
                if i_pdb != None and is_pdb_ind_within_ext_limit(a_pdb, b_pdb, i_pdb, extend):
                    ii = (int(i_pdb.seqnum), i_pdb.icode)
                    if ii not in pdb_chain_dict[pdb_id][chain_id][loop]:
                        pdb_chain_dict[pdb_id][chain_id][loop].append(ii)

        add_adjacent_residue_info(pdb_id, loop, pdb_chain_adj_dict)

    return pdb_chain_dict, pdb_chain_adj_dict

def generate_partial_pdb_files_for_single_pdb(partial_pdbx_dir, pdb_id, chain_id_seqnum_dict, adj_res_dict, isPymol=True, model_id='1'):

    fp = open(os.path.join(pdbx_dir, pdb_id + '.cif'))
    fw = {}
    for chain_id in chain_id_seqnum_dict:
        fw[chain_id] = {}
        for loop_name in chain_id_seqnum_dict[chain_id]:
            fw[chain_id][loop_name] = open(os.path.join(partial_pdbx_dir, loop_name + '.cif'), 'w')
    # fw = open(os.path.join(partial_pdbx_dir, pdb_id + '.cif'), 'w')

    lines = fp.readlines()
    fp.close()

    for chain_id in fw:
        for loop_name in fw[chain_id]:
            for i in range(4):
                fw[chain_id][loop_name].write(lines[i])

    # for i in range(4):
    #     fw.write(lines[i])
    # fw.write('loop_\n')

    total_lines = len(lines)
    in_poly_seq_section = False
    in_atom_site_section = False
    for i in range(4, total_lines):
        line = lines[i]

        if line.startswith('_atom_site.'):
            if in_atom_site_section == False:
                # fw.write('loop_\n')
                for chain_id in fw:
                    for loop_name in fw[chain_id]:
                        fw[chain_id][loop_name].write('loop_\n')
                in_atom_site_section = True
            # fw.write(line)
            for chain_id in fw:
                for loop_name in fw[chain_id]:
                    fw[chain_id][loop_name].write(line)

        elif in_atom_site_section == True:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                decom = line.strip().split()
                
                model_no = decom[20]
                chain_id = decom[18] #line[79]
                seqnum = int(decom[16]) #line[71:75].strip()
                icode = decom[9] #line[32].strip()
                if icode == "?":
                    icode = ""

                if model_no != model_id:    #do not take from other model number
                    continue
                if chain_id not in chain_id_seqnum_dict and chain_id not in adj_res_dict:
                    continue

                if chain_id in chain_id_seqnum_dict:
                    for loop_name in chain_id_seqnum_dict[chain_id]:
                        if (seqnum, icode) in chain_id_seqnum_dict[chain_id][loop_name]:
                            fw[chain_id][loop_name].write(line)

                if chain_id in adj_res_dict:
                    for loop_name in adj_res_dict[chain_id]:
                        loop_chain_id = loop_name.strip().split(':')[0].strip().split('_')[1]
                        
                        if adj_res_dict[chain_id][loop_name] == None:
                            fw[loop_chain_id][loop_name].write(line)
                        
                        else:
                            for s, e in adj_res_dict[chain_id][loop_name]:
                                if s <= seqnum and seqnum <= e:
                                    fw[loop_chain_id][loop_name].write(line)
                                    break

                # if seqnum not in chain_id_seqnum_list[chain_id]:
                #     continue

            elif line.startswith('#'):
                in_atom_site_section = False
                for chain_id in fw:
                    for loop_name in fw[chain_id]:
                        fw[chain_id][loop_name].write(line)

            # fw.write(line)

        elif not isPymol and line.startswith('_entity_poly_seq.'):
            if in_poly_seq_section == False:
                # fw.write('loop_\n')
                for chain_id in fw:
                    for loop_name in fw[chain_id]:
                        fw[chain_id][loop_name].write('loop_\n')
                in_poly_seq_section = True
            # fw.write(line)
            for chain_id in fw:
                for loop_name in fw[chain_id]:
                    fw[chain_id][loop_name].write(line)
        elif not isPymol and in_poly_seq_section == True:
            if line.startswith('#'):
                in_poly_seq_section = False
            # fw.write(line)
            for chain_id in fw:
                for loop_name in fw[chain_id]:
                    fw[chain_id][loop_name].write(line)

    
    # fw.close()
    for chain_id in chain_id_seqnum_dict:
        for loop_name in chain_id_seqnum_dict[chain_id]:
            fw[chain_id][loop_name].close()

def _generate_partial_pdb_files_worker(p):
    generate_partial_pdb_files_for_single_pdb(*p)

def prepare_partial_pdbs(partial_pdbx_dir, loop_node_list_str, loop_cif_extension):
    create_directory(partial_pdbx_dir)

    loop_list = []
    for loop_name in loop_node_list_str:
        partial_pdb_fname = os.path.join(partial_pdbx_dir, loop_name + '.cif')
        if not os.path.isfile(partial_pdb_fname):
            loop_list.append(loop_name)

    if len(loop_list) == 0:
        logger.info('Using existing partial pdb files in ' + partial_pdbx_dir[base_path_len:] + '.')
        print('')
        return

    logger.info('Generating partial pdb files (loop.cif).')
    start_time = time.time()    
    # remove_all_from_dir(partial_pdbx_dir)

    pdb_chain_dict, pdb_chain_adj_dict = generate_pdb_chain_dict(loop_list, loop_cif_extension)

    parameter_list = []
    for pdb_id in pdb_chain_dict:
        if pdb_id in pdb_chain_adj_dict:
            parameter_list.append((partial_pdbx_dir, pdb_id, pdb_chain_dict[pdb_id], pdb_chain_adj_dict[pdb_id]))
        else:
            parameter_list.append((partial_pdbx_dir, pdb_id, pdb_chain_dict[pdb_id], []))
    pool = mp.Pool(number_of_multiprocess)
    pool.map(_generate_partial_pdb_files_worker, parameter_list)

    wait_for_certain_time_according_to_wait_factor(len(parameter_list))
    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')
