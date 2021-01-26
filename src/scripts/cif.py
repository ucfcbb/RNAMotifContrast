import sys
import os
import logging
import random

sys.path.append('../../')
from config import *
sys.path.append(lib_dir)
sys.path.append(scripts_dir)
from utils import *
from classes import *

# def get_nonpolymer_entities_name(pdb_fn):
#     """get the residue level nomenclature mapping for non-polymer entities in the pdb (dict[chain_id]=[,,,])"""
#     ret = {}
#     in_section = False
#     single_entry_flag = False
#     pdb_fp = open(pdb_fn)

#     for line in pdb_fp.readlines():
#         if line.startswith("_pdbx_nonpoly_scheme"):
#             in_section = True

#             ##for single entry
#             pieces = line.strip().split()
#             if len(pieces) > 1:
#                 single_entry_flag = True

#             # if single_entry_flag is True and line.startswith("_pdbx_nonpoly_scheme.auth_seq_num"):
#             #     seqnum = pieces[1]
#             if single_entry_flag is True and line.startswith("_pdbx_nonpoly_scheme.pdb_seq_num"):
#                 seqnum = pieces[1]
#             # if single_entry_flag is True and line.startswith("_pdbx_nonpoly_scheme.auth_mon_id"):
#             #     resname = pieces[1]
#             if single_entry_flag is True and line.startswith("_pdbx_nonpoly_scheme.pdb_mon_id"):
#                 resname = pieces[1]
#             if single_entry_flag is True and line.startswith("_pdbx_nonpoly_scheme.pdb_strand_id"):
#                 chain_id = pieces[1]
#             if single_entry_flag is True and line.startswith("_pdbx_nonpoly_scheme.pdb_ins_code"):
#                 icode = pieces[1]
#                 if icode == ".":
#                     icode = ""
#                 if chain_id not in ret:
#                     ret[chain_id] = []
#                 index = Chainindex(chain_id, seqnum, icode)
#                 residue = Residue(index, resname)
#                 if residue not in ret[chain_id]:
#                     ret[chain_id].append(residue)
#                 # if resname not in ret[chain_id]:
#                 #     ret[chain_id].append(resname)
#             # ignoring the icode, because it's not required here. We're just checking the name to continue while selecting the residues.
#             ##for single entry

#             continue
#         elif line.startswith("#") and in_section is True:
#             in_section = False
#             break
#         elif in_section is True and line.startswith(";"):
#             continue
#         elif in_section is True:
#             decom = line.strip().split()
#             resname = decom[6]  #taking non auth part (pdb_mon_id)
#             # resname = decom[7]  #taking auth part (auth_mon_id)
#             chain_id = decom[8]
#             seqnum = decom[4]   # (pdb_seq_num)
#             # seqnum = decom[5] # (auth_seq_num)
#             icode = decom[9]
#             if icode == ".":
#                 icode = ""
#             # ignoring the icode, because it's not required here. We're just checking the name to continue while selecting the residues.
#             if chain_id not in ret:
#                 ret[chain_id] = []

#             index = Chainindex(chain_id, seqnum, icode)
#             residue = Residue(index, resname)
#             if residue not in ret[chain_id]:
#                 ret[chain_id].append(residue)
#             # if resname not in ret[chain_id]:
#             #     ret[chain_id].append(resname)
#     pdb_fp.close()
#     return ret

# def get_residue_with_model(pdb_fn, nonpolymer_dict, model_id):
#     """get the residues (dict[chain_id]=[,,,])"""
#     ret = {}

#     current_index = Chainindex('\0', -65535, '\0')
#     #chain_ter = set() #no use here
#     in_section = False
#     pdb_fp = open(pdb_fn)
#     for line in pdb_fp.readlines():
#         if line.startswith("ATOM") or line.startswith("HETATM"):
#             # print(line)
#             in_section = True
#             decom = line.strip().split()
#             model_no = decom[20]
#             if model_no != model_id:    #do not take from other model number
#                 continue
#             chain_id = decom[18] #line[79]
#             if ret.get(chain_id) is None:
#                 ret[chain_id] = []
#             resname = decom[17] #line[20:23].strip()
#             #seqnum = line[28:32].strip()
#             seqnum = decom[16] #line[71:75].strip()
#             icode = decom[9] #line[32].strip()
#             if icode == "?":
#                 icode = ""
#             index = Chainindex(chain_id, seqnum, icode)

#             if chain_id in nonpolymer_dict.keys() and Residue(index, resname) in nonpolymer_dict[chain_id]:
#                 #print resname
#                 continue

#             if index != current_index:
#                 ret[chain_id].append(Residue(index, resname))
#                 current_index = copy.deepcopy(index)
#         elif in_section is True and line.startswith("#"):
#             break
#     pdb_fp.close()
#     return ret

# def get_residue(pdb_fn, nonpolymer_dict):
#     return get_residue_with_model(pdb_fn, nonpolymer_dict, "1")

# def get_missing_residue(pdb_fn):
#     if not os.path.isfile(pdb_fn):
#         return {}
#     """get the missing residues in the pdb (dict[chain_id]=[,,,])"""
#     ret = {}
#     in_section = False
#     pdb_fp = open(pdb_fn)

#     for line in pdb_fp.readlines():
#         if line.startswith("_pdbx_unobs_or_zero_occ_residues"):
#             in_section = True
#             continue
#         elif line.startswith("#") and in_section is True:
#             in_section = False
#             break
#         elif in_section is True and line.startswith(";"):
#             continue
#         elif in_section is True:
#             decom = line.strip().split()
#             model = decom[1] #line[3]
#             if model != " " and model != "1":
#                 continue
#             resname = decom[5] #line[11:14].strip()
#             chain_id = decom[4] #line[9]
#             seqnum = decom[6] #line[14:18].strip()
#             icode = decom[7] #line[19].strip()
#             if icode == "?":
#                 icode = ""
#             if chain_id not in ret:
#                 ret[chain_id] = []

#             index = Chainindex(chain_id, seqnum, icode)
#             ret[chain_id].append(Residue(index, resname))
#     pdb_fp.close()
#     return ret

def insert_missing_residue(residue_dict, missing_residue_dict):
    """insert the misssing residues back"""
    ret = {}
    # print "len: " + str((residue_dict["T"])) + ", " + str((missing_residue_dict["T"]))
    for chain_id in residue_dict:
        # print "\n\nchain id: " + chain_id
        if chain_id not in missing_residue_dict:
            ret[chain_id] = residue_dict[chain_id]
        else:
            i = j = 0
            residue_list = residue_dict[chain_id]
            missing_residue_list = missing_residue_dict[chain_id]
            # print "now: " + str((residue_list)) + ", " + str((missing_residue_list))
            # sys.exit()
            ret[chain_id] = []
            while i < len(residue_list) and j < len(missing_residue_list):
                # print "comparing " + str(residue_list[i]) + " and " + str(missing_residue_list[j])
                if residue_list[i].index < missing_residue_list[j].index:
                    ret[chain_id].append(residue_list[i])
                    i += 1
                elif residue_list[i].index > missing_residue_list[j].index:     #conditional problem here
                    ret[chain_id].append(missing_residue_list[j])
                    j += 1
                else:
                    ret[chain_id].append(residue_list[i])
                    i += 1
                    j += 1
            if i == len(residue_list):
                ret[chain_id] += missing_residue_list[j:]
            if j == len(missing_residue_list):
                ret[chain_id] += residue_list[i:]
    return ret

# def get_modified_residue(pdb_fn):
#     """get modified residues in the pdb (dict[Chainindex]=[(resname1, resnam2),])"""
#     ret = {}
#     in_section = False
#     pdb_fp = open(pdb_fn)

#     for line in pdb_fp.readlines():
#         if line.startswith("_pdbx_struct_mod_residue"):
#             in_section = True
#             continue
#         elif line.startswith("#") and in_section is True:
#             in_section = False
#             break
#         elif in_section is True and line.startswith(";"):
#             continue
#         elif in_section is True:
#             decom = line.strip().split()
#             resname1 = decom[5] #taking auth part
#             chain_id = decom[4]
#             seqnum = decom[6]
#             icode = decom[7]
#             if icode == "?":
#                 icode = ""
#             resname2 = decom[8]
#             if resname2 == "" or resname2 == "?":
#                 resname2 = resname1
#             index = Chainindex(chain_id, seqnum, icode)
#             ret[index] = (resname1, resname2)
#     pdb_fp.close()
#     return ret

def replace_modified_residue(pdb_id, chain_id, residue_list, modified_residue_dict, residue_abbreviation_dict, amino_acid_list):
    """replace the modified residues to one letter residues"""
    ret = []
    # print "residue list: " + str(residue_list)
    for residue in residue_list:
        # print "residue_index: " + str(residue.index) + ", residue_symbol: " + str(residue.symbol)
        # print residue.index in modified_residue_dict
        # if residue.index in modified_residue_dict:
        #     print residue.symbol == modified_residue_dict[residue.index][0]
        if residue.index in modified_residue_dict and residue.symbol == modified_residue_dict[residue.index][0]:
            if modified_residue_dict[residue.index][1] in residue_abbreviation_dict:
                residue.symbol = residue_abbreviation_dict[modified_residue_dict[residue.index][1]]
            else:
                logger.warning(pdb_id + '\tchain id: ' + chain_id + ', ' + str(residue) + ': PARENT "' + modified_residue_dict[residue.index][1] + '" NOT FOUND.')
                continue
        elif residue.symbol in residue_abbreviation_dict:
            # print "previous SYMBOL: " + residue.symbol + ", "
            residue.symbol = residue_abbreviation_dict[residue.symbol]
            # print "current SYMBOL: " + residue.symbol
        else:
            #missed it somehow
            if residue.symbol in amino_acid_list:
                logger.warning(pdb_id + '\tchain id: ' + chain_id + ', ' + str(residue) + ': NOT FOUND; FOUND in Amino List.')
            else:
                # logger.warning(pdb_id + '\tchain id: ' + chain_id + ', ' + str(residue) + ': NOT FOUND.')
                residue.symbol = 'X'
                ret.append(residue)
            continue
            #not added to the list
        ret.append(residue)
    return ret

# def replace_unknown_letter_in_ref(ref_seq, residue_seq, pdb_id):
def replace_unknown_letter_in_ref(ref_seq, residue_seq):
    """
    replace unknown letters in the reference sequence, such as "X" or "N"
    :return: the updated reference sequence which only contains "ACGU"
    """

    if len(ref_seq) != len(residue_seq):
        return None
    ref_seq_list = list(ref_seq)
    for i in range(len(ref_seq_list)):
        if ref_seq_list[i] not in ["A", "C", "G", "U"] and residue_seq[i] in ["A", "C", "G", "U"]:
            ref_seq_list[i] = residue_seq[i]
        elif ref_seq_list[i] not in ["A", "C", "G", "U"] and residue_seq[i] not in ["A", "C", "G", "U"]:
            
            # if ref_seq_list[i] not in unknown_letter_dict:
            #     unknown_letter_dict[ref_seq_list[i]] = []
            #     unknown_letter_count[ref_seq_list[i]] = 0
            # if pdb_id not in unknown_letter_dict[ref_seq_list[i]]:
            #     unknown_letter_dict[ref_seq_list[i]].append(pdb_id)
            # unknown_letter_count[ref_seq_list[i]] += 1

            # ref_seq_list[i] = random.choice(["A", "C", "G", "U"])
            ref_seq_list[i] = "X"

    ref_seq_replaced = "".join(ref_seq_list)
    return ref_seq_replaced






def load_pdb_data_from_file(pdb_fname):
    nonpolymer_entities = {}
    residues = {}
    missing_residues = {}
    modified_residues = {}

    fp = open(pdb_fname)
    lines = fp.readlines()
    fp.close()

    # in_a_section = False
    has_multiple_rows = False

    total_lines = len(lines)
    line_no = 0
    while line_no < total_lines:
        line = lines[line_no]

        if line.startswith('loop_'):
            has_multiple_rows = True

        elif line.startswith('#'):
            # in_a_section = False
            has_multiple_rows = False

        elif line.startswith('_pdbx_nonpoly_scheme.'):
            line_no, nonpolymer_entities = load_nonpolymers(lines, line_no, total_lines, has_multiple_rows)
            # in_a_section = False
            has_multiple_rows = False

        elif line.startswith('_atom_site.'):
            line_no, residues = load_residues(lines, line_no, total_lines, has_multiple_rows)
            has_multiple_rows = False

        elif line.startswith('_pdbx_unobs_or_zero_occ_residues.'):
            line_no, missing_residues = load_missing_residues(lines, line_no, total_lines, has_multiple_rows)
            has_multiple_rows = False

        elif line.startswith('_pdbx_struct_mod_residue.'):
            line_no, modified_residues = load_modified_residues(lines, line_no, total_lines, has_multiple_rows)
            has_multiple_rows = False

        line_no += 1

    filter_nonpoly_dict(residues, nonpolymer_entities)

    return residues, missing_residues, modified_residues

def filter_nonpoly_dict(residues, nonpolymer_entities):
    for chain_id in nonpolymer_entities:
        if chain_id in residues:
            for res in nonpolymer_entities[chain_id]:
                if res in residues[chain_id]:
                    residues[chain_id].remove(res)
            if len(residues[chain_id]) == 0:
                residues.pop(chain_id)

def loop_through_section_end(lines, line_no, total_lines):
    while line_no < total_lines and not lines[line_no].startswith('#'):
        line_no += 1

    return line_no

def load_nonpolymers(lines, line_no, total_lines, has_multiple_rows):
    non_poly_data = {}
    
    seqnum = ''
    resname = ''
    chain_id = ''
    icode = ''

    #### TODO: take auth item if exists
    index_dict = {}
    # index_dict['pdb_seq_num'] = 4
    index_dict['auth_seq_num'] = 5
    # index_dict['pdb_mon_id'] = 6
    index_dict['auth_mon_id'] = 7
    index_dict['pdb_strand_id'] = 8
    index_dict['pdb_ins_code'] = 9

    required_data_count = 0
    item_counter = 0
    while line_no < total_lines:
        line = lines[line_no]

        if line.startswith('#'):
            break
        if line.startswith(';'):
            line_no += 1
            continue

        pieces = line.strip().split()
        if has_multiple_rows == False:
            
            # if line.startswith('_pdbx_nonpoly_scheme.pdb_seq_num'):
            if line.startswith('_pdbx_nonpoly_scheme.auth_seq_num'):
                seqnum = pieces[1]
                required_data_count += 1

            # elif line.startswith('_pdbx_nonpoly_scheme.pdb_mon_id'):
            elif line.startswith('_pdbx_nonpoly_scheme.auth_mon_id'):
                resname = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_nonpoly_scheme.pdb_strand_id'):
                chain_id = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_nonpoly_scheme.pdb_ins_code'):
                icode = pieces[1]
                if icode == '.':
                    icode = ''
                required_data_count += 1

            if required_data_count == 4:
                if chain_id not in non_poly_data:
                    non_poly_data[chain_id] = []

                index = Chainindex(chain_id, seqnum, icode)
                residue = Residue(index, resname)
                if residue not in non_poly_data[chain_id]:
                    non_poly_data[chain_id].append(residue)

                line_no = loop_through_section_end(lines, line_no, total_lines)
                break
        else:
            if line.startswith('_'):
                item = line.strip().split('.')[1]
                if item in index_dict:
                    index_dict[item] = item_counter
                item_counter += 1

            else:
                pieces = line.strip().split()
                # resname = pieces[index_dict['pdb_mon_id']]  #taking non auth part (pdb_mon_id)
                resname = pieces[index_dict['auth_mon_id']]  #taking auth part (auth_mon_id)

                chain_id = pieces[index_dict['pdb_strand_id']]
                
                # seqnum = pieces[index_dict['pdb_seq_num']]   # (pdb_seq_num)
                seqnum = pieces[index_dict['auth_seq_num']] # (auth_seq_num)

                icode = pieces[index_dict['pdb_ins_code']]
                if icode == ".":
                    icode = ""
                # ignoring the icode, because it's not required here. We're just checking the name to continue while selecting the residues.
                if chain_id not in non_poly_data:
                    non_poly_data[chain_id] = []

                index = Chainindex(chain_id, seqnum, icode)
                residue = Residue(index, resname)

                if residue not in non_poly_data[chain_id]:
                    non_poly_data[chain_id].append(residue)

        line_no += 1

    return line_no, non_poly_data

def load_residues(lines, line_no, total_lines, has_multiple_rows):
    residue_data = {}

    model_id = '1'
    model_no = ''
    chain_id = ''
    resname = ''
    seqnum = ''
    icode = ''
    current_index = Chainindex('\0', -65535, '\0')

    index_dict = {}
    index_dict['pdbx_PDB_ins_code'] = 9
    index_dict['auth_seq_id'] = 16
    index_dict['auth_comp_id'] = 17
    index_dict['auth_asym_id'] = 18
    index_dict['pdbx_PDB_model_num'] = 20

    item_counter = 0
    while line_no < total_lines:
        line = lines[line_no]

        if line.startswith('#'):
            break
        if line.startswith(';'):
            line_no += 1
            continue

        pieces = line.strip().split()
        if has_multiple_rows == True:
            if line.startswith('_'):
                item = line.strip().split('.')[1]
                if item in index_dict:
                    index_dict[item] = item_counter
                item_counter += 1

            if line.startswith('ATOM') or line.startswith('HETATM'):
                pieces = line.strip().split()
                model_no = pieces[index_dict['pdbx_PDB_model_num']]
                
                if model_no == model_id:    #do not take from other model number
                    
                    chain_id = pieces[index_dict['auth_asym_id']] #line[79]
                    if chain_id not in residue_data:
                        residue_data[chain_id] = []

                    resname = pieces[index_dict['auth_comp_id']] #line[20:23].strip()   # inconsistent with nonpoly data (auth/?)
                    #seqnum = line[28:32].strip()
                    seqnum = pieces[index_dict['auth_seq_id']] #line[71:75].strip()    # inconsistent with nonpoly data (auth/?)
                    icode = pieces[index_dict['pdbx_PDB_ins_code']] #line[32].strip()
                    
                    if icode == '?':
                        icode = ''
                    index = Chainindex(chain_id, seqnum, icode)

                    residue = Residue(index, resname)
                    # if not (chain_id in nonpolymer_data and residue in nonpolymer_data[chain_id]):
                        #print resname
                        # line_no += 1
                        # continue

                    if index != current_index:
                        residue_data[chain_id].append(residue)
                        current_index = copy.deepcopy(index)
            
        else:
            pass
            # it supposed to have multiple rows always

        line_no += 1

    return line_no, residue_data

def load_missing_residues(lines, line_no, total_lines, has_multiple_rows):
    missing_residue_data = {}

    model_id = '1'
    model_no = ''
    chain_id = ''
    resname = ''
    seqnum = ''
    icode = ''

    index_dict = {}
    index_dict['PDB_model_num'] = 1
    index_dict['auth_asym_id'] = 4
    index_dict['auth_comp_id'] = 5
    index_dict['auth_seq_id'] = 6
    index_dict['PDB_ins_code'] = 7

    required_data_count = 0
    item_counter = 0
    while line_no < total_lines:
        line = lines[line_no]

        if line.startswith('#'):
            break
        if line.startswith(';'):
            line_no += 1
            continue

        pieces = line.strip().split()
        if has_multiple_rows == False:
            
            if line.startswith('_pdbx_unobs_or_zero_occ_residues.PDB_model_num'):
                model_no = pieces[1]
                if model_id != model_no:
                    line_no = loop_through_section_end(lines, line_no, total_lines)
                    break
                    # return line_no, {}

            elif line.startswith('_pdbx_unobs_or_zero_occ_residues.auth_asym_id'):
                chain_id = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_unobs_or_zero_occ_residues.auth_comp_id'):
                resname = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_unobs_or_zero_occ_residues.auth_seq_id'):
                seqnum = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_unobs_or_zero_occ_residues.PDB_ins_code'): 
                icode = pieces[1]
                if icode == '?':
                    icode = ''
                required_data_count += 1

            if required_data_count == 4:
                if chain_id not in missing_residue_data:
                    missing_residue_data[chain_id] = []

                index = Chainindex(chain_id, seqnum, icode)
                residue = Residue(index, resname)

                if residue not in missing_residue_data[chain_id]:
                    missing_residue_data[chain_id].append(residue)

                line_no = loop_through_section_end(lines, line_no, total_lines)
                break
            
        else:
            if line.startswith('_'):
                item = line.strip().split('.')[1]
                if item in index_dict:
                    index_dict[item] = item_counter
                item_counter += 1

            else:
                model_no = pieces[index_dict['PDB_model_num']]
                
                if model_id == model_no:
                    chain_id = pieces[index_dict['auth_asym_id']] #line[9]
                    resname = pieces[index_dict['auth_comp_id']] #line[11:14].strip()
                    seqnum = pieces[index_dict['auth_seq_id']] #line[14:18].strip()
                    icode = pieces[index_dict['PDB_ins_code']] #line[19].strip()

                    if icode == '?':
                        icode = ''
                    if chain_id not in missing_residue_data:
                        missing_residue_data[chain_id] = []

                    index = Chainindex(chain_id, seqnum, icode)
                    residue = Residue(index, resname)

                    if residue not in missing_residue_data[chain_id]:
                        missing_residue_data[chain_id].append(residue)

        line_no += 1

    return line_no, missing_residue_data

def load_modified_residues(lines, line_no, total_lines, has_multiple_rows):
    modified_residue_data = {}

    chain_id = ''
    resname1 = ''
    resname2 = ''
    seqnum = ''
    icode = ''

    index_dict = {}
    index_dict['auth_asym_id'] = 4
    index_dict['auth_comp_id'] = 5
    index_dict['auth_seq_id'] = 6
    index_dict['PDB_ins_code'] = 7
    index_dict['parent_comp_id'] = 8

    required_data_count = 0
    item_counter = 0
    while line_no < total_lines:
        line = lines[line_no]

        if line.startswith('#'):
            break
        if line.startswith(';'):
            line_no += 1
            continue

        pieces = line.strip().split()
        if has_multiple_rows == False:

            if line.startswith('_pdbx_struct_mod_residue.auth_asym_id'):
                chain_id = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_struct_mod_residue.auth_comp_id'):
                resname1 = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_struct_mod_residue.auth_seq_id'):
                seqnum = pieces[1]
                required_data_count += 1

            elif line.startswith('_pdbx_struct_mod_residue.PDB_ins_code'): 
                icode = pieces[1]
                if icode == '?':
                    icode = ''
                required_data_count += 1

            elif line.startswith('_pdbx_struct_mod_residue.parent_comp_id'):
                resname2 = pieces[1]
                if resname2 == '' or resname2 == '?':
                    resname2 = resname1
                required_data_count += 1

            if required_data_count == 5:
                index = Chainindex(chain_id, seqnum, icode)
                modified_residue_data[index] = (resname1, resname2)

                line_no = loop_through_section_end(lines, line_no, total_lines)
                break
            
        else:
            if line.startswith('_'):
                item = line.strip().split('.')[1]
                if item in index_dict:
                    index_dict[item] = item_counter
                item_counter += 1

            else:
                chain_id = pieces[index_dict['auth_asym_id']]
                resname1 = pieces[index_dict['auth_comp_id']] #taking auth part
                seqnum = pieces[index_dict['auth_seq_id']]
                icode = pieces[index_dict['PDB_ins_code']]
                if icode == '?':
                    icode = ''
                resname2 = pieces[index_dict['parent_comp_id']]

                if resname2 == '' or resname2 == '?':
                    resname2 = resname1

                index = Chainindex(chain_id, seqnum, icode)
                modified_residue_data[index] = (resname1, resname2)

        line_no += 1

    return line_no, modified_residue_data
