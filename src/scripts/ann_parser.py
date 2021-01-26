import glob
import os
import sys
import logging

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from classes import *
from utils import *
# from util_cluster import *
# from residue import *

def parseFR3D(fr3d_fn, detailed=False, model_id="1"):
    ret = []

    if not os.path.isfile(fr3d_fn):
        logger.warning(fr3d_fn[base_path_len:] + ' not found.')
        return ret

    fr3d_fp = open(fr3d_fn)
    pdb_id = os.path.basename(fr3d_fn)[:-5]

    for line in fr3d_fp.readlines():
        pieces = line.strip().split(',')
        for i, piece in enumerate(pieces):
            pieces[i] = piece.strip('"')
        if len(pieces) < 6:
            continue

        if pieces[0].strip().split('|')[1].strip() != model_id or pieces[5].strip().split('|')[1].strip() != model_id:
            continue

        intera_type = pieces[1]
        direction = pieces[2]

        if len(intera_type) == 0 and len(direction) == 0:
            continue

        if intera_type.startswith('c'):
            orientation = 'cis'
            edges = intera_type[1] + '/' + intera_type[2]
        elif intera_type.startswith('t'):
            orientation = 'trans'
            edges = intera_type[1] + '/' + intera_type[2]
        else:
            orientation = ''
        
        if direction == 's35':
            direction = 'upward'
        elif direction == 's53':
            direction = 'downward'
        elif direction == 's55':
            direction = 'outward'
        elif direction == 's33':
            direction = 'inward'
        else:
            direction = ''

        t_op1 = t_op2 = ''

        subpieces = pieces[0].strip().split('|')
        chain_id1 = subpieces[2].strip()
        residue1 = subpieces[3].strip()
        seqnum1 = subpieces[4].strip()
        icode1 = ''
        if len(subpieces) > 7:
            icode1 = subpieces[7].strip()
        if len(subpieces) > 8:
            t_op1 = subpieces[8].strip()

        subpieces = pieces[5].strip().split('|')
        chain_id2 = subpieces[2].strip()
        residue2 = subpieces[3].strip()
        seqnum2 = subpieces[4].strip()
        icode2 = ''
        if len(subpieces) > 7:
            icode2 = subpieces[7].strip()
        if len(subpieces) > 8:
            t_op2 = subpieces[8].strip()

        ind1 = Chainindex(chain_id1, seqnum1, icode1)
        ind2 = Chainindex(chain_id2, seqnum2, icode2)

        index1 = min(ind1, ind2)
        index2 = max(ind1, ind2)

        bp = residue1 + '-' + residue2

        if index1 == index2 and t_op1 != t_op2:
            continue

        if len(orientation) > 0:
            if index1 != ind1:
                edges, bp = get_inverse_bp_info(edges, bp)

            if detailed:
                interact_info = (index1, index2, edges, orientation, bp)
            else:    
                interact_info = (index1, index2, edges, orientation)

            ret.append(interact_info)

        elif len(direction) > 0:
            if index1 != ind1:
                direction, bp = get_inverse_stk_info(direction, bp)

            if detailed:
                interact_info = (index1, index2, direction, bp)
            else:
                interact_info = (index1, index2, direction)

            ret.append(interact_info)

    fr3d_fp.close()
    
    return sorted(list(set(ret)))

def parseDSSR(dssr_fn, detailed=False):
    ret = []
    if not os.path.isfile(dssr_fn):
        logger.warning(dssr_fn[base_path_len:] + ' not found.')
        return ret

    dssr_fp = open(dssr_fn)
    pdb_id = os.path.basename(dssr_fn)[:-8]
    flag_read_bp = False
    flag_read_st = False

    for line in dssr_fp.readlines():
        if line.startswith('List of ') and line.strip().endswith(' base pairs'):
            flag_read_bp = True

        elif line.startswith('List of ') and line.strip().endswith(' non-pairing interactions'):
            flag_read_st = True

        elif flag_read_bp == True and line.strip() == '****************************************************************************':
            flag_read_bp = False

        elif flag_read_st == True and line.strip() == '****************************************************************************':
            flag_read_st = False

        elif flag_read_bp == True and len(line.strip()) > 0:
            pieces = line.strip().split()

            if len(pieces) == 8:     #len(pieces) = 7 for title
                temp = pieces[1].strip().split('.')
                chain_id1 = temp[0].strip() #+ '.' + temp[1][1:]
                seqnum1, icode1 = getDSSRseqnum(temp[1].strip())
                temp = pieces[2].strip().split('.')
                chain_id2 = temp[0].strip() #+ '.' + temp[1][1:]
                seqnum2, icode2 = getDSSRseqnum(temp[1].strip()) #temp[1][1:]
                orientation = 'cis' if pieces[6][0] == 'c' else 'trans'
                edges = pieces[6][1] + '/' + pieces[6][2]

                if not seqnum1.lstrip('-').isdigit() or not seqnum2.lstrip('-').isdigit():  #allow negative numbers
                    fpt = open('dssr_parser_log.txt', 'a')
                    fpt.write('Bad seqnum found. ' + pdb_id + ': '+ seqnum1 + ', ' + seqnum2 + '\n')
                    fpt.close()
                    continue

                if edges[0] == '.' or edges[2] == '.':
                    continue

                if detailed:
                    bp = pieces[3][0]+'-'+pieces[3][2]
                    ret.append((Chainindex(chain_id1, seqnum1, icode1), Chainindex(chain_id2, seqnum2, icode2), edges, orientation, bp.upper()))
                else:
                    ret.append((Chainindex(chain_id1, seqnum1, icode1), Chainindex(chain_id2, seqnum2, icode2), edges, orientation))

        elif flag_read_st == True and len(line.strip()) > 0:
            pieces = line.strip().split()
            
            if len(pieces) > 0 and 'stacking:' in line:
                if 'forward' in line:
                    direction = 'upward'
                elif 'backward' in line:
                    direction = 'downward'
                elif 'inward' in line:
                    direction = 'inward'
                elif 'outward' in line:
                    direction = 'outward'

                temp = pieces[1].split('.')
                chain_id1 = temp[0] #+ '.' + temp[1][1:]
                resname1, seqnum1, icode1 = getDSSRseqnum(temp[1], True) #temp[1][1:]
                temp = pieces[2].split('.')
                chain_id2 = temp[0] #+ '.' + temp[1][1:]
                resname2, seqnum2, icode2 = getDSSRseqnum(temp[1], True) #temp[1][1:]

                if not seqnum1.lstrip('-').isdigit() or not seqnum2.lstrip('-').isdigit():  #allow negative numbers
                    fpt = open('dssr_parser_log.txt', 'a')
                    fpt.write('Bad seqnum found. ' + pdb_id + ': '+ seqnum1 + ', ' + seqnum2 + '\n')
                    fpt.close()
                    continue

                if detailed:
                    bp = resname1+'-'+resname2
                    ret.append((Chainindex(chain_id1, seqnum1, icode1), Chainindex(chain_id2, seqnum2, icode2), direction, bp.upper()))
                else:
                    ret.append((Chainindex(chain_id1, seqnum1, icode1), Chainindex(chain_id2, seqnum2, icode2), direction))

            # elif len(pieces) > 0 and "H-bonds" in line:
            #     direction = "hbond"
            #     ret.append((Chainindex(pieces[1][0], pieces[1][3:], ""), Chainindex(pieces[2][0], pieces[2][3:], ""), direction))
            #break
            #print line.strip()

    dssr_fp.close()

    return ret
    #fp = open(os.path.join(dssr_dir, dssr_fn+'.prs'), 'w')
    # fp.close()

def parseMergedAnnotation(merged_fn, detailed=False):
    ret = []
    if not os.path.isfile(merged_fn):
        logger.warning(merged_fn[base_path_len:] + ' not found.')
        return ret

    merged_fp = open(merged_fn)
    pdb_id = os.path.basename(merged_fn)[:-7]

    for line in merged_fp.readlines():
        if len(line.strip()) > 0:
            pieces = line.strip().split('\t')

            index1 = Chainindex.from_str_index(pieces[0].strip())
            index2 = Chainindex.from_str_index(pieces[1].strip())
            bp = pieces[2].strip()

            if len(pieces) == 5:    #bp
                edges = pieces[3].strip()
                orientation = pieces[4].strip()

                if detailed:
                    ret.append((index1, index2, edges, orientation, bp.upper()))
                else:
                    ret.append((index1, index2, edges, orientation))

            elif len(pieces) == 4:  #stack
                direction = pieces[3].strip()

                if detailed:
                    ret.append((index1, index2, direction, bp.upper()))
                else:
                    ret.append((index1, index2, direction))

            # elif len(pieces) > 0 and "H-bonds" in line:
            #     direction = "hbond"
            #     ret.append((Chainindex(pieces[1][0], pieces[1][3:], ""), Chainindex(pieces[2][0], pieces[2][3:], ""), direction))
            #break
            #print line.strip()

    merged_fp.close()

    return ret

def get_inverse_bp_info(edges, bp):
    edges = edges[2] + edges[1] + edges[0]
    bp1, bp2 = bp.split('-')
    bp = bp2 + '-' + bp1

    return edges, bp

def get_inverse_stk_info(direction, bp):
    rev_direction = {'upward':'downward', 'downward':'upward', 'inward':'inward', 'outward':'outward'}
    direction = rev_direction[direction]
    bp1, bp2 = bp.split('-')
    bp = bp2 + '-' + bp1

    return direction, bp

