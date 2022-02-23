import sys
import os
import copy
import re

class Chainindex():
    """the index of residue in the PDB chain"""
    def __init__(self, chain_id, seqnum, icode):
        """set icode to empty if no icode"""
        self.chain_id = chain_id
        self.seqnum = int(seqnum)   # seqnum is treated as integer
        self.icode = icode

    def __repr__(self):
        ret = ''
        if self.chain_id.isalpha():
            ret += self.chain_id+str(self.seqnum)
        else:
            ret += "'"+self.chain_id+"'"+str(self.seqnum)

        return ret if self.icode == "" else ret+"."+self.icode

    def __lt__(self, other):
        if isinstance(other, self.__class__) and self.chain_id == other.chain_id:
            return (self.seqnum < other.seqnum) or (self.seqnum == other.seqnum and self.icode < other.icode)
        else:
            if isinstance(other, self.__class__) and self.chain_id < other.chain_id:
                return True
            return False

    def __gt__(self, other):
        if isinstance(other, self.__class__) and self.chain_id == other.chain_id:
            return (self.seqnum > other.seqnum) or (self.seqnum == other.seqnum and self.icode > other.icode)
        else:
            if isinstance(other, self.__class__) and self.chain_id > other.chain_id:
                return True
            return False

    def __eq__(self, other):
        if isinstance(other, self.__class__) and self.chain_id == other.chain_id:
            return (self.seqnum == other.seqnum) and (self.icode == other.icode)
        else:
            return False

    def __ne__(self, other):
            return not self.__eq__(other)

    def __hash__(self):
        return hash((self.chain_id, self.seqnum, self.icode))

    @classmethod
    def from_mca_index(cls, mca_str):
        """construct Seqindex class from a string (MCA format)"""
        if mca_str.startswith("'"):
            chain_id = mca_str[1]
            i = 3
        else:
            chain_id = mca_str[0]
            i = 1

        if mca_str[-1].isalpha():
            icode = mca_str[-1]
            j = len(mca_str)-2
        else:
            icode = ""
            j = len(mca_str)

        seqnum = mca_str[i:j]
        return cls(chain_id, seqnum, icode)

    @classmethod
    def from_dssr_index(cls, mca_str):
        """construct Seqindex class from a string (DSSR format)"""
        if mca_str.startswith("'"):
            chain_id = mca_str[1]
            i = 3
        else:
            chain_id = mca_str[0]
            i = 1

        if mca_str[-1].isalpha():
            icode = mca_str[-1]
            j = len(mca_str)-2
        else:
            icode = ""
            j = len(mca_str)

        seqnum = mca_str[i:j]
        return cls(chain_id, seqnum, icode)

    @classmethod
    def from_str_index(cls, str):
        """construct Seqindex class from a string (DSSR format)"""
        if str.startswith("'"):
            chain_id = str[1:].strip().split("'")[0]
            i = len(chain_id)+2
        else:
            chain_id = re.split('-?(\d+)',str)[0]
            i = len(chain_id)

        if str[-1].isalpha():
            icode = str[-1]
            j = len(str)-2
        else:
            icode = ""
            j = len(str)

        seqnum = str[i:j]
        return cls(chain_id, seqnum, icode)


class Residue():
    def __init__(self, index, symbol):
        self.index = copy.deepcopy(index)
        self.symbol = symbol

    def __repr__(self):
        return str(self.index)+"("+self.symbol+")"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.index == other.index and self.symbol == other.symbol
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.index) ^ hash(self.symbol)

class StackTNode:
    def __init__(self, i, j, size, children, has_cano = False):
        self.i = i
        self.j = j
        self.size = size
        self.children = children
        self.has_cano = has_cano

    def __repr__(self):
        return "(%d, %d)<%d>" % (self.i, self.j, self.size)

    def set_cano_status(self, status=True):
        self.has_cano = status

    def display(self, lvl=0):
        # print '.'*3*lvl + str(self)
        if self.children:
            for node in self.children:
                node.display(lvl=lvl+1)

    def fail_cWW_stack_len_threshold(self, cutoff):
        new_cutoff = cutoff
        if self.has_cano == False:
            # setting max value to remove continuous cWW without canonical as cWW_stack
            new_cutoff = 100000
        return self.size < new_cutoff

    def filter_cWW_stack(self, stack_size_cutoff):
        # check if there is any children that fails the threshold
        has_failing_children = True
        while has_failing_children: # len(filter(lambda x: x.fail_cWW_stack_len_threshold(stack_size_cutoff), self.children)) != 0:
            new_children = []
            has_failing_children = False
            for node in self.children:
                #if node.size<stack_size_cutoff:
                if node.fail_cWW_stack_len_threshold(stack_size_cutoff):
                    new_children += node.children
                    has_failing_children = True
                else:
                    new_children.append(node)
            self.children = new_children

        for node in self.children:
            node.filter_cWW_stack(stack_size_cutoff)

    def get_loop_regions(self):
        start = self.i+self.size-1
        loop_regions = []
        for r in self.children:
            loop_regions.append((start, r.i))
            start = r.j
        loop_regions.append((start, self.j-self.size+1))
        return loop_regions

    def loop(self, ret):
        loop_regions = self.get_loop_regions()

        if self.i != -1:
            suffix = ""
            # if self.has_cano == False:
            #     suffix = "_nc-cWW_stack"
            if len(self.children) == 0:
                ret.append((loop_regions, "HL"+suffix))
            elif len(self.children) == 1:
                ret.append((loop_regions, "IL"+suffix))
            elif len(self.children) > 1:
                ret.append((loop_regions, "ML"+suffix))

        for node in self.children:
            node.loop(ret)

    # def loop(self, ret):
    #     start = self.i+self.size
    #     loop_regions = []
    #     for r in self.children:
    #         if start <= r.i-1:
    #             loop_regions.append((start, r.i-1))
    #         start = r.j+1
    #     if start <= self.j-self.size:
    #         loop_regions.append((start, self.j-self.size))
    #
    #     if self.i != -1:
    #         if len(self.children) == 0:
    #             ret.append((loop_regions, "HL"))
    #         elif len(self.children) == 1:
    #             ret.append((loop_regions, "IL"))
    #         elif len(self.children) > 1:
    #             ret.append((loop_regions, "ML"))
    #
    #     for node in self.children:
    #         node.loop(ret)

class PDB_FASTA_Index_Converter:
    def __init__(self, mapping_file_dir, mapping_file_name):
        # self.pdb_id = pdb_id
        # self.mapping_file_dir = mapping_file_dir
        # self.mapping_file_name = mapping_file_name
        fp = open(os.path.join(mapping_file_dir, mapping_file_name))
        chain = mapping_file_name.strip().split(".")[0].strip().split("_")[1]
        self.pdb_to_fasta_map = {}
        self.fasta_to_pdb_map = {}
        for line in fp.readlines():
            pieces = line.strip().split("\t")
            pdb_ind = Chainindex.from_str_index(pieces[0].strip())
            # if "'" in pdb_ind:
            #     pdb_ind = pdb_ind[len(chain)+2:]
            # else:
            #     pdb_ind = pdb_ind[len(chain):]
            fasta_ind = pieces[1].strip()
            self.pdb_to_fasta_map[pdb_ind] = fasta_ind
            self.fasta_to_pdb_map[fasta_ind] = pdb_ind
        fp.close()
    def convert_PDBindx_To_FASTAindx(self, ind):
        if ind in self.pdb_to_fasta_map:
            return self.pdb_to_fasta_map[ind]
        return None
    def convert_FASTAindx_To_PDBindx(self, ind):
        if ind in self.fasta_to_pdb_map:
            return self.fasta_to_pdb_map[ind]
        return None

class Node:
    """loop format: A:1-2_3-4"""
    def __init__(self, chain, regions):
        self.chain = chain
        self.regions = regions

    def __eq__(self, other):
        return self.chain == other.chain and set(self.regions) == set(other.regions)

    def __ne__(self, other):
        return self.chain != other.chain or set(self.regions) != set(other.regions)

    def __lt__(self, other):
        if str(self) < str(other):
            return True
        False

    def __gt__(self, other):
        if str(self) > str(other):
            return True
        False

    def __repr__(self):
        int_regions = []
        for region in self.regions:
            s, e = region.strip().split("-")
            int_regions.append((int(s), int(e)))
        sorted_regions = []
        for s, e in sorted(int_regions):
            sorted_regions.append(str(s) + "-" + str(e))
            
        return self.chain+":"+"_".join(sorted_regions)
        # return self.chain+":"+"_".join(sorted(self.regions))

    def __hash__(self):
        return 11*hash(self.chain) + hash(frozenset(self.regions))

class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2

    def __eq__(self, other):
        # return hash(self.node1) ^ hash(self.node2) == hash(other.node1) ^ hash(other.node2)
        return (self.node1 == other.node1 and self.node2 == other.node2) or (self.node1 == other.node2 and self.node2 == other.node1)

    def __lt__(self, other):
        if self.node1 == other.node1:
            return self.node2 < other.node2
        else:
            return self.node1 < other.node1

    def __gt__(self, other):
        if self.node1 == other.node1:
            return self.node2 > other.node2
        else:
            return self.node1 > other.node1

    def __hash__(self):
        return hash(self.node1) + hash(self.node2)

    def __repr__(self):
        return str(self.node1)+"<->"+str(self.node2)