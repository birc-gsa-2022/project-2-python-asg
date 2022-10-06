#############################################
# Suffixtree Naive Approach 

#############################################
# Libraries

import sys
from collections import deque
from cigar import edits_to_cigar
from align import get_edits

#############################################
# Classes

class Node(object):
    def __init__(self, start, end, suffix_order_leaf=None):
        self.start = start  # start_index of node (with respect to ref-string position).
        self.end = end  # end_index of node (with respect to ref-string position).
        self.out = {}  # dict containing children of node (named character corresponfing to start_index)
        self.suffix_order_leaf = suffix_order_leaf  # order of longest suffixes (longest=0)
        
#############################################
# Functions

def SuffixTree(string):
    ''' Function for building a suffix tree (naive approach). 
    Starts from left (beginning with longest suffix (whole string)). For each
    iteration follows path as far as possible and adds branch ($ if all path existed; 
    new 'long' node/subtree if first letter didnt exist in root node).
    
    Example:
    ref = 'abaaba'
    tree = SuffixTree(ref)
    read = 'aba'
    print([t for t in bf_order(tree)])
    #>>> [[None, None, None], [0, 1, None], [1, 3, None], [6, 7, 6], [3, 7, 2], [1, 3, None], [6, 7, 5], [6, 7, 4], [3, 7, 1], [6, 7, 3], [3, 7, 0]]
    '''
    if string == '' or string == None:
        return None

    string += '$'  # add sentinal to string.
    tree = Node(None,None)  # create root.
    string_length = len(string)
    
    count = 0  # enables tracking of longest suffix (since we iterate left->right each iteration will continuously add longest->shortest suffix).
    for i in range(string_length):  # loop through all suffixes.
        #print(count)
        current = tree  # set root as starting point.
        j = i # set i as startpoint for j and use k (see below) to increment j inside (nodes) while loop.
        while j < string_length:  # from the root walk down as far as possible.
            if string[j] in current.out:  # if a child contains 'first' letter, go that direction.
                child = current.out[string[j]]  
                length = child.end-child.start
                k = j+1  
                # if value of child contains multible letters go throug 
                # all and see if they match, if end of value is reach (node
                # gets exhausted) proceed to next child (if it exists).
                while k-j < length and string[k] == string[child.start+k-j]:
                    k += 1
                num = k-j
                if num == length:
                    current = child 
                    j = k
               # if node only contains some of the letters we split/branch it at last matching index.
                else:
                    start_edited = child.start + num
                    branch = Node(child.start, start_edited)  # create branch node.
                    branch.out[string[k]] = Node(k,string_length, count)  # add new child to branch.
                    branch.out[string[start_edited]] = child  # add existing child to branch.
                    child.start = start_edited  # edit start position of existing child.
                    current.out[string[j]] = branch  # replace existing node with new branched node.
                    count+=1
            # if no children contains first letter of suffix, node containing whole suffix is added (e.g. first step).
            else:  
                current.out[string[j]] = Node(j, string_length, count)
                count+=1
    return tree

def bf_order(tree):
    '''Breath-first traversal using queue.
    Returns list containing all Node() values [start, end, suffix_order_leaf].
    
    Example:
    tree = SuffixTree('abab')
    print([t for t in bf_order(tree)])
    #>>> [[None, None, None], [0, 2, None], [1, 2, None], [4, 5, 4], [4, 5, 2], [2, 5, 0], [4, 5, 3], [2, 5, 1]]
    
    '''
    if tree == '' or tree == None:
        return None

    queue = deque([tree])

    while queue:
        if type(queue[-1]) is type(Node(None,None,None)):
            tmp = queue.pop()
            queue.appendleft([tmp.start,tmp.end,tmp.suffix_order_leaf])
            for sub in tmp.out:
                queue.appendleft(tmp.out[sub])
        elif isinstance(queue[-1], list) == True:
            yield queue.pop()
        elif queue[-1] == None:
            queue.pop()
            
    
def match_seq(tree, ref, read):
    '''Matches if read/pattern exists in tree and returns below subtree.
    
    Example:
    ref = 'mississippi'
    tree = SuffixTree(ref)
    read = 'ss'
    subtree = match_seq(tree, ref, read)
    print([t for t in bf_order(subtree) if t[2] != None])
    #>>> [[8, 12, 5], [5, 12, 2]]
    '''
    if tree == '' or tree == None:
        return None
    if ref == '' or ref == None:
        return None
    if read == '' or read == None:
        return None

    current = tree
    seq=''
   
    i=0
    while i < len(read):
        # print(i, len(read))
        if read[i] not in current.out:
            return None
        child = current.out[read[i]]
        lab = ref[child.start:child.end]
        j=i+1
        while j-i < len(lab) and j < len(read) and read[j] == lab[j-i]:
            seq+=lab[j-i]
            j+=1
        if j-i == len(lab):
            current = child
            i = j  
        elif seq == read or j == len(read):
            return child
        else:
            return None
    else:
        return current


def read_fasta():
    # load input:
    inFile = sys.argv[1]
    with open(inFile,'r') as f:
        lines = f.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list


def read_fastq():
    inFile = sys.argv[2]  
    with open(inFile,'r') as f:
        lines = f.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('@'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list


################################################################
# Code:
    
# if __name__ == '__main__':

#     fasta_recs = read_fasta()
#     fastq_recs = read_fastq()
    
#     for fa_rec in fasta_recs:
#         ref = fa_rec[1]
#         tree = SuffixTree(ref)
#         for fq_rec in fastq_recs:
#             read = fq_rec[1]
#             subtree = match_seq(tree, ref, read)
#             matches = [t for t in bf_order(subtree) if t[2] != None]
#             for match in matches:
#                 match = match[2]
#                 read_name = fq_rec[0]
#                 read_seq = fq_rec[1]
#                 edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
#                 cigar = edits_to_cigar(edits[2])
#                 output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
#                 print('\t'.join(output))
          
        
################################################################

