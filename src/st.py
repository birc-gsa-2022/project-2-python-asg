#############################################
#############################################
#############################################
# Data structures


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
    def __init__(self, start, end, length, lable=None):
        self.start = start  # start_index of node (with respect to ref-string position).
        self.end = end  # end_index of node (with respect to ref-string position).
        self.length = length  # length of node  (end_index - start_index).
        self.out = {}  # dict containing children of node (named character corresponfing to start_index)
        self.lable = lable  # order of longest paths (longest=0)
        
#############################################
# Functions

def repString(string:str, start:int, end:int):
    '''Function for representing a string (so we dont
    have to store the whole thing).
    
    '''
    if string == '' or string == None:
        return []
    return string[start:end]


def SuffixTree(string):
    ''' Function for building a suffix tree (naive approach). 
    Starts from left (beginning with longest suffix = the whole string). For each
    iteration follows path as far as possible and adds branch ($ if all path existed; 
    new 'long' node/subtree if first letter didnt exist in root node).
    
    Example:
    ref = 'abaaba'
    tree = SuffixTree(ref)
    read = 'aba'
    print([t for t in bf_order(tree)])
    
    '''
    if string == '' or string == None:
        return []

    string += '$'  # add sentinal to string.
    tree = Node(None,None,None)  # create root.
    count = 0  # enables tracking of longest path (since we iterate left->right each iteration will continuously add longest->shortest path).
    for i in range(len(string)):  # loop through all suffixes.
        current = tree  # set root as starting point.
        j = i # set i as startpoint for j and use k (see below) to increment j inside (nodes) while loop.
        while j < len(string):  # from the root walk down as far as possible.
            if string[j] in current.out:  # if a child contains 'first' letter, go that direction.
                child = current.out[string[j]]  
                val = repString(string, child.start, child.end)
                k = j+1  
                # if value of child contains multible letters go throug 
                # all and see if they match, if end of value is reach (node
                # gets exhausted) proceed to next child (if it exists).
                while k-j < len(val) and string[k] == val[k-j]:
                    k += 1
                if k-j == len(val):
                    current = child 
                    j = k
               # if node only contains some of the letters we split/branch it at last matching index.
                else:   
                    branch = Node(child.start, child.start + k-j, child.start+k-j-child.start)  # create branch node.
                    branch.out[string[k]] = Node(k,len(string), len(string)-k, count)  # add new child to branch.
                    branch.out[val[k-j]] = child  # add existing child to branch.
                    child.start = child.start + k-j  # edit start position of existing child.
                    child.length = child.end - child.start  # edit length of existing child.
                    current.out[string[j]] = branch  # replace existing node with new branched node.
                    count+=1
            # if no children contains first letter of suffix, node containing whole suffix is added (e.g. first run/step).
            else:  
                current.out[string[j]] = Node(j, len(string), len(string)-j, count)
                count+=1
    return tree

def bf_order(tree):
    '''Breath-first traversal using queue.
    Returns list containing all Node() values [start, end, length, lable].
    
    Example:
    tree = SuffixTree('abab')
    print([t for t in bf_order(tree)])
    [[None, None, None, None], [0, 2, 2, None], [1, 2, 1, None], [4, 5, 1, 4],
     [4, 5, 1, 2], [2, 5, 3, 0], [4, 5, 1, 3], [2, 5, 3, 1]]
    
    '''
    if tree == '' or tree == None:
        return []

    queue = deque([tree])
    while queue:
        if type(queue[-1]) is type(Node(None,None,None)):
            tmp = queue.pop()
            queue.appendleft([tmp.start,tmp.end,tmp.length,tmp.lable])
            for sub in tmp.out:
                queue.appendleft(tmp.out[sub])
        elif isinstance(queue[-1], list) == True:
            yield queue.pop()
        elif queue[-1] == None:
            queue.pop()
            
    
def match_seq(tree, ref, read):
    '''Matches if read/pattern exists in tree and returns below subtree.
    
    Example:
    ref = 'abaaba'
    tree = SuffixTree(ref)
    read = 'aba'
    subtree = match_seq(tree, ref, read)
    print([t for t in bf_order(subtree)])
    '''
    if tree == '' or tree == None:
        return []
    if ref == '' or ref == None:
        return []
    if read == '' or read == None:
        return []

    current = tree
    seq=''
   
    i=0
    while i < len(read):
        if read[i] not in current.out:
            return None
        child = current.out[read[i]]
        lab = repString(ref, child.start, child.end)
        j=i+1
        while j-i < len(lab) and j < len(read) and read[j] == lab[j-i]:
            seq+=lab[j-i]
            j+=1
        if j-i == len(lab):
            current = child
            i = j  
        if seq == read or j == len(read):
            return child
    else:
        return None
    

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
    
if __name__ == '__main__':
    
    fasta_recs = read_fasta()
    fastq_recs = read_fastq()
    
    for fa_rec in fasta_recs:
        ref = fa_rec[1]
        tree = SuffixTree(ref)
        for fq_rec in fastq_recs:
            read = fq_rec[1]
            subtree = match_seq(tree, ref, read)
            matches = [t for t in bf_order(subtree) if t[3] != None]
            for match in matches:
                match = match[3]
                read_name = fq_rec[0]
                read_seq = fq_rec[1]
                edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
                cigar = edits_to_cigar(edits[2])
                output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
                print('\t'.join(output))
        
################################################################