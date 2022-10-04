################################################################
# libraries:
import time
import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

################################################################
# functions:

from st import SuffixTree
from st import bf_order
from st import match_seq
from naive import naive_algorithm
from SEQsimulator import simulate_string
from SEQsimulator import get_exact_read

def suffixtree_read_mapper(ref, read):
    tree = SuffixTree(ref)
    subtree = match_seq(tree, ref, read)
    matches = [t for t in bf_order(subtree) if t[3] != None]
    all_matches = [t[3] for t in matches]
    return all_matches

################################################################
# tests:
  
# Test suffixtree-algorithm vs naive-algorithm for same result:
for i in range(500000+1):
    print('Iteration nr: ', i+1)
    ref = simulate_string(random.randint(30,90))
    read = get_exact_read(ref, random.randint(1,20))
    if suffixtree_read_mapper(ref, read).sort() != naive_algorithm(ref,read).sort():
        print('Algorithm mistake!')
        break
    if i == 500000:
        print('DONE')


# Runtimes for the tree construction (varying ref lengths):
ref_lengths = [250,500,750,1000,1250,1500,1750]
runtimes = []
for idx in range(7):
    print('Iteration nr: ', idx+1) 
    replicate = []
    for j in range(10):
        ref = simulate_string(ref_lengths[idx])
        start_time = time.time()
        SuffixTree(ref)
        end_time = time.time()
        replicate.append(end_time-start_time)
    runtimes.append(np.mean(replicate))
# plot running times:
for i in range(len(runtimes)):
    runtimes[i] = runtimes[i]/ref_lengths[i]
fig, ax = plt.subplots()
sns.lineplot(x=ref_lengths, y=runtimes, ax=ax)
plt.xlabel('ref length')
plt.ylabel('runtime (s)')
plt.tight_layout()
plt.show()


# Runtimes for mapping (varying read lengths):
ref_lengths = [250,500,750,1000,1250]
read_lengths_10 = [10]*5
read_lengths_20 = [20]*5
read_lengths_30 = [30]*5
read_lengths_40 = [40]*5
read_lengths_50 = [50]*5
runtimes_10 = []
runtimes_20 = []
runtimes_30 = []
runtimes_40 = []
runtimes_50 = []
for idx in range(5):
    print('Iteration nr: ', idx+1)
    runtimes_10_replicate = []
    runtimes_20_replicate = []
    runtimes_30_replicate = []
    runtimes_40_replicate = []
    runtimes_50_replicate = []
    
    for i in range(10):
        ref = simulate_string(ref_lengths[idx])
        tree = SuffixTree(ref)
        
        # Dont know why it is nessesary to read run this chunk of code, but 
        # if i dont the first the first runtime is affected. Mayde something 
        # to do with loading the modules??
        read = get_exact_read(ref, 10)
        subtree = match_seq(tree, ref, read)
        matches = [t for t in bf_order(subtree) if t[3] != None]
        
        read = get_exact_read(ref, read_lengths_20[idx])
        start_time = time.time()
        subtree = match_seq(tree, ref, read)
        matches = [t for t in bf_order(subtree) if t[3] != None]
        end_time = time.time()
        runtimes_20_replicate.append(end_time-start_time)
        
        read = get_exact_read(ref, read_lengths_30[idx])
        start_time = time.time()
        subtree = match_seq(tree, ref, read)
        matches = [t for t in bf_order(subtree) if t[3] != None]
        end_time = time.time()
        runtimes_30_replicate.append(end_time-start_time)
        
        read = get_exact_read(ref, read_lengths_40[idx])
        start_time = time.time()
        subtree = match_seq(tree, ref, read)
        matches = [t for t in bf_order(subtree) if t[3] != None]
        end_time = time.time()
        runtimes_40_replicate.append(end_time-start_time)
        
        read = get_exact_read(ref, read_lengths_50[idx])
        start_time = time.time()
        subtree = match_seq(tree, ref, read)
        matches = [t for t in bf_order(subtree) if t[3] != None]
        end_time = time.time()
        runtimes_50_replicate.append(end_time-start_time)
        
        read = get_exact_read(ref, read_lengths_10[idx])
        start_time = time.time()
        subtree = match_seq(tree, ref, read)
        matches = [t for t in bf_order(subtree) if t[3] != None]
        end_time = time.time()
        runtimes_10_replicate.append(end_time-start_time)
        
    runtimes_10.append(np.mean(runtimes_10_replicate))
    runtimes_20.append(np.mean(runtimes_20_replicate))
    runtimes_30.append(np.mean(runtimes_30_replicate))
    runtimes_40.append(np.mean(runtimes_40_replicate))
    runtimes_50.append(np.mean(runtimes_50_replicate))

# plot running times:
fig, ax = plt.subplots()
sns.lineplot(x=ref_lengths, y=runtimes_50, ax=ax, label='read length = 50')
sns.lineplot(x=ref_lengths, y=runtimes_40, ax=ax, label='read length = 40')
sns.lineplot(x=ref_lengths, y=runtimes_30, ax=ax, label='read length = 30')
sns.lineplot(x=ref_lengths, y=runtimes_20, ax=ax, label='read length = 20')
sns.lineplot(x=ref_lengths, y=runtimes_10, ax=ax, label='read length = 10')
plt.xlabel('ref length')
plt.ylabel('runtime (s)')
plt.tight_layout()
plt.show()

