import os

import csv
sieve = []
with open('/Users/hinklet/Downloads/134316_mascot_k562_sieve.csv','r') as f:
    reader = csv.reader(f)
    for row in reader:
        sieve.append(row)

del sieve[0]

sieve_proteins = []
peptide_list = []
for info in sieve:
    if info[-2] not in peptide_list:
        sieve_proteins.append(info[6])
        peptide_list.append(info[-2])


sieve_proteins_set = set(sieve_proteins).union()

sieve_proteins_list = list(sieve_proteins_set)

from collections import Counter

num_peps_sieve = Counter(sieve_proteins)

num_peps = [num_peps_sieve[x] for x in sieve_proteins_list]

one_hit_wonder_sieve = [x for x in sieve_proteins_list if num_peps_sieve[x]<=1]

non_1h1_sieve = [x for x in sieve_proteins_list if num_peps_sieve[x]>1]


pi = []
with open('output/qvalues_idwl_134316_pep.csv','r') as f:
    reader = csv.reader(f)
    for row in reader:
        pi.append(row)

del pi[0]
del pi[-1]

pi_proteins = [x[0] for x in pi if float(x[2])<=.01]
pi_peptides = [len(x[6:]) for x in pi]

one_hit_wonder_pi = [pi_proteins[x] for x in range(len(pi_peptides)) if pi_peptides[x]<=1]

non_1h1_pi = [pi_proteins[x] for x in range(len(pi_peptides)) if pi_peptides[x]>1]

print 'Number of Sieve Proteins Passing = '+str(len(sieve_proteins_list))
print 'Number of 1 Hit Wonder Sieve Proteins Passing = '+str(len(one_hit_wonder_sieve))
print 'Number of PI Proteins Passing = '+str(len(pi_proteins))
print 'Number of 1 Hit Wonder PI Proteins Passing = '+str(len(one_hit_wonder_pi))

non_1h1_intersect = list(set(non_1h1_pi)&set(non_1h1_sieve))

reg_intersect = list(set(pi_proteins)&set(sieve_proteins_list))

print 'Number of Intersecting Proteins = '+str(len(reg_intersect))
print 'Number of Intersecting Non 1 Hit Wonder Proteins = '+str(len(non_1h1_intersect))