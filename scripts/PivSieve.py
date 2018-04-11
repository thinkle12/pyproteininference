import os

import csv
sieve = []
with open('/Users/hinklet/random_analysis/k562_analysis/sieve/154787_sieve_lum1_comet.csv','r') as f:
    reader = csv.reader(f)
    for row in reader:
        sieve.append(row)

del sieve[0]

sieve_proteins = []
peptide_list = []
for info in sieve:
    if info[-2] not in peptide_list:
        sieve_proteins.append(info[6].split('|')[0])
        peptide_list.append(info[-2])


sieve_proteins_set = set(sieve_proteins).union()

sieve_proteins_list = list(sieve_proteins_set)

sieve_proteins_list = [x for x in sieve_proteins_list]

from collections import Counter

num_peps_sieve = Counter(sieve_proteins)

num_peps = [num_peps_sieve[x] for x in sieve_proteins_list]

one_hit_wonder_sieve = [x for x in sieve_proteins_list if num_peps_sieve[x]<=1]

non_1h1_sieve = [x for x in sieve_proteins_list if num_peps_sieve[x]>1]


pi = []
with open('/Users/hinklet/random_analysis/k562_analysis/pi_output/qvalues_leads_idwl_154787_pep.csv','r') as f:
    reader = csv.reader(f)
    for row in reader:
        pi.append(row)

del pi[0]
del pi[-1]

pi_proteins = [x[0].split('|')[0] for x in pi if float(x[2])<=.01]
pi_peptides = [len(x[6:]) for x in pi if float(x[2])<=.01]

one_hit_wonder_pi = [pi_proteins[x] for x in range(len(pi_peptides)) if pi_peptides[x]<=1]

non_1h1_pi = [pi_proteins[x] for x in range(len(pi_peptides)) if pi_peptides[x]>1]

print 'Number of Sieve Proteins Passing = '+str(len(sieve_proteins_list))
print 'Number of 1 Hit Wonder Sieve Proteins Passing = '+str(len(one_hit_wonder_sieve))
print 'Number of PI Proteins Passing = '+str(len(pi_proteins))
print 'Number of 1 Hit Wonder PI Proteins Passing = '+str(len(one_hit_wonder_pi))

non_1h1_intersect = list(set(non_1h1_pi)&set(non_1h1_sieve))

reg_intersect = list(set(pi_proteins)&set(sieve_proteins_list))

print 'Number of Intersecting Proteins = '+str(len(reg_intersect))
# print 'Number of Intersecting Non 1 Hit Wonder Proteins = '+str(len(non_1h1_intersect))
#
# one_hit_wonder_intersect = list(set(one_hit_wonder_sieve)&set(one_hit_wonder_pi))
#
# print 'Number of Intersecting 1 hit wonder proteins = '+str(len(one_hit_wonder_intersect))

from matplotlib_venn import venn2, venn2_circles
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

pp_ven = PdfPages('pi_vs_sieve_Venn_Diagram.pdf')

# Subset sizes
s = (
   len(pi_proteins)-len(reg_intersect),  # Ab
   len(sieve_proteins_list)-len(reg_intersect),  # aB
   len(reg_intersect),  # AB
)

v = venn2(subsets=s, set_labels=('PI    ', '    Sieve'))

# Subset labels
lbl1 = v.get_label_by_id('10')
lbl2 = v.get_label_by_id('01')
lbl3 = v.get_label_by_id('11')

lbl1.set_text('PI'+'\n'+'but'+'\n'+'not' +'\n'+'Sieve'+'\n'+str(len(pi_proteins)-len(reg_intersect)))
x1, y1 = lbl1.get_position()
lbl1.set_position((x1-.2, y1+0))

lbl2.set_text('Sieve'+'\n'+'but'+'\n'+'not'+'\n'+'PI'+'\n'+str(len(sieve_proteins_list)-len(reg_intersect)))
x2, y2 = lbl2.get_position()
lbl2.set_position((x1+1.25, y1+0))

lbl3.set_text('Sieve'+'\n'+'and'+'\n'+'PI'+'\n'+str(len(reg_intersect)))

# Subset colors
v.get_patch_by_id('10').set_color('#66ccff')
v.get_patch_by_id('01').set_color('#66b3ff')
v.get_patch_by_id('11').set_color('#ff9933')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)

#Border styles
c = venn2_circles(subsets=s, linestyle='solid')
plt.title('Venn Diagram Comparing Passing Proteins in Sieve and PI')
pp_ven.savefig()
plt.show()
plt.close()
pp_ven.close()


pp_ven2 = PdfPages('pi_vs_sieve_Venn_Diagram_non_1h1.pdf')

# Subset sizes
s = (
   len(non_1h1_pi)-len(non_1h1_intersect),  # Ab
   len(non_1h1_sieve)-len(non_1h1_intersect),  # aB
   len(non_1h1_intersect),  # AB
)

v = venn2(subsets=s, set_labels=('PI    ', '    Sieve'))

# Subset labels
lbl1 = v.get_label_by_id('10')
lbl2 = v.get_label_by_id('01')
lbl3 = v.get_label_by_id('11')

lbl1.set_text('PI'+'\n'+'but'+'\n'+'not' +'\n'+'Sieve'+'\n'+str(len(non_1h1_pi)-len(non_1h1_intersect)))
x1, y1 = lbl1.get_position()
lbl1.set_position((x1-.2, y1+0))

lbl2.set_text('Sieve'+'\n'+'but'+'\n'+'not'+'\n'+'PI'+'\n'+str(len(non_1h1_sieve)-len(non_1h1_intersect)))
x2, y2 = lbl2.get_position()
lbl2.set_position((x1+1.25, y1+0))

lbl3.set_text('Sieve'+'\n'+'and'+'\n'+'PI'+'\n'+str(len(non_1h1_intersect)))

# Subset colors
v.get_patch_by_id('10').set_color('#66ccff')
v.get_patch_by_id('01').set_color('#66b3ff')
v.get_patch_by_id('11').set_color('#ff9933')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)

#Border styles
c = venn2_circles(subsets=s, linestyle='solid')
plt.title('Venn Diagram Comparing Passing Proteins in Sieve and PI'+'\n'+'For Non One Hit Wonder Proteins')
pp_ven2.savefig()
plt.show()
plt.close()
pp_ven2.close()