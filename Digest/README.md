# in silico trypsin digestion

Usage: python trypsin.py --input input_filename --output output_filename --miss 1

Note:

1. The script is written under Python 2.7
2. Biopython is prerequisite
3. The trypsin digestion script shared here follows proline rule, which means it does not cut lysine (K) or arginine (R) if they are followed by proline (P).
4. The program outputs a text format file which contains protein ID in the first column and corresponding tryptic peptides in the second column.

How to use:

1. Copy your fasta file and this script to same folder.
2. Open a command terminal and cd to this folder.
3. Type: python trypsin.py --input input_filename --output output_filename --miss 1

Three arguments are:

-input: name of your fasta file which contains protein sequences to be digested

-output output txt file name

-miss: number of allowed miss cleavage site, choose value from 0, 1, 2.

Speed

To digest Uniprot Human reference proteome (including canonical proteins and isoforms) which contains 88,717 proteins, the program finish within 22 seconds.
