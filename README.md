# ProteinInference

Usage: python Command_Line_PI_runner.py 


Note:

	1.	The script is written under Python 2.7.
	2.	A Local version of GLPK/GLPSOL is required.
	3.	Biopython, numpy, and matplotlib (if doing any plotting output) are necessary.
	4.	Tabbed Target and Decoy results from Percolator are required.
	5.	The Protein Inference code does Parsimony - Finds the minimal set of proteins given the peptides - and groups based on these lead proteins
	6.	FDR and Q value calculation is also based on the lead proteins of the groups.
	7.	Naming conventions are biased towards outputting Swissprot (reviewed) over Trembl (unreviewed)

How to use:

	1.	Download the ProteinInference package into a specific folder
	2.	Open a command terminal and cd to this folder.
	3.	Type: python Command_Line_PI_runner.py  -t example_percolator_target_psm.txt -d example_percolator_decoy_psm.txt -o example_output.csv -db example_db.fasta -mc 2 -dt trypsin -pl 7 -pep .9 -sm "multiplicative_log" -st "q_value" -ppk -gp "glpk" -fdr .01 -ex ["q_value_comma_sep"]

The 16 arguments are:

-t: name of the target percolator file which contains results from the search of interest

-d: name of the decoy percolator file which contains results from the search of interest

-db: name of your fasta file which contains protein sequences to be digested

-o: base name for the output file - name will be changed based on what parameters are chosen

-mc: number of missed cleavages that were allowed on the search from mascot or comet

-dt: the digestion type for the search i.e.: trypsin

-pl: the peptide length to restrict by i.e.: 7, meaning peptides of length less than 7 will not be used in PI

-pep: the posterior error probability value to restrict by i.e.: .9, meaning pep values greater than .9 will not be used in PI

-qv: the q value to restrict value by i.e.: .2, meaning q values greater than .2 will not be used in PI

-sm: the score method to choose, one of the following: 'best_peptide_per_protein', 'iterative_downweighted_log', 'multiplicative_log',' downweighted_multiplicative_log', 'downweighted_version2', 'top_two_combined','geometric_mean'

-st: the score type to use, can be either ‘q_value’ or ‘pep_value’

-ppk: whether or not to run protein picker, which is a target-decoy analysis used for removing lower scoring targets/decoys between pairs of like targets-decoys

-gp: how to do grouping can be one of the following: ‘simple_subsetting’, ‘glpk’; recommended to use ‘glpk’

-fdr: fdr value to use as a cutoff, only useful if export type is not a q value export type

-roc: filename of pdf if wanting to output a rock curve

-ex: a list of export types that are wanted i.e.: ["q_value_comma_sep”, “q_value_leads”, ”all”]

Speed

Running the full PI pipeline on a large K562 run from Lumos 1 finishes in about 5 minutes. Most of the time is spent doing the in silicon tryptic digestion of the search database.
