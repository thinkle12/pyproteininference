# ProteinInference

Usage: python Command_Line_PI_runner.py --help


__Note:__

1. The script is written under Python 2.7.
2. A Local version of __GLPK/GLPSOL__ is required.
3. __Biopython__, __numpy__, and __matplotlib__ (if doing any plotting output) are necessary.
4. Tabbed Target and Decoy results from __Percolator__ are required.
5. The Protein Inference code does __Parsimony__ - Finds the minimal set of proteins given the peptides - and groups based on these lead proteins
6. FDR and Q value calculation is also based on the lead proteins of the groups.
7. Naming conventions are biased towards outputting Swissprot (reviewed) over Trembl (unreviewed)

__Command Line Usage:__

1. Download the ProteinInference package into a specific folder
2. Open a command terminal and cd to this folder.
3. Type: 
```bash
python Command_Line_PI_runner.py  -t example_percolator_target_psm.txt -d example_percolator_decoy_psm.txt -o example_output.csv -db example_db.fasta -mc 2 -dt trypsin -pl 7 -pep .9 -sm "multiplicative_log" -st "q_value" -ppk -gp "glpk" -fdr .01 -ex ["q_value_comma_sep"]
```


__The 16 arguments are:__

__-t:__ name of the target percolator file which contains results from the search of interest

__-d:__ name of the decoy percolator file which contains results from the search of interest

__-d:__ name of your fasta file which contains protein sequences to be digested

__-o:__ base name for the output file - name will be changed based on what parameters are chosen

__-mc:__ number of missed cleavages that were allowed on the search from mascot or comet

__-dt:__ the digestion type for the search i.e.: trypsin

__-pl:__ the peptide length to restrict by i.e.: 7, meaning peptides of length less than 7 will not be used in PI

__-pep:__ the posterior error probability value to restrict by i.e.: .9, meaning pep values greater than .9 will not be used in PI

__-qv:__ the q value to restrict value by i.e.: .2, meaning q values greater than .2 will not be used in PI

__-sm:__ the score method to choose, one of the following: 'best_peptide_per_protein', 'iterative_downweighted_log', 'multiplicative_log',' downweighted_multiplicative_log', 'downweighted_version2', 'top_two_combined','geometric_mean'

__-st:__ the score type to use, can be either ‘q_value’ or ‘pep_value’

__-ppk:__ whether or not to run protein picker, which is a target-decoy analysis used for removing lower scoring targets/decoys between pairs of like targets-decoys

__-gp:__ how to do grouping can be one of the following: ‘simple_subsetting’, ‘glpk’; recommended to use ‘glpk’

__-fdr:__ fdr value to use as a cutoff, only useful if export type is not a q value export type

__-roc:__ filename of pdf if wanting to output a rock curve

__-ex:__ a list of export types that are wanted i.e.: ["q_value_comma_sep”, “q_value_leads”, ”all”]

__Speed__

Running the full PI pipeline on a large K562 run from Lumos 1 finishes in about 5 minutes. Most of the time is spent doing the in silicon tryptic digestion of the search database.


#Example Python Runner

```python
import ProteinInference
from Digest import insilicodigest

#Do in silico trypsin digestion
digest = insilicodigest.InSilicoDigest(database_path='data/UniprotKBConcat1708_HUMAN.fasta',
                                       num_miss_cleavs=2,
                                       digest_type='trypsin')
digest.execute()

#Initiate the reader...
#Input for now is a target percolator output and a decoy percolator output
pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file='data/159260_Bioplex2_b10090_percolator_target_psm.txt',
                                                          decoy_file='/data/159260_Bioplex2_b10090_percolator_decoy_psm.txt')
#Execeute the reader instance, this loads the data into the reader class
pep_and_prot_data.execute()

#Next create a data store which is a class that stores all data for all steps of the PI process
#Each method and each class calls from this data class to gather information for analyses
data = ProteinInference.datastore.DataStore(pep_and_prot_data)

#Here restrict the data to having peptides with length 7 or greater and a pep of less than .9
restrict = ProteinInference.datastore.RestrictMainData(data, peptide_length=7, posterior_error_prob_threshold=.9,q_value_threshold=None)
restrict.execute()

#Here generate the pre score data using 'Q' values
score_setup = ProteinInference.datastore.PreScoreQValue(data)
score_setup.execute()

#Here we do scoring
score = ProteinInference.scoring.DownweightedMultiplicativeLog(data_class=data)
score.execute()

#Here we run protein picker
picker = ProteinInference.picker.StandardPicker(data_class=data)
picker.execute()

#Run GLPK to generate the minimal list of proteins that account for the peptides
#Running GLPK consists of 3 classes, setup, runner, and grouper which need to be run in succession
glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data,glpkin_filename='glpkinout/glpkout_example.mod')
glpksetup.execute()
runglpk = ProteinInference.grouping.GlpkRunner(path_to_glpsol = '/gne/research/apps/protchem/glpk/bin/glpsol',glpkin = 'glpkinout/glpkout_example.mod',glpkout = 'glpkinout/glpkout_example.sol',file_override = False)
runglpk.execute()
group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft', glpksolution_filename='glpkinout/glpkout_example.sol')
group.execute()

#Next run fdrcalc on the data....
fdr = ProteinInference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=.01)
fdr.execute()

#Also run Q value calculation
q = ProteinInference.fdrcalc.QValueCalculation(data_class=data)
q.execute()

#Write the output to a csv...
qval_out_csep = ProteinInference.export.CsvOutCommaSepQValues(data_class=data, filename_out='output/qvalues_csep_dwml_159260_Bioplex2_b10090_q_value.csv')
qval_out_csep.execute()

#Generate a Roc Plot
roc = ProteinInference.benchmark.RocPlot(data_class=data)
roc.execute(pdf='plots/dwml_159260_Bioplex2_b10090_q_value.pdf')

```