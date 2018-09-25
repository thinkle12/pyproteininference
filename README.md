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
3. Have the following files in current working directory: example_percolator_target_psm.txt, example_percolator_decoy_psm.txt, example_db.fasta, Protein_Inference_Params.yaml
4. The percolator files need to be generated from percolator - PercolatorAnalysis package will generate properly formatted files, fasta file should be from uniprot - fasta file needs to be the same or highly similar to the fasta file used to search the data via mascot/comet, yaml parameter file must be formatted properly like the example params.
5. Type: 
```bash
python Command_Line_PI_runner.py  -t example_percolator_target_psm.txt -d example_percolator_decoy_psm.txt -o output_dir/ -db example_db.fasta -ym Protein_Inference_Params.yaml
```


__The 6 arguments are:__

__-t:__ name of the target percolator file which contains results from the search of interest

__-d:__ name of the decoy percolator file which contains results from the search of interest

__-db:__ name of your fasta file which contains protein sequences to be digested

__-o:__ ProteinInference Result Directory to write to - Name of file will be determined by parameters selected in yaml file and searchID

__-roc:__ filename of pdf if wanting to output a rock curve

__-yaml:__ filename of the Protein Inference Yaml Parameter


__Speed__

Running the full PI pipeline on a large K562 run from Lumos 1 finishes in about 5 minutes. Most of the time is spent doing the in silicon tryptic digestion of the search database.


#Example Python Runner

```python
import ProteinInference
from Digest import insilicodigest
import os
import tempfile
import yaml


tag = 'yourtag'

target_files = 'percolator_target_output.txt'
decoy_files = 'percolator_decoy_output.txt'

yaml_params = "parameters/Protein_Inference_Params.yaml"

database = 'data/UniprotKBConcat1708_HUMAN.fasta'
output_dir = 'output/'

# Intermediate files writen to tmp
temp_dir = tempfile.gettempdir()
write_dir_input = temp_dir

with open(yaml_params, 'r') as stream:
    yaml_parameteres_for_digest = yaml.load(stream)

# Do in silico digest....
digest = insilicodigest.InSilicoDigest(database_path=database, num_miss_cleavs=int(yaml_parameteres_for_digest['Parameters']['Missed_Cleavages']), digest_type=yaml_parameteres_for_digest['Parameters']['Digest_Type'])
digest.execute()

#Initiate the reader...
#Input for now is a target percolator output and a decoy percolator output
pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file=target_files,
                                                          decoy_file=decoy_files,
                                                           digest_class=digest,
                                                           yaml_param_file=yaml_params)

#Execeute the reader instance, this loads the data into the reader class
pep_and_prot_data.execute()
#Next create a data store which is a class that stores all data for all steps of the PI process
#Each method and each class calls from this data class to gather information for analyses
data = ProteinInference.datastore.DataStore(pep_and_prot_data)

#Here restrict the data to having peptides with length 7 or greater
if data.yaml_params['Parameters']['Restrict_Pep']:
    pep_restrict = float(data.yaml_params['Parameters']['Restrict_Pep'])
else:
    pep_restrict = None
if data.yaml_params['Parameters']['Restrict_Q']:
    q_restrict = float(data.yaml_params['Parameters']['Restrict_Q'])
else:
    q_restrict = None

if data.yaml_params['Parameters']['Restrict_Peptide_Length']:
    pl_restrict = int(data.yaml_params['Parameters']['Restrict_Peptide_Length'])
else:
    pl_restrict = None


print 'restricting data'
restrict = ProteinInference.datastore.RestrictMainData(data,peptide_length=pl_restrict,posterior_error_prob_threshold=pep_restrict,q_value_threshold=q_restrict)
restrict.execute()

#Here generate the pre score data using 'PEP' or 'Q' values
if data.yaml_params['Parameters']['Score_Type'] == 'pep_value':
    score_setup = ProteinInference.datastore.PreScorePepValue(data)
if data.yaml_params['Parameters']['Score_Type'] == 'q_value':
    score_setup = ProteinInference.datastore.PreScoreQValue(data)

#Execute score setup...
score_setup.execute()

score_method = data.yaml_params['Parameters']['Score_Method']

#Here select scoring
if score_method=='best_peptide_per_protein':
    score = ProteinInference.scoring.BestPeptidePerProtein(data_class=data)
if score_method=='iterative_downweighted_log':
    score = ProteinInference.scoring.IterativeDownweightedLog(data_class=data)
if score_method=='multiplicative_log':
    score = ProteinInference.scoring.MultiplicativeLog(data_class=data)
if score_method=='downweighted_multiplicative_log':
    score = ProteinInference.scoring.DownweightedMultiplicativeLog(data_class=data)
if score_method=='downweighted_version2':
    score = ProteinInference.scoring.DownweightedVersion2(data_class=data)
if score_method=='top_two_combined':
    score = ProteinInference.scoring.TopTwoCombined(data_class=data)
if score_method=='geometric_mean':
    score = ProteinInference.scoring.GeometricMeanLog(data_class=data)

#Execute scoring...
score.execute()

#Run protein picker on the data
if data.yaml_params['Parameters']['Picker']:
    picker = ProteinInference.picker.StandardPicker(data_class=data)
    picker.execute()
else:
    pass

grouping_type = data.yaml_params['Parameters']['Group']

#Run simple group subsetting
if grouping_type=='simple_subsetting':
    group = ProteinInference.grouping.SimpleSubsetting(data_class=data)
    group.execute()

#Run GLPK setup, runner, grouper...
if grouping_type=='glpk':
    if grouping_type == 'glpk':
        glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data, glpkin_filename=os.path.join(write_dir_input,
                                                                                                      'glpkin_' + data.search_id + '.mod'))
        glpksetup.execute()
        glpkrun = ProteinInference.grouping.GlpkRunner(path_to_glpsol=data.yaml_params['Parameters']['GLPK_Path'],
                                                       glpkin=os.path.join(write_dir_input,
                                                                           'glpkin_' + data.search_id + '.mod'),
                                                       glpkout=os.path.join(write_dir_input,
                                                                            'glpkin_' + data.search_id + '.sol'),
                                                       file_override=False)
        glpkrun.execute()
        group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft',
                                                      glpksolution_filename=os.path.join(write_dir_input,
                                                                                         'glpkin_' + data.search_id + '.sol'))
        group.execute()

if grouping_type=='multi_subsetting':
    group = ProteinInference.grouping.MultiSubsetting(data_class=data)
    group.execute()


# fdr.execute()
q = ProteinInference.fdrcalc.QValueCalculation(data_class=data)
q.execute()

export_type = data.yaml_params['Parameters']['Export']

#Write the output to a csv...
if 'q_value' in export_type:
    export = ProteinInference.export.CsvOutLeadsQValues(data_class=data,filename_out=output_dir+tag+'_'+'q_value_leads_'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
```