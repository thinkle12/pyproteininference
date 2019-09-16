#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""


import protein_inference
from digest import insilicodigest
import os
import tempfile
import yaml

# target_file = "data/159260_Bioplex2_b10090_percolator_target_psm.txt"
# decoy_file = "data/159260_Bioplex2_b10090_percolator_decoy_psm.txt"
# target_file = '/Users/hinklet/PythonPackages/PercolatorAnalysis/percolator_output/002402_percolator_target_psm_all_psms.txt'
# decoy_file = '/Users/hinklet/PythonPackages/PercolatorAnalysis/percolator_output/002402_percolator_decoy_psm_all_psms.txt'

tag = 'brd_combined'

data_dir = '/Users/hinklet/random_analysis/brd_analysis/percolator_output/'

data_files = os.listdir(data_dir)

target_files = [os.path.join(data_dir,x) for x in data_files if 'target' in x]
decoy_files = [os.path.join(data_dir,x) for x in data_files if 'decoy' in x]

# target_files = '159260_Bioplex2_b10090_percolator_target_psm.txt'
# decoy_files = '159260_Bioplex2_b10090_percolator_decoy_psm.txt'

yaml_params = "parameters/Protein_Inference_Params.yaml"
#database = "/Users/hinklet/random_analysis/shigella_GP_p1/Mouse_Shigella_up17_18_properformat.fasta"
#output_dir = "/Users/hinklet/random_analysis/shigella_GP_p1/protein_inference_output/"

database = 'data/UniprotKBConcat1708_HUMAN.fasta'
output_dir = '/Users/hinklet/random_analysis/brd_analysis/pi_output/'

# if "Bioplex" in target_file:
#     tag = '_'.join(target_file.split('/')[-1].split('_')[:3])
#     print tag
# if "Bioplex" not in target_file:
#     tag = target_file.split('/')[-1].split('_')[0]
#     print tag

temp_dir = tempfile.gettempdir()
write_dir_input = temp_dir

# Have to load the params here...
# Because we have to pass to in silico digest first...
# I really dont like this... but its essential given
# That we are going to start not storing all possible protein column...
with open(yaml_params, 'r') as stream:
    yaml_parameteres_for_digest = yaml.load(stream)

# Do in silico digest....
# Need digest before reader... ugh
digest = insilicodigest.PyteomicsDigest(database_path=database, num_miss_cleavs=int(yaml_parameteres_for_digest['Parameters']['Missed_Cleavages']), digest_type=yaml_parameteres_for_digest['Parameters']['Digest_Type'])
digest.execute()

#Initiate the reader...
#Input for now is a target percolator output and a decoy percolator output
pep_and_prot_data = protein_inference.reader.PercolatorRead(target_file=target_files,
                                                          decoy_file=decoy_files,
                                                           digest_class=digest,
                                                           yaml_param_file=yaml_params)

#Execeute the reader instance, this loads the data into the reader class
pep_and_prot_data.execute()
#Next create a data store which is a class that stores all data for all steps of the PI process
#Each method and each class calls from this data class to gather information for analyses
data = protein_inference.datastore.DataStore(pep_and_prot_data)

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
restrict = protein_inference.datastore.RestrictMainData(data,peptide_length=pl_restrict,posterior_error_prob_threshold=pep_restrict,q_value_threshold=q_restrict)
restrict.execute()

#Here generate the pre score data using 'PEP' or 'Q' values
if data.yaml_params['Parameters']['Score_Type'] == 'pep_value':
    score_setup = protein_inference.datastore.PreScorePepValue(data)
if data.yaml_params['Parameters']['Score_Type'] == 'q_value':
    score_setup = protein_inference.datastore.PreScoreQValue(data)

#Execute score setup...
score_setup.execute()

score_method = data.yaml_params['Parameters']['Score_Method']

#Here select scoring
if score_method=='best_peptide_per_protein':
    score = protein_inference.scoring.BestPeptidePerProtein(data_class=data)
if score_method=='iterative_downweighted_log':
    score = protein_inference.scoring.IterativeDownweightedLog(data_class=data)
if score_method=='multiplicative_log':
    score = protein_inference.scoring.MultiplicativeLog(data_class=data)
if score_method=='downweighted_multiplicative_log':
    score = protein_inference.scoring.DownweightedMultiplicativeLog(data_class=data)
if score_method=='downweighted_version2':
    score = protein_inference.scoring.DownweightedVersion2(data_class=data)
if score_method=='top_two_combined':
    score = protein_inference.scoring.TopTwoCombined(data_class=data)
if score_method=='geometric_mean':
    score = protein_inference.scoring.GeometricMeanLog(data_class=data)

#Execute scoring...
score.execute()

#Run protein picker on the data
if data.yaml_params['Parameters']['Picker']:
    picker = protein_inference.picker.StandardPicker(data_class=data)
    picker.execute()
else:
    pass

try:
    os.mkdir('glpkinout/')
except OSError:
    pass

grouping_type = data.yaml_params['Parameters']['Group']

#Run simple group subsetting
if grouping_type=='simple_subsetting':
    group = protein_inference.grouping.SimpleSubsetting(data_class=data)
    group.execute()

server_glpk_path = '/gne/research/apps/protchem/glpk/bin/glpsol'

#Run GLPK setup, runner, grouper...
if grouping_type=='glpk':
    if grouping_type == 'glpk':
        glpksetup = protein_inference.grouping.GlpkSetup(data_class=data, glpkin_filename=os.path.join(write_dir_input,
                                                                                                      'glpkin_' + data.search_id + '.mod'))
        glpksetup.execute()
        glpkrun = protein_inference.grouping.GlpkRunner(path_to_glpsol=data.yaml_params['Parameters']['GLPK_Path'],
                                                       glpkin=os.path.join(write_dir_input,
                                                                           'glpkin_' + data.search_id + '.mod'),
                                                       glpkout=os.path.join(write_dir_input,
                                                                            'glpkin_' + data.search_id + '.sol'),
                                                       file_override=False)
        glpkrun.execute()
        group = protein_inference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft',
                                                      glpksolution_filename=os.path.join(write_dir_input,
                                                                                         'glpkin_' + data.search_id + '.sol'))
        group.execute()

if grouping_type=='multi_subsetting':
    group = protein_inference.grouping.MultiSubsetting(data_class=data)
    group.execute()


# #Next run fdrcalc on the data....
# fdr = protein_inference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=float(data.yaml_params['Parameters']['FDR']))
# fdr.execute()
q = protein_inference.fdrcalc.QValueCalculation(data_class=data)
q.execute()
#Finally we have our output restricted data...
# restricted = data.fdr_restricted_grouped_scored_proteins
#Print the len of restricted data... which is how many protein groups pass FDR threshold
# print 'Number of Proteins passing an FDR of'+str(data.yaml_params['Parameters']['FDR'])+' = '+str(len(restricted))

export_type = data.yaml_params['Parameters']['Export']

#Write the output to a csv...
if 'leads' in export_type:
    export = protein_inference.export.CsvOutLeads(data_class=data,filename_out=output_dir+'Custom'+'_'+'leads'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'all' in export_type:
    export = protein_inference.export.CsvOutAll(data_class=data,filename_out=output_dir+'Custom'+'_'+'all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'comma_sep' in export_type:
    export = protein_inference.export.CsvOutCommaSep(data_class=data, filename_out=output_dir+'Custom'+'_'+'comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value_comma_sep' in export_type:
    export = protein_inference.export.CsvOutCommaSepQValues(data_class=data,filename_out=output_dir+'Custom'+'_'+'q_value_comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value' in export_type:
    export = protein_inference.export.CsvOutLeadsQValues(data_class=data,filename_out=output_dir+tag+'_'+'q_value_leads_other'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value' in export_type:
    export = protein_inference.export.CsvOutAllQValues(data_class=data,filename_out=output_dir+tag+'_'+'q_value_all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value' in export_type:
    export = protein_inference.export.CsvOutLeadsQValuesLong(data_class=data,filename_out=output_dir + tag + '_long_' + 'q_value_all' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
    export.execute()


print 'Protein Inference Finished'
# if args.roc_curve:
#     roc = protein_inference.benchmark.RocPlot(data_class=data)
#     roc.execute(pdf=args.roc_curve)


