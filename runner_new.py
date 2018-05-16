#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""


import ProteinInference
from Digest import insilicodigest
import os
import tempfile

# target_file = "data/159260_Bioplex2_b10090_percolator_target_psm.txt"
# decoy_file = "data/159260_Bioplex2_b10090_percolator_decoy_psm.txt"
# target_file = '/Users/hinklet/PythonPackages/PercolatorAnalysis/percolator_output/002402_percolator_target_psm_all_psms.txt'
# decoy_file = '/Users/hinklet/PythonPackages/PercolatorAnalysis/percolator_output/002402_percolator_decoy_psm_all_psms.txt'

data_dir = '/Users/hinklet/random_analysis/shigella_GP_p1/percolator_output/'

data_files = os.listdir(data_dir)

target_files = [os.path.join(data_dir,x) for x in data_files if 'target' in x]
decoy_files = [os.path.join(data_dir,x) for x in data_files if 'decoy' in x]

yaml_params = "parameters/Protein_Inference_Params.yaml"
database = "/Users/hinklet/random_analysis/shigella_GP_p1/Mouse_Shigella_up17_18_properformat.fasta"
output_dir = "/Users/hinklet/random_analysis/shigella_GP_p1/protein_inference_output/"

# if "Bioplex" in target_file:
#     tag = '_'.join(target_file.split('/')[-1].split('_')[:3])
#     print tag
# if "Bioplex" not in target_file:
#     tag = target_file.split('/')[-1].split('_')[0]
#     print tag

temp_dir = tempfile.gettempdir()
write_dir_input = temp_dir

#Initiate the reader...
#Input for now is a target percolator output and a decoy percolator output
pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file=target_files,
                                                          decoy_file=decoy_files)

#Execeute the reader instance, this loads the data into the reader class
pep_and_prot_data.execute()
#Next create a data store which is a class that stores all data for all steps of the PI process
#Each method and each class calls from this data class to gather information for analyses
data = ProteinInference.datastore.DataStore(pep_and_prot_data,yaml_param_file=yaml_params)

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

#Do in silico digest....
digest = insilicodigest.InSilicoDigest(database_path=database, num_miss_cleavs=int(data.yaml_params['Parameters']['Missed_Cleavages']), digest_type=data.yaml_params['Parameters']['Digest_Type'])
digest.execute()

try:
    os.mkdir('glpkinout/')
except OSError:
    pass

grouping_type = data.yaml_params['Parameters']['Group']

#Run simple group subsetting
if grouping_type=='simple_subsetting':
    group = ProteinInference.grouping.SimpleSubsetting(data_class=data)
    group.execute()

server_glpk_path = '/gne/research/apps/protchem/glpk/bin/glpsol'

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


# #Next run fdrcalc on the data....
# fdr = ProteinInference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=float(data.yaml_params['Parameters']['FDR']))
# fdr.execute()
q = ProteinInference.fdrcalc.QValueCalculation(data_class=data)
q.execute()
#Finally we have our output restricted data...
# restricted = data.fdr_restricted_grouped_scored_proteins
#Print the len of restricted data... which is how many protein groups pass FDR threshold
# print 'Number of Proteins passing an FDR of'+str(data.yaml_params['Parameters']['FDR'])+' = '+str(len(restricted))

export_type = data.yaml_params['Parameters']['Export']

#Write the output to a csv...
if 'leads' in export_type:
    export = ProteinInference.export.CsvOutLeads(data_class=data,filename_out=output_dir+'Custom'+'_'+'leads'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'all' in export_type:
    export = ProteinInference.export.CsvOutAll(data_class=data,filename_out=output_dir+'Custom'+'_'+'all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'comma_sep' in export_type:
    export = ProteinInference.export.CsvOutCommaSep(data_class=data, filename_out=output_dir+'Custom'+'_'+'comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value_comma_sep' in export_type:
    export = ProteinInference.export.CsvOutCommaSepQValues(data_class=data,filename_out=output_dir+'Custom'+'_'+'q_value_comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value' in export_type:
    export = ProteinInference.export.CsvOutLeadsQValues(data_class=data,filename_out=output_dir+'GP_p1'+'_'+'q_value_leads_other'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value' in export_type:
    export = ProteinInference.export.CsvOutAllQValues(data_class=data,filename_out=output_dir+'GP_p1'+'_'+'q_value_all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()


print 'Protein Inference Finished'
# if args.roc_curve:
#     roc = ProteinInference.benchmark.RocPlot(data_class=data)
#     roc.execute(pdf=args.roc_curve)


