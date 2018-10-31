#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""

import py_protein_inference
import argparse
from digest import insilicodigest
import os
import yaml


parser = argparse.ArgumentParser(description='Protein Inference')
parser.add_argument("-t","--target", dest="target", required=True,
                    help="Input target psm output from percolator", metavar="FILE")
parser.add_argument("-d","--decoy", dest="decoy", required=True,
                    help="Input decoy psm output from percolator", metavar="FILE")
parser.add_argument("-o","--output", dest="dir_name", required=True,
                    help="py_protein_inference Result Directory to write to - Name of file will be determined by parameters selected and searchID", metavar="FILE")
parser.add_argument("-db","--database", dest="database", required=True,
                    help="Provide the database used in the MS search", metavar="FILE")
parser.add_argument("-roc","--roc_curve_filename", dest="roc_curve", required=False,
                    help="Provide the a filename with .pdf extension for roc curve output", metavar="FILE")
parser.add_argument("-ym","--yaml_params", dest="yaml_params", required=False,
                    help="Provide a Protein Inference Yaml Parameter File... If none given, default parameters will be ran", metavar="FILE")
args = parser.parse_args()

if "Bioplex" in args.target:
    tag = '_'.join(args.target.split('/')[-1].split('_')[:3])
    print tag
if "Bioplex" not in args.target:
    tag = args.target.split('/')[-1].split('_')[0]
    print tag


# Need to read in yaml params here for digest...
# Digest needs to be done first to determine alternative proteins
# Need to determine all potentially possible proteins because comet loader
# Now does not load all possible proteins
with open(args.yaml_params, 'r') as stream:
    yaml_parameteres_for_digest = yaml.load(stream)


# Do in silico digest....
digest = insilicodigest.InSilicoDigest(database_path=args.database,
                                       num_miss_cleavs=int(yaml_parameteres_for_digest['Parameters']['Missed_Cleavages']),
                                       digest_type=yaml_parameteres_for_digest['Parameters']['Digest_Type'])
digest.execute()

#Initiate the reader...
#Input for now is a target percolator output and a decoy percolator output
pep_and_prot_data = py_protein_inference.reader.PercolatorRead(target_file=args.target,
                                                          decoy_file=args.decoy,
                                                           yaml_param_file=args.yaml_params,
                                                           digest_class=digest)

#Execeute the reader instance, this loads the data into the reader class
pep_and_prot_data.execute()
#Next create a data store which is a class that stores all data for all steps of the PI process
#Each method and each class calls from this data class to gather information and store data for analyses
data = py_protein_inference.datastore.DataStore(pep_and_prot_data)

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
restrict = py_protein_inference.datastore.RestrictMainData(data,peptide_length=pl_restrict,posterior_error_prob_threshold=pep_restrict,q_value_threshold=q_restrict)
restrict.execute()

#Here generate the pre score data using 'PEP' or 'Q' values
if data.yaml_params['Parameters']['Score_Type'] == 'pep_value':
    score_setup = py_protein_inference.datastore.PreScorePepValue(data)
if data.yaml_params['Parameters']['Score_Type'] == 'q_value':
    score_setup = py_protein_inference.datastore.PreScoreQValue(data)

#Execute score setup...
score_setup.execute()

score_method = data.yaml_params['Parameters']['Score_Method']

#Here select scoring
if score_method=='best_peptide_per_protein':
    score = py_protein_inference.scoring.BestPeptidePerProtein(data_class=data)
if score_method=='iterative_downweighted_log':
    score = py_protein_inference.scoring.IterativeDownweightedLog(data_class=data)
if score_method=='multiplicative_log':
    score = py_protein_inference.scoring.MultiplicativeLog(data_class=data)
if score_method=='downweighted_multiplicative_log':
    score = py_protein_inference.scoring.DownweightedMultiplicativeLog(data_class=data)
if score_method=='downweighted_version2':
    score = py_protein_inference.scoring.DownweightedVersion2(data_class=data)
if score_method=='top_two_combined':
    score = py_protein_inference.scoring.TopTwoCombined(data_class=data)
if score_method=='geometric_mean':
    score = py_protein_inference.scoring.GeometricMeanLog(data_class=data)

#Execute scoring...
score.execute()

#Run protein picker on the data
if data.yaml_params['Parameters']['Picker']:
    picker = py_protein_inference.picker.StandardPicker(data_class=data)
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
    group = py_protein_inference.grouping.SimpleSubsetting(data_class=data)
    group.execute()

server_glpk_path = '/gne/research/apps/protchem/glpk/bin/glpsol'

#Run GLPK setup, runner, grouper...
if grouping_type=='glpk':
    glpksetup = py_protein_inference.grouping.GlpkSetup(data_class=data,glpkin_filename='glpkinout/glpkin_'+tag+'.mod')
    glpksetup.execute()
    glpkrun = py_protein_inference.grouping.GlpkRunner(path_to_glpsol = data.yaml_params['Parameters']['GLPK_Path'],glpkin='glpkinout/glpkin_'+tag+'.mod',glpkout='glpkinout/glpkout_'+tag+'.sol',file_override=False)
    glpkrun.execute()
    group = py_protein_inference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft', glpksolution_filename='glpkinout/glpkout_'+tag+'.sol')
    group.execute()

if grouping_type=='multi_subsetting':
    group = py_protein_inference.grouping.MultiSubsetting(data_class=data)
    group.execute()


#Next run fdrcalc on the data....
fdr = py_protein_inference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=float(data.yaml_params['Parameters']['FDR']))
fdr.execute()
q = py_protein_inference.fdrcalc.QValueCalculation(data_class=data)
q.execute()
#Finally we have our output restricted data...
restricted = data.fdr_restricted_grouped_scored_proteins
#Print the len of restricted data... which is how many protein groups pass FDR threshold
print 'Number of Proteins passing an FDR of'+str(data.yaml_params['Parameters']['FDR'])+' = '+str(len(restricted))

export_type = data.yaml_params['Parameters']['Export']

#Write the output to a csv...
if 'leads' in export_type:
    export = py_protein_inference.export.CsvOutLeads(data_class=data,filename_out=args.dir_name+tag+'_'+'leads'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'all' in export_type:
    export = py_protein_inference.export.CsvOutAll(data_class=data,filename_out=args.dir_name+tag+'_'+'all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'comma_sep' in export_type:
    export = py_protein_inference.export.CsvOutCommaSep(data_class=data, filename_out=args.dir_name+tag+'_'+'comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value_comma_sep' in export_type:
    export = py_protein_inference.export.CsvOutCommaSepQValues(data_class=data,filename_out=args.dir_name+tag+'_'+'q_value_comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value_leads' in export_type:
    export = py_protein_inference.export.CsvOutLeadsQValues(data_class=data,filename_out=args.dir_name+tag+'_'+'q_value_leads'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()
if 'q_value_all' in export_type:
    export = py_protein_inference.export.CsvOutAllQValues(data_class=data,filename_out=args.dir_name+tag+'_'+'q_value_all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
    export.execute()


print 'Protein Inference Finished' 
# if args.roc_curve:
#     roc = ProteinInference.benchmark.RocPlot(data_class=data)
#     roc.execute(pdf=args.roc_curve)


