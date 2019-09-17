#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""

import protein_inference
from protein_inference import in_silico_digest
from protein_inference.parameters import ProteinInferenceParameter
import os
import argparse


parser = argparse.ArgumentParser(description='Protein Inference')
parser.add_argument("-t","--target", dest="target", required=True,
                    help="Input target psm output from percolator", metavar="FILE")
parser.add_argument("-d","--decoy", dest="decoy", required=True,
                    help="Input decoy psm output from percolator", metavar="FILE")
parser.add_argument("-o","--output", dest="dir_name", required=True,
                    help="protein_inference Result Directory to write to - Name of file will be determined by parameters selected and searchID", metavar="FILE")
parser.add_argument("-db","--database", dest="database", required=True,
                    help="Provide the database used in the MS search", metavar="FILE")
parser.add_argument("-roc","--roc_curve_filename", dest="roc_curve", required=False,
                    help="Provide the a filename with .pdf extension for roc curve output", metavar="FILE")
parser.add_argument("-ym","--yaml_params", dest="yaml_params", required=False,
                    help="Provide a Protein Inference Yaml Parameter File... If none given, default parameters will be ran", metavar="FILE")
args = parser.parse_args()

### STEP 1: Load parameter file ###
### STEP 1: Load parameter file ###
### STEP 1: Load parameter file ###
protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=args.yaml_params)


### STEP 2: Start with running an In Silico Digestion ###
### STEP 2: Start with running an In Silico Digestion ###
### STEP 2: Start with running an In Silico Digestion ###
digest = in_silico_digest.InSilicoDigest(database_path=args.database,parameter_file_object=protein_inference_parameters, id_splitting=True)
digest.digest_fasta_database()


### STEP 3: Read PSM Data ###
### STEP 3: Read PSM Data ###
### STEP 3: Read PSM Data ###
pep_and_prot_data = protein_inference.reader.PercolatorReader(target_file=args.target,
                                                              decoy_file=args.decoy,
                                                              parameter_file_object=protein_inference_parameters,
                                                              digest_class=digest)
pep_and_prot_data.read_psms()

### STEP 4: Initiate the datastore class ###
### STEP 4: Initiate the datastore class ###
### STEP 4: Initiate the datastore class ###
data = protein_inference.datastore.DataStore(pep_and_prot_data)

### Step 5: Restrict the PSM data
### Step 5: Restrict the PSM data
### Step 5: Restrict the PSM data
data.restrict_psm_data(parameter_file_object=protein_inference_parameters)

### Step 6: Generate protein scoring input
### Step 6: Generate protein scoring input
### Step 6: Generate protein scoring input
data.create_scoring_input(score_input = protein_inference_parameters.score_type)

### Step 7: Remove non unique peptides if running exclusion
### Step 7: Remove non unique peptides if running exclusion
### Step 7: Remove non unique peptides if running exclusion
if protein_inference_parameters.inference_type=="exclusion":
    # This gets ran if we run exclusion...
    data.exclude_non_distinguishing_peptides(digest_class=digest)


### STEP 8: Score our PSMs given a score method
### STEP 8: Score our PSMs given a score method
### STEP 8: Score our PSMs given a score method
score = protein_inference.scoring.Score(data_class=data)
score.score_psms(score_method=protein_inference_parameters.score_method)


### STEP 9: Run protein picker on the data
### STEP 9: Run protein picker on the data
### STEP 9: Run protein picker on the data
if protein_inference_parameters.picker:
    data.protein_picker()
else:
    pass


### STEP 10: Apply Inference
### STEP 10: Apply Inference
### STEP 10: Apply Inference
inference_type = protein_inference_parameters.inference_type

# For parsimony... Run GLPK setup, runner, grouper...
if inference_type == 'parsimony':
    group = protein_inference.inference.Parsimony(data_class=data, digest_class=digest)
    group.infer_proteins()

if inference_type == "inclusion":
    group = protein_inference.inference.Inclusion(data_class=data, digest_class=digest)
    group.infer_proteins()

if inference_type == "exclusion":
    group = protein_inference.inference.Exclusion(data_class=data, digest_class=digest)
    group.infer_proteins()

### STEP 11: Run FDR and Q value Calculations
### STEP 11: Run FDR and Q value Calculations
### STEP 11: Run FDR and Q value Calculations
data.set_based_fdr(false_discovery_rate=float(protein_inference_parameters.fdr))
data.calculate_q_values()

#Print the len of restricted data... which is how many protein groups pass FDR threshold
print('Number of Proteins passing an FDR of'+str(protein_inference_parameters.fdr)+' = '+str(len(data.fdr_restricted_grouped_scored_proteins)))


### STEP 12: Export to CSV
### STEP 12: Export to CSV
### STEP 12: Export to CSV
export_type = protein_inference_parameters.export
export = protein_inference.export.Export(data_class = data)
if 'leads' in export_type:
    export.csv_export_leads_restricted(filename_out=args.dir_name+protein_inference_parameters.tag+'_'+'leads'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
if 'all' in export_type:
    export.csv_export_all_restricted(filename_out=args.dir_name+protein_inference_parameters.tag+'_'+'all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
if 'comma_sep' in export_type:
    export.csv_export_comma_sep_restricted(filename_out=args.dir_name+protein_inference_parameters.tag+'_'+'comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
if 'q_value_comma_sep' in export_type:
    export.csv_export_q_value_comma_sep(filename_out=args.dir_name+protein_inference_parameters.tag+'_'+'q_value_comma_sep'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
if 'q_value' in export_type:
    export.csv_export_q_value_leads(filename_out=args.dir_name+protein_inference_parameters.tag+'_'+'q_value_leads'+'_'+data.short_score_method+'_'+data.score_type+'.csv')
if 'q_value_all' in export_type:
    export.csv_export_q_value_all(filename_out=args.dir_name+protein_inference_parameters.tag+'_'+'q_value_all'+'_'+data.short_score_method+'_'+data.score_type+'.csv')


print('Protein Inference Finished')





