#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""
import csv
import tempfile
from unittest import TestCase
from pkg_resources import resource_filename

import protein_inference
from protein_inference import in_silico_digest
from protein_inference.parameters import ProteinInferenceParameter
import os

TEST_DATABASE = resource_filename('protein_inference', '../tests/data/test_database.fasta')
TARGET_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_target.txt')
DECOY_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_decoy.txt')
PARAMETER_FILE = resource_filename('protein_inference', '../tests/data/test_params_inclusion.yaml')
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')

LEAD_OUTPUT_FILE = resource_filename('protein_inference', '../tests/output/test_inclusion_q_value_leads_ml_posterior_error_prob.csv')
ALL_OUTPUT_FILE = resource_filename('protein_inference', '../tests/output/test_inclusion_q_value_all_ml_posterior_error_prob.csv')

IDENTIFIER_INDEX = 0
SCORE_INDEX = 1
Q_VALUE_INDEX = 2
GROUP_ID_INDEX = 5
PEPTIDES_INDEX = 6


class TestLoadInclusionWorkflow(TestCase):

    def test_workflow(self):

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        # TODO check to make sure proper parameters are being loaded
        self.assertEqual(protein_inference_parameters.digest_type, 'trypsin')
        self.assertEqual(protein_inference_parameters.export, 'q_value')
        self.assertEqual(protein_inference_parameters.fdr, 0.01)
        self.assertEqual(protein_inference_parameters.glpk_path, 'glpsol')
        self.assertEqual(protein_inference_parameters.missed_cleavages, 3)
        self.assertEqual(protein_inference_parameters.picker, True)
        self.assertEqual(protein_inference_parameters.restrict_pep, .9)
        self.assertEqual(protein_inference_parameters.restrict_peptide_length, 7)
        self.assertEqual(protein_inference_parameters.restrict_q, .9)
        self.assertEqual(protein_inference_parameters.score_method, 'multiplicative_log')
        self.assertEqual(protein_inference_parameters.score, 'posterior_error_prob')
        self.assertEqual(protein_inference_parameters.score_type, 'multiplicative')
        self.assertEqual(protein_inference_parameters.decoy_symbol, '##')
        self.assertEqual(protein_inference_parameters.isoform_symbol, '-')
        self.assertEqual(protein_inference_parameters.reviewed_identifier_symbol, 'sp|')
        self.assertEqual(protein_inference_parameters.inference_type, 'inclusion')
        self.assertEqual(protein_inference_parameters.tag, 'test_inclusion')
        self.assertEqual(protein_inference_parameters.grouping_type, 'subset_peptides')


        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = in_silico_digest.InSilicoDigest(database_path=TEST_DATABASE,parameter_file_object=protein_inference_parameters, id_splitting=True)
        digest.digest_fasta_database()

        # TODO dont test digest... test this in a separate test


        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data = protein_inference.reader.GenericReader(target_file=TARGET_FILE,
                                                                      decoy_file=DECOY_FILE,
                                                                      parameter_file_object=protein_inference_parameters,
                                                                      digest_class=digest)
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 16)


        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data = protein_inference.datastore.DataStore(pep_and_prot_data)

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters)

        self.assertEqual(len(data.main_data_restricted), 15)

        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        data.create_scoring_input(score_input = protein_inference_parameters.score)

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
            export.csv_export_leads_restricted(filename_out=OUTPUT_DIR + protein_inference_parameters.tag +'_' +'leads' +'_' + data.short_score_method +'_' + data.score + '.csv')
        if 'all' in export_type:
            export.csv_export_all_restricted(filename_out=OUTPUT_DIR + protein_inference_parameters.tag +'_' +'all' +'_' + data.short_score_method +'_' + data.score + '.csv')
        if 'comma_sep' in export_type:
            export.csv_export_comma_sep_restricted(filename_out=OUTPUT_DIR + protein_inference_parameters.tag +'_' +'comma_sep' +'_' + data.short_score_method +'_' + data.score + '.csv')
        if 'q_value_comma_sep' in export_type:
            export.csv_export_q_value_comma_sep(filename_out=OUTPUT_DIR + protein_inference_parameters.tag +'_' +'q_value_comma_sep' +'_' + data.short_score_method +'_' + data.score + '.csv')
        if 'q_value' in export_type:
            export.csv_export_q_value_leads(filename_out=OUTPUT_DIR + protein_inference_parameters.tag +'_' +'q_value_leads' +'_' + data.short_score_method +'_' + data.score + '.csv')
        if 'q_value' in export_type:
            export.csv_export_q_value_all(filename_out=OUTPUT_DIR + protein_inference_parameters.tag +'_' +'q_value_all' +'_' + data.short_score_method +'_' + data.score + '.csv')


        print('Protein Inference Finished')


        lead_output = []
        with open(LEAD_OUTPUT_FILE, 'r') as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=',')
            for row in reader:
                lead_output.append(row)

        del lead_output[0]

        protein_groups = data.protein_group_objects

        for i in range(len(protein_groups)):
            lead_protein = protein_groups[i].proteins[0]
            self.assertEqual(lead_protein.identifier, lead_output[i][IDENTIFIER_INDEX])
            self.assertAlmostEqual(lead_protein.score, float(lead_output[i][SCORE_INDEX]))
            self.assertEqual(protein_groups[i].q_value, float(lead_output[i][Q_VALUE_INDEX]))
            self.assertEqual(protein_groups[i].number_id, int(lead_output[i][GROUP_ID_INDEX]))
            self.assertEqual(lead_protein.peptides, set(lead_output[i][PEPTIDES_INDEX:]))

