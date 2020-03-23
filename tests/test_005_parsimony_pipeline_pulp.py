#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""
import csv
import tempfile
import unittest
from unittest import TestCase
from pkg_resources import resource_filename

import protein_inference
from protein_inference import in_silico_digest
from protein_inference.parameters import ProteinInferenceParameter
import os
import logging


TEST_DATABASE = resource_filename('protein_inference', '../tests/data/test_database.fasta')
TARGET_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_target.txt')
DECOY_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_decoy.txt')
PARAMETER_FILE = resource_filename('protein_inference', '../tests/data/test_params_parsimony_pulp.yaml')
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')
GLPKINOUT_PATH = resource_filename('protein_inference', '../tests/glpkinout/')

LEAD_OUTPUT_FILE = resource_filename('protein_inference', '../tests/output/test_parsimony_q_value_leads_ml_posterior_error_prob.csv')
ALL_OUTPUT_FILE = resource_filename('protein_inference', '../tests/output/test_parsimony_q_value_all_ml_posterior_error_prob.csv')

IDENTIFIER_INDEX = 0
SCORE_INDEX = 1
Q_VALUE_INDEX = 2
GROUP_ID_INDEX = 5
PEPTIDES_INDEX = 6

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("protein_inference.tests.test_001_parsimony_pipeline.py")


class TestLoadParsimonyPulpWorkflow(TestCase):

    @unittest.skip("Skipping Pulp Test, No CBC executable in build env")
    def test_workflow_parsimony_pulp(self):

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
        self.assertEqual(protein_inference_parameters.inference_type, 'parsimony')
        self.assertEqual(protein_inference_parameters.tag, 'test_parsimony')
        self.assertEqual(protein_inference_parameters.grouping_type, 'shared_peptides')
        self.assertEqual(protein_inference_parameters.max_identifiers_peptide_centric, 5)
        self.assertEqual(protein_inference_parameters.lp_solver, 'pulp')



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

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data = protein_inference.datastore.DataStore(pep_and_prot_data, digest_class=digest)

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters)

        self.assertEqual(len(data.main_data_restricted), 26)

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
            group.infer_proteins(glpkinout_directory=None, skip_running_glpk=None)

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
        logger.info('Number of Proteins passing an FDR of '+str(protein_inference_parameters.fdr)+' = '+str(len(data.fdr_restricted_grouped_scored_proteins)))


        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        export_type = protein_inference_parameters.export
        export = protein_inference.export.Export(data_class = data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)


        logger.info('Protein Inference Finished')

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
            # TODO also, make the pulp inference code cleaner... TY
            self.assertEqual(lead_protein.peptides, set(lead_output[i][PEPTIDES_INDEX:]))



