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
import logging


TEST_DATABASE = resource_filename(
    "protein_inference", "../tests/data/test_database.fasta"
)
TARGET_FILE = resource_filename(
    "protein_inference", "../tests/data/test_perc_data_target.txt"
)
DECOY_FILE = resource_filename(
    "protein_inference", "../tests/data/test_perc_data_decoy.txt"
)
PARAMETER_FILE = resource_filename(
    "protein_inference", "../tests/data/test_params_parsimony_glpk.yaml"
)
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')
for sub_dir in ["leads", "all", "peptides", "psms", "psm_ids"]:
    if not os.path.exists(os.path.join(OUTPUT_DIR, sub_dir)):
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))

GLPKINOUT_PATH = resource_filename("protein_inference", "../tests/glpkinout/")
SKIP_RUNNING_GLPK = True

IDENTIFIER_INDEX = 0
SCORE_INDEX = 1
Q_VALUE_INDEX = 2
GROUP_ID_INDEX = 5
PEPTIDES_INDEX = 6

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("protein_inference.tests.test_020_parsimony_shared_peptides_assignment.py")


class TestLoadParsimonyGlpkWorkflowSharedPeptideReassignment(TestCase):

    def test_workflow_parsimony_glpk_shared_peptides_best(self):

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        # Set the shared_peptides param to "best"
        protein_inference_parameters.shared_peptides = "best"


        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = in_silico_digest.InSilicoDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data = protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest_class=digest,
            append_alt_from_db=False,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data = protein_inference.datastore.DataStore(
            pep_and_prot_data, digest_class=digest
        )

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data.restrict_psm_data()

        self.assertEqual(len(data.main_data_restricted), 26)

        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        data.create_scoring_input()

        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters.inference_type == protein_inference.inference.Inference.EXCLUSION:
            # This gets ran if we run exclusion...
            data.exclude_non_distinguishing_peptides()

        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        score = protein_inference.scoring.Score(data_class=data)
        score.score_psms(score_method=protein_inference_parameters.protein_score)

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
        if inference_type == protein_inference.inference.Inference.PARSIMONY:
            group = protein_inference.inference.Parsimony(
                data_class=data, digest_class=digest
            )
            group.infer_proteins(
                glpkinout_directory=GLPKINOUT_PATH, skip_running_glpk=SKIP_RUNNING_GLPK
            )

        if inference_type == protein_inference.inference.Inference.INCLUSION:
            group = protein_inference.inference.Inclusion(
                data_class=data, digest_class=digest
            )
            group.infer_proteins()

        if inference_type == protein_inference.inference.Inference.EXCLUSION:
            group = protein_inference.inference.Exclusion(
                data_class=data, digest_class=digest
            )
            group.infer_proteins()

        lead_protein_psms = [x[0].psms for x in data.grouped_scored_proteins]

        # Get the PSM Identifiers from our lead proteins
        psm_id_list = []
        for protein in lead_protein_psms:
            for psm in protein:
                psm_id_list.append(psm.identifier)

        # Make sure the there are no duplicate psms in our lead proteins
        self.assertEqual(len(psm_id_list), len(set(psm_id_list)))

    def test_workflow_parsimony_glpk_shared_peptides_all(self):

        ## Do the same analysis but put shared peptides everywhere

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters_all = ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        # Set the shared_peptides param to "best"
        protein_inference_parameters_all.shared_peptides = "all"


        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest_all = in_silico_digest.InSilicoDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters_all.digest_type,
            missed_cleavages=protein_inference_parameters_all.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters_all.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters_all.restrict_peptide_length,
            id_splitting=True,
        )
        digest_all.digest_fasta_database()

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data_all = protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters_all,
            digest_class=digest_all,
            append_alt_from_db=False,
        )
        pep_and_prot_data_all.read_psms()

        self.assertEqual(len(pep_and_prot_data_all.psms), 27)

        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data_all = protein_inference.datastore.DataStore(
            pep_and_prot_data_all, digest_class=digest_all
        )

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data_all.restrict_psm_data()

        self.assertEqual(len(data_all.main_data_restricted), 26)

        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        data_all.create_scoring_input()

        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters_all.inference_type == protein_inference.inference.Inference.EXCLUSION:
            # This gets ran if we run exclusion...
            data_all.exclude_non_distinguishing_peptides()

        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        score_all = protein_inference.scoring.Score(data_class=data_all)
        score_all.score_psms(score_method=protein_inference_parameters_all.protein_score)

        ### STEP 9: Run protein picker on the data
        ### STEP 9: Run protein picker on the data
        ### STEP 9: Run protein picker on the data
        if protein_inference_parameters_all.picker:
            data_all.protein_picker()
        else:
            pass

        ### STEP 10: Apply Inference
        ### STEP 10: Apply Inference
        ### STEP 10: Apply Inference
        inference_type = protein_inference_parameters_all.inference_type

        # For parsimony... Run GLPK setup, runner, grouper...
        if inference_type == protein_inference.inference.Inference.PARSIMONY:
            group = protein_inference.inference.Parsimony(
                data_class=data_all, digest_class=digest_all
            )
            group.infer_proteins(
                glpkinout_directory=GLPKINOUT_PATH, skip_running_glpk=SKIP_RUNNING_GLPK
            )

        if inference_type == protein_inference.inference.Inference.INCLUSION:
            group = protein_inference.inference.Inclusion(
                data_class=data_all, digest_class=digest_all
            )
            group.infer_proteins()

        if inference_type == protein_inference.inference.Inference.EXCLUSION:
            group = protein_inference.inference.Exclusion(
                data_class=data_all, digest_class=digest_all
            )
            group.infer_proteins()


        lead_protein_psms = [x[0].psms for x in data_all.grouped_scored_proteins]

        # Get the PSM Identifiers from our lead proteins
        psm_id_list = []
        for protein in lead_protein_psms:
            for psm in protein:
                psm_id_list.append(psm.identifier)

        # Make sure the there are duplicate psms in our lead protein objects
        self.assertNotEqual(len(psm_id_list), len(set(psm_id_list)))
