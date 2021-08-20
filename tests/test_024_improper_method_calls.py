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

import py_protein_inference
from py_protein_inference import in_silico_digest
from py_protein_inference.parameters import ProteinInferenceParameter
import os
import logging


TEST_DATABASE = resource_filename("py_protein_inference", "../tests/data/test_database.fasta")
TARGET_FILE = resource_filename("py_protein_inference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("py_protein_inference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("py_protein_inference", "../tests/data/test_params_peptide_centric.yaml")

OUTPUT_DIR = tempfile.gettempdir()


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("py_protein_inference.tests.test_024_improper_method_calls.py")


class TestImproperMethodCalls(TestCase):
    def test_incorrect_method_calls(self):

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

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
        pep_and_prot_data = py_protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
            append_alt_from_db=False,
        )

        # Validate that we will get a ValueError if we try to load the DataStore object WITHOUT running read_psms from Reader class
        with self.assertRaises(ValueError):
            data = py_protein_inference.datastore.DataStore(pep_and_prot_data, digest=digest)

        # Now read the psms
        pep_and_prot_data.read_psms()

        # Now load data object properly...
        data = py_protein_inference.datastore.DataStore(pep_and_prot_data, digest=digest)

        # Try to run restrict_psm_data with empty data...
        data.main_data_form = []
        with self.assertRaises(ValueError):
            data.restrict_psm_data()

        # Try to run create_scoring_input with empty data...
        with self.assertRaises(ValueError):
            data.create_scoring_input()

        # Reload data object properly again
        data = py_protein_inference.datastore.DataStore(pep_and_prot_data, digest=digest)

        # Run data restrict
        data.restrict_psm_data()

        # Try to run exclude_non_distinguishing_peptides without creating scoring input
        with self.assertRaises(ValueError):
            data.exclude_non_distinguishing_peptides()

        # Try to initialize a Score object without running create_scoring_input
        with self.assertRaises(ValueError):
            score = py_protein_inference.scoring.Score(data=data)

        # Finally intiate scoring input properly
        data.create_scoring_input()

        # Try to run picker without scored data...
        with self.assertRaises(ValueError):
            data.protein_picker()

        # Try to run inference without scored data...
        with self.assertRaises(ValueError):
            group = py_protein_inference.inference.PeptideCentric(data=data, digest=digest)

        with self.assertRaises(ValueError):
            group = py_protein_inference.inference.Parsimony(data=data, digest=digest)

        with self.assertRaises(ValueError):
            group = py_protein_inference.inference.Exclusion(data=data, digest=digest)

        with self.assertRaises(ValueError):
            group = py_protein_inference.inference.Inclusion(data=data, digest=digest)

        with self.assertRaises(ValueError):
            group = py_protein_inference.inference.FirstProtein(data=data, digest=digest)

        # Score psms properly
        score = py_protein_inference.scoring.Score(data=data)
        score.score_psms(score_method=protein_inference_parameters.protein_score)

        # Run picker properly
        data.protein_picker()

        # Try to run q value calculation without grouped scored proteins...
        with self.assertRaises(ValueError):
            data.calculate_q_values()

        # Run inference properly
        group = py_protein_inference.inference.PeptideCentric(data=data, digest=digest)
        group.infer_proteins()

        # Run q value calculation...
        data.calculate_q_values()

        # Export as normal
        export_type = protein_inference_parameters.export
        export = py_protein_inference.export.Export(data=data)
        export.export_to_csv(directory=os.path.join(OUTPUT_DIR, "leads"), export_type=export_type)
