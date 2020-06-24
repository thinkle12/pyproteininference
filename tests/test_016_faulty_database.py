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

TEST_DATABASE_ORIGINAL = resource_filename(
    "protein_inference", "../tests/data/test_database.fasta"
)
TEST_DATABASE_MISSING_PEP = resource_filename(
    "protein_inference", "../tests/data/test_database_missing_peptide.fasta"
)
TEST_DATABASE_MISSING_PROT = resource_filename(
    "protein_inference", "../tests/data/test_database_missing_protein.fasta"
)
TARGET_FILE = resource_filename(
    "protein_inference", "../tests/data/test_perc_data_target.txt"
)
DECOY_FILE = resource_filename(
    "protein_inference", "../tests/data/test_perc_data_decoy.txt"
)
PARAMETER_FILE = resource_filename(
    "protein_inference", "../tests/data/test_params_inclusion.yaml"
)
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')
for sub_dir in ["leads", "all", "peptides", "psms", "psm_ids"]:
    if not os.path.exists(os.path.join(OUTPUT_DIR, sub_dir)):
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))

LEAD_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/leads/test_inclusion_q_value_leads_ml_posterior_error_prob.csv",
)
ALL_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/all/test_inclusion_q_value_all_ml_posterior_error_prob.csv",
)
PEPTIDE_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/peptides/test_inclusion_q_value_leads_peptides_ml_posterior_error_prob.csv",
)
PSM_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/psms/test_inclusion_q_value_leads_psms_ml_posterior_error_prob.csv",
)
PSM_ID_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/psm_ids/test_inclusion_q_value_leads_psm_ids_ml_posterior_error_prob.csv",
)

IDENTIFIER_INDEX = 0
SCORE_INDEX = 1
Q_VALUE_INDEX = 2
GROUP_ID_INDEX = 5
PEPTIDES_INDEX = 6

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(
    "protein_inference.tests.test_016_test_faulty_database_pipeline"
)


class TestFaultyDatabasePipeline(TestCase):
    def test_faulty_database_pipeline(self):

        protein_inference_parameters = ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        # Perform digest on each database
        digest_normal = in_silico_digest.InSilicoDigest(
            database_path=TEST_DATABASE_ORIGINAL,
            parameter_file_object=protein_inference_parameters,
            id_splitting=True,
        )
        digest_normal.digest_fasta_database()

        digest_missing_pep = in_silico_digest.InSilicoDigest(
            database_path=TEST_DATABASE_MISSING_PEP,
            parameter_file_object=protein_inference_parameters,
            id_splitting=True,
        )
        digest_missing_pep.digest_fasta_database()

        digest_missing_prot = in_silico_digest.InSilicoDigest(
            database_path=TEST_DATABASE_MISSING_PROT,
            parameter_file_object=protein_inference_parameters,
            id_splitting=True,
        )
        digest_missing_prot.digest_fasta_database()

        # Read PSM data with each of the digests with varying degrees of missingness
        pep_and_prot_data = protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest_class=digest_normal,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        pep_and_prot_data_miss_pep = protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest_class=digest_missing_pep,
        )
        pep_and_prot_data_miss_pep.read_psms()

        self.assertEqual(len(pep_and_prot_data_miss_pep.psms), 27)

        pep_and_prot_data_miss_prot = protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest_class=digest_missing_prot,
        )
        pep_and_prot_data_miss_prot.read_psms()

        self.assertEqual(len(pep_and_prot_data_miss_prot.psms), 27)

        # Load the digests and reader objects into datastore objects
        data_normal = protein_inference.datastore.DataStore(
            pep_and_prot_data, digest_class=digest_normal
        )

        data_miss_pep = protein_inference.datastore.DataStore(
            pep_and_prot_data_miss_pep, digest_class=digest_missing_pep
        )

        data_miss_prot = protein_inference.datastore.DataStore(
            pep_and_prot_data_miss_prot, digest_class=digest_missing_prot
        )

        for i in range(len(data_normal.main_data_form)):
            self.assertEqual(
                data_normal.main_data_form[0].identifier,
                data_miss_pep.main_data_form[0].identifier,
            )
            self.assertEqual(
                data_normal.main_data_form[0].non_flanking_peptide,
                data_miss_pep.main_data_form[0].non_flanking_peptide,
            )
            self.assertEqual(
                data_normal.main_data_form[0].pepvalue,
                data_miss_pep.main_data_form[0].pepvalue,
            )
            self.assertEqual(
                data_normal.main_data_form[0].percscore,
                data_miss_pep.main_data_form[0].percscore,
            )
            self.assertListEqual(
                data_normal.main_data_form[0].possible_proteins,
                data_miss_pep.main_data_form[0].possible_proteins,
            )
            self.assertEqual(
                data_normal.main_data_form[0].psm_id,
                data_miss_pep.main_data_form[0].psm_id,
            )
            self.assertEqual(
                data_normal.main_data_form[0].qvalue,
                data_miss_pep.main_data_form[0].qvalue,
            )
            self.assertEqual(
                data_normal.main_data_form[0].stripped_peptide,
                data_miss_pep.main_data_form[0].stripped_peptide,
            )

        for i in range(len(data_normal.main_data_form)):
            self.assertEqual(
                data_normal.main_data_form[0].identifier,
                data_miss_prot.main_data_form[0].identifier,
            )
            self.assertEqual(
                data_normal.main_data_form[0].non_flanking_peptide,
                data_miss_prot.main_data_form[0].non_flanking_peptide,
            )
            self.assertEqual(
                data_normal.main_data_form[0].pepvalue,
                data_miss_prot.main_data_form[0].pepvalue,
            )
            self.assertEqual(
                data_normal.main_data_form[0].percscore,
                data_miss_prot.main_data_form[0].percscore,
            )
            self.assertListEqual(
                data_normal.main_data_form[0].possible_proteins,
                data_miss_prot.main_data_form[0].possible_proteins,
            )
            self.assertEqual(
                data_normal.main_data_form[0].psm_id,
                data_miss_prot.main_data_form[0].psm_id,
            )
            self.assertEqual(
                data_normal.main_data_form[0].qvalue,
                data_miss_prot.main_data_form[0].qvalue,
            )
            self.assertEqual(
                data_normal.main_data_form[0].stripped_peptide,
                data_miss_prot.main_data_form[0].stripped_peptide,
            )

        stripped_peptides = [x.stripped_peptide for x in data_normal.main_data_form]
        possible_proteins = [x.possible_proteins for x in data_normal.main_data_form]

        proteins = set([y for x in possible_proteins for y in x])

        # Loop over all peptides from the normal data and make sure the peptide a key in peptide_to_protein dict
        for peps in stripped_peptides:
            self.assertIn(peps, digest_normal.peptide_to_protein_dictionary.keys())
            self.assertIn(peps, digest_missing_pep.peptide_to_protein_dictionary.keys())
            self.assertIn(
                peps, digest_missing_prot.peptide_to_protein_dictionary.keys()
            )
            self.assertSetEqual(
                digest_normal.peptide_to_protein_dictionary[peps],
                digest_missing_pep.peptide_to_protein_dictionary[peps],
            )
            self.assertSetEqual(
                digest_normal.peptide_to_protein_dictionary[peps],
                digest_missing_prot.peptide_to_protein_dictionary[peps],
            )

        for prots in proteins:
            self.assertIn(prots, digest_normal.protein_to_peptide_dictionary.keys())
            self.assertIn(
                prots, digest_missing_pep.protein_to_peptide_dictionary.keys()
            )
            self.assertIn(
                prots, digest_missing_prot.protein_to_peptide_dictionary.keys()
            )

        # Check to make sure datastore dictionaries are equal
        self.assertDictEqual(
            data_normal.peptide_to_protein_dictionary(),
            data_miss_pep.peptide_to_protein_dictionary(),
        )
        self.assertDictEqual(
            data_normal.peptide_to_protein_dictionary(),
            data_miss_prot.peptide_to_protein_dictionary(),
        )

        self.assertDictEqual(
            data_normal.protein_to_peptide_dictionary(),
            data_miss_pep.protein_to_peptide_dictionary(),
        )
        self.assertDictEqual(
            data_normal.protein_to_peptide_dictionary(),
            data_miss_prot.protein_to_peptide_dictionary(),
        )
