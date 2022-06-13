import os
import tempfile
from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")
TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_heuristic.yaml")
OUTPUT_DIR = tempfile.gettempdir()
PDF_FILENAME = os.path.join(OUTPUT_DIR, "test_pyproteininference.pdf")
# OUTPUT_DIR = resource_filename('pyproteininference', '../tests/output/')
for sub_dir in ["leads", "all", "peptides", "psms", "psm_ids"]:
    if not os.path.exists(os.path.join(OUTPUT_DIR, sub_dir)):
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))


class TestHeuristicWorkflow(TestCase):
    def test_workflow_heuristic_with_params(self):

        hp = pyproteininference.heuristic.HeuristicPipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=TARGET_FILE,
            decoy_files=DECOY_FILE,
            combined_files=None,
            target_directory=None,
            decoy_directory=None,
            combined_directory=None,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
            pdf_filename=PDF_FILENAME,
        )

        hp.execute()

        self.assertListEqual(hp.selected_methods, ["inclusion"])

        self.assertIsInstance(hp.selected_datastores["inclusion"], pyproteininference.datastore.DataStore)

        result1 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.05,
            upper_empirical_threshold=2,
            lower_empirical_threshold=1,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result1, ["parsimony"])

        result2 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.1,
            upper_empirical_threshold=0.5,
            lower_empirical_threshold=0.25,
            pdf_filename=PDF_FILENAME,
        )
        self.assertListEqual(result2, ["inclusion"])

    def test_workflow_heuristic_without_params(self):

        hp = pyproteininference.heuristic.HeuristicPipeline(
            database_file=TEST_DATABASE,
            target_files=TARGET_FILE,
            decoy_files=DECOY_FILE,
            combined_files=None,
            target_directory=None,
            decoy_directory=None,
            combined_directory=None,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
            pdf_filename=PDF_FILENAME,
        )

        hp.execute()

        self.assertListEqual(hp.selected_methods, ["inclusion"])

        self.assertIsInstance(hp.selected_datastores["inclusion"], pyproteininference.datastore.DataStore)

        result1 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.05,
            upper_empirical_threshold=2,
            lower_empirical_threshold=1,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result1, ["parsimony"])

        result2 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.1,
            upper_empirical_threshold=0.5,
            lower_empirical_threshold=0.25,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result2, ["inclusion"])

    def test_workflow_heuristic_optimal_export(self):

        hp = pyproteininference.heuristic.HeuristicPipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=TARGET_FILE,
            decoy_files=DECOY_FILE,
            combined_files=None,
            target_directory=None,
            decoy_directory=None,
            combined_directory=None,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
            pdf_filename=PDF_FILENAME,
            output_type="optimal",
        )

        hp.execute()

        self.assertListEqual(hp.selected_methods, ["inclusion"])

        self.assertIsInstance(hp.selected_datastores["inclusion"], pyproteininference.datastore.DataStore)

        result1 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.05,
            upper_empirical_threshold=2,
            lower_empirical_threshold=1,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result1, ["parsimony"])

        result3 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.1,
            upper_empirical_threshold=0.5,
            lower_empirical_threshold=0.25,
            pdf_filename=PDF_FILENAME,
        )
        self.assertListEqual(result3, ["inclusion"])

    def test_workflow_heuristic_with_different_thresholds(self):

        hp = pyproteininference.heuristic.HeuristicPipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=TARGET_FILE,
            decoy_files=DECOY_FILE,
            combined_files=None,
            target_directory=None,
            decoy_directory=None,
            combined_directory=None,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
            pdf_filename=PDF_FILENAME,
        )

        hp.execute()

        # Both inclusion/exclusion passing high upper threshold but neither peptide-centric or
        # parsimony passing lower threshold
        result1 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.05,
            upper_empirical_threshold=10,
            lower_empirical_threshold=0.1,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result1, ["inclusion", "exclusion"])

        # Neither inclusion/exclusion passing low upper threshold but peptide-centric and
        # parsimony both passing high lower threshold
        result2 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.1,
            upper_empirical_threshold=2,
            lower_empirical_threshold=2,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result2, ["parsimony", "peptide_centric"])

        # No methods passing thresholds
        result3 = hp.determine_optimal_inference_method(
            false_discovery_rate_threshold=0.1,
            upper_empirical_threshold=0.00001,
            lower_empirical_threshold=0.000001,
            pdf_filename=PDF_FILENAME,
        )

        self.assertListEqual(result3, ["inclusion"])
