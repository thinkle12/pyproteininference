import os
import tempfile
from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")
TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_parsimony_pulp.yaml")
OUTPUT_DIR = tempfile.gettempdir()
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
        )

        hp.execute(skip_plot=True)

        self.assertEqual(hp.selected_method, "parsimony")

        self.assertIsInstance(hp.selected_datastore, pyproteininference.datastore.DataStore)

        result1 = hp.determine_optimal_inference_method(empirical_threshold=0.5)

        self.assertEqual(result1, "exclusion")

        result2 = hp.determine_optimal_inference_method(empirical_threshold=1.2)

        self.assertEqual(result2, "exclusion")

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
        )

        hp.execute(skip_plot=True)

        self.assertEqual(hp.selected_method, "parsimony")

        self.assertIsInstance(hp.selected_datastore, pyproteininference.datastore.DataStore)

        result1 = hp.determine_optimal_inference_method(empirical_threshold=0.5)

        self.assertEqual(result1, "exclusion")

        result2 = hp.determine_optimal_inference_method(empirical_threshold=1.2)

        self.assertEqual(result2, "exclusion")
