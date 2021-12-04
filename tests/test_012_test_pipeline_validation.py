import os
import shutil
import tempfile
from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_inclusion.yaml")
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('pyproteininference', '../tests/output/')

TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
COMBINED_FILE = resource_filename("pyproteininference", "../tests/data/combined_files/test_combined_data.txt")

temp_dir = tempfile.gettempdir()

TARGET_DIRECTORY = os.path.join(temp_dir, "target_directory")
DECOY_DIRECTORY = os.path.join(temp_dir, "decoy_directory")
if not os.path.exists(TARGET_DIRECTORY):
    os.makedirs(TARGET_DIRECTORY)
if not os.path.exists(DECOY_DIRECTORY):
    os.makedirs(DECOY_DIRECTORY)

shutil.copyfile(str(TARGET_FILE), os.path.join(TARGET_DIRECTORY, "target_file.txt"))
shutil.copyfile(DECOY_FILE, os.path.join(DECOY_DIRECTORY, "decoy_file.txt"))
COMBINED_DIRECTORY = resource_filename("pyproteininference", "../tests/data/combined_files")


class TestPipelineValidation(TestCase):
    def test_validation(self):

        # Test target and decoy file
        pipeline1 = pyproteininference.pipeline.ProteinInferencePipeline(
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

        pipeline1.execute()

        self.assertEqual(len(pipeline1.data.main_data_form), 27)

        # Test combined file
        pipeline2 = pyproteininference.pipeline.ProteinInferencePipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=None,
            decoy_files=None,
            combined_files=COMBINED_FILE,
            target_directory=None,
            decoy_directory=None,
            combined_directory=None,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
        )

        pipeline2.execute()

        self.assertEqual(len(pipeline2.data.main_data_form), 27)

        # Test target and decoy directories
        pipeline3 = pyproteininference.pipeline.ProteinInferencePipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=None,
            decoy_files=None,
            combined_files=None,
            target_directory=TARGET_DIRECTORY,
            decoy_directory=DECOY_DIRECTORY,
            combined_directory=None,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
        )

        self.assertTrue(os.path.exists(pipeline3.target_files[0]))
        self.assertTrue(os.path.exists(pipeline3.decoy_files[0]))

        pipeline3.execute()

        self.assertEqual(len(pipeline3.data.main_data_form), 27)

        # Test combined directory
        pipeline4 = pyproteininference.pipeline.ProteinInferencePipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=None,
            decoy_files=None,
            combined_files=None,
            target_directory=None,
            decoy_directory=None,
            combined_directory=COMBINED_DIRECTORY,
            output_directory=OUTPUT_DIR,
            id_splitting=True,
        )

        self.assertTrue(os.path.exists(pipeline4.combined_files[0]))

        pipeline4.execute()

        self.assertEqual(len(pipeline4.data.main_data_form), 27)

        # Test Proper error reporting
        with self.assertRaises(ValueError):
            pyproteininference.pipeline.ProteinInferencePipeline(
                parameter_file=PARAMETER_FILE,
                database_file=TEST_DATABASE,
                target_files=TARGET_FILE,
                decoy_files=DECOY_FILE,
                combined_files=None,
                target_directory=TARGET_DIRECTORY,
                decoy_directory=DECOY_DIRECTORY,
                combined_directory=None,
                output_directory=OUTPUT_DIR,
                id_splitting=True,
            )

        with self.assertRaises(ValueError):
            pyproteininference.pipeline.ProteinInferencePipeline(
                parameter_file=PARAMETER_FILE,
                database_file=TEST_DATABASE,
                target_files=None,
                decoy_files=None,
                combined_files=COMBINED_FILE,
                target_directory=None,
                decoy_directory=None,
                combined_directory=COMBINED_DIRECTORY,
                output_directory=OUTPUT_DIR,
                id_splitting=True,
            )

        with self.assertRaises(ValueError):
            pyproteininference.pipeline.ProteinInferencePipeline(
                parameter_file=PARAMETER_FILE,
                database_file=TEST_DATABASE,
                target_files=TARGET_FILE,
                decoy_files=DECOY_FILE,
                combined_files=COMBINED_FILE,
                target_directory=TARGET_DIRECTORY,
                decoy_directory=DECOY_DIRECTORY,
                combined_directory=COMBINED_DIRECTORY,
                output_directory=OUTPUT_DIR,
                id_splitting=True,
            )
