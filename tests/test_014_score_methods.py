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

temp_dir = tempfile.gettempdir()


class TestScoreMethods(TestCase):
    def test_scoring_methods(self):
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters = pyproteininference.parameters.ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        digest = pyproteininference.in_silico_digest.PyteomicsDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        pep_and_prot_data = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            combined_files=None,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
            append_alt_from_db=False,
        )
        pep_and_prot_data.read_psms()

        # STEP 4: Initiate the datastore class #
        # STEP 4: Initiate the datastore class #
        # STEP 4: Initiate the datastore class #
        data = pyproteininference.datastore.DataStore(pep_and_prot_data, digest=digest)

        # Step 5: Restrict the PSM data
        # Step 5: Restrict the PSM data
        # Step 5: Restrict the PSM data
        data.restrict_psm_data()

        # Step 6: Generate protein scoring input
        # Step 6: Generate protein scoring input
        # Step 6: Generate protein scoring input
        data.create_scoring_input()

        # Step 7: Remove non unique peptides if running exclusion
        # Step 7: Remove non unique peptides if running exclusion
        # Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters.inference_type == pyproteininference.inference.Inference.EXCLUSION:
            # This gets ran if we run exclusion...
            data.exclude_non_distinguishing_peptides()

        # STEP 8: Score our PSMs given a score method
        # STEP 8: Score our PSMs given a score method
        # STEP 8: Score our PSMs given a score method
        score = pyproteininference.scoring.Score(data=data)
        score.score_psms(score_method="best_peptide_per_protein")
        self.assertEqual(data.scored_proteins[0].score, 3.5e-06)

        score.score_psms(score_method="iterative_downweighted_log")
        self.assertEqual(data.scored_proteins[0].score, 64.16418132258529)

        score.score_psms(score_method="multiplicative_log")
        self.assertEqual(data.scored_proteins[0].score, 82.89306334778564)

        score.score_psms(score_method="downweighted_multiplicative_log")
        self.assertEqual(data.scored_proteins[0].score, 38.42870136439401)

        score.score_psms(score_method="downweighted_version2")
        self.assertEqual(data.scored_proteins[0].score, 29.723342219506584)

        score.score_psms(score_method="top_two_combined")
        self.assertEqual(data.scored_proteins[0].score, 25.523842819277114)

        score.score_psms(score_method="geometric_mean")
        self.assertEqual(data.scored_proteins[0].score, 11.790572089099232)
