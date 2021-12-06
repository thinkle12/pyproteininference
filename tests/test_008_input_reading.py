from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference
from pyproteininference import in_silico_digest
from pyproteininference.parameters import ProteinInferenceParameter

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")
TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_peptide_centric.yaml")

TARGET_FILE_ADDITIVE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target_additive.txt")
DECOY_FILE_ADDITIVE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy_additive.txt")

TARGET_FILE_MULTIPLICATIVE = resource_filename(
    "pyproteininference", "../tests/data/test_perc_data_target_multiplicative.txt"
)
DECOY_FILE_MULTIPLICATIVE = resource_filename(
    "pyproteininference", "../tests/data/test_perc_data_decoy_multiplicative.txt"
)

PARAMETER_FILE_ADDITIVE = resource_filename(
    "pyproteininference", "../tests/data/test_params_additive_custom_score.yaml"
)
PARAMETER_FILE_MULTIPLICATIVE = resource_filename(
    "pyproteininference", "../tests/data/test_params_multiplicative_custom_score.yaml"
)


class TestReader(TestCase):
    def test_generic_and_percolator_readers(self):

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        digest = in_silico_digest.PyteomicsDigest(
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
            parameter_file_object=protein_inference_parameters,
            digest=digest,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        pep_and_prot_data_perc = pyproteininference.reader.PercolatorReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
        )
        pep_and_prot_data_perc.read_psms()

        self.assertEqual(len(pep_and_prot_data_perc.psms), 27)

        # Test to see if the values between the generic reader and the Percolator reader are identical
        for k in range(len(pep_and_prot_data.psms)):
            cur_generic = pep_and_prot_data.psms[k]
            cur_perc = pep_and_prot_data_perc.psms[k]
            self.assertEqual(cur_generic.identifier, cur_perc.identifier)
            self.assertEqual(cur_generic.pepvalue, cur_perc.pepvalue)
            self.assertEqual(cur_generic.percscore, cur_perc.percscore)
            self.assertEqual(cur_generic.possible_proteins, cur_perc.possible_proteins)
            self.assertEqual(cur_generic.psm_id, cur_perc.psm_id)
            self.assertEqual(cur_generic.qvalue, cur_perc.qvalue)

        # Test to see if we can read a list of files
        pep_and_prot_data = pyproteininference.reader.GenericReader(
            target_file=[TARGET_FILE, TARGET_FILE],
            decoy_file=[DECOY_FILE, DECOY_FILE],
            parameter_file_object=protein_inference_parameters,
            digest=digest,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 54)

        pep_and_prot_data_perc = pyproteininference.reader.PercolatorReader(
            target_file=[TARGET_FILE, TARGET_FILE],
            decoy_file=[DECOY_FILE, DECOY_FILE],
            parameter_file_object=protein_inference_parameters,
            digest=digest,
        )
        pep_and_prot_data_perc.read_psms()

        self.assertEqual(len(pep_and_prot_data_perc.psms), 54)

        # Test to see if the values between the generic reader and the Percolator reader are identical
        for k in range(len(pep_and_prot_data.psms)):
            cur_generic = pep_and_prot_data.psms[k]
            cur_perc = pep_and_prot_data_perc.psms[k]
            self.assertEqual(cur_generic.identifier, cur_perc.identifier)
            self.assertEqual(cur_generic.pepvalue, cur_perc.pepvalue)
            self.assertEqual(cur_generic.percscore, cur_perc.percscore)
            self.assertEqual(cur_generic.possible_proteins, cur_perc.possible_proteins)
            self.assertEqual(cur_generic.psm_id, cur_perc.psm_id)
            self.assertEqual(cur_generic.qvalue, cur_perc.qvalue)

    def test_generic_functionality(self):

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters_add = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE_ADDITIVE)

        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        digest = in_silico_digest.PyteomicsDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters_add.digest_type,
            missed_cleavages=protein_inference_parameters_add.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters_add.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters_add.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        pep_and_prot_data_add = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE_ADDITIVE,
            decoy_file=DECOY_FILE_ADDITIVE,
            parameter_file_object=protein_inference_parameters_add,
            digest=digest,
        )
        pep_and_prot_data_add.read_psms()

        self.assertEqual(len(pep_and_prot_data_add.psms), 27)

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters_mult = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE_MULTIPLICATIVE)

        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        digest = in_silico_digest.PyteomicsDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters_mult.digest_type,
            missed_cleavages=protein_inference_parameters_mult.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters_mult.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters_mult.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        pep_and_prot_data_mult = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE_MULTIPLICATIVE,
            decoy_file=DECOY_FILE_MULTIPLICATIVE,
            parameter_file_object=protein_inference_parameters_mult,
            digest=digest,
            append_alt_from_db=False,
        )
        pep_and_prot_data_mult.read_psms()

        self.assertEqual(len(pep_and_prot_data_mult.psms), 27)
