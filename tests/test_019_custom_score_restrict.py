from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference
from pyproteininference import in_silico_digest
from pyproteininference.parameters import ProteinInferenceParameter

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")

TARGET_FILE_MULTIPLICATIVE = resource_filename(
    "pyproteininference", "../tests/data/test_perc_data_target_multiplicative.txt"
)
DECOY_FILE_MULTIPLICATIVE = resource_filename(
    "pyproteininference", "../tests/data/test_perc_data_decoy_multiplicative.txt"
)
PARAMETER_FILE_MULTIPLICATIVE = resource_filename(
    "pyproteininference", "../tests/data/test_params_multiplicative_custom_score.yaml"
)

TARGET_FILE_ADDITIVE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target_additive.txt")
DECOY_FILE_ADDITIVE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy_additive.txt")
PARAMETER_FILE_ADDITIVE = resource_filename(
    "pyproteininference", "../tests/data/test_params_additive_custom_score.yaml"
)


class TestCustomRestrict(TestCase):
    def test_override_params_and_mult_restrict(self):

        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE_MULTIPLICATIVE)

        # Redefine some of the param options...
        protein_inference_parameters.restrict_pep = 0.9
        protein_inference_parameters.restrict_q = 0.1
        protein_inference_parameters.restrict_custom = 0.0009

        digest = in_silico_digest.PyteomicsDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        pep_and_prot_data = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE_MULTIPLICATIVE,
            decoy_file=DECOY_FILE_MULTIPLICATIVE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        data = pyproteininference.datastore.DataStore(pep_and_prot_data, digest=digest)

        data.restrict_psm_data()

        # Make sure we restrict by custom score
        self.assertEqual(len(data.main_data_restricted), 11)

        # Make sure other restrict params get set to None
        self.assertIsNone(protein_inference_parameters.restrict_pep)
        self.assertIsNone(protein_inference_parameters.restrict_q)

    def test_override_params_and_add_restrict(self):

        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE_ADDITIVE)

        # Redefine some of the param options...
        protein_inference_parameters.restrict_pep = 0.9
        protein_inference_parameters.restrict_q = 0.1
        protein_inference_parameters.restrict_custom = 6

        digest = in_silico_digest.PyteomicsDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        pep_and_prot_data = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE_ADDITIVE,
            decoy_file=DECOY_FILE_ADDITIVE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        data = pyproteininference.datastore.DataStore(pep_and_prot_data, digest=digest)

        data.restrict_psm_data()

        # Make sure custom score gets restricted
        self.assertEqual(len(data.main_data_restricted), 9)

        # Make sure parameters have been set to None
        self.assertIsNone(protein_inference_parameters.restrict_pep)
        self.assertIsNone(protein_inference_parameters.restrict_q)
