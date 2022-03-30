from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference
from pyproteininference import in_silico_digest
from pyproteininference.parameters import ProteinInferenceParameter

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")
TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_exclusion.yaml")


class TestMappingRecovery(TestCase):
    def test_mapping_recovery(self):

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
            append_alt_from_db=False,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        # STEP 4: Initiate the datastore class #
        # STEP 4: Initiate the datastore class #
        # STEP 4: Initiate the datastore class #
        data = pyproteininference.datastore.DataStore(pep_and_prot_data, digest=digest)

        # Step 5: Restrict the PSM data
        # Step 5: Restrict the PSM data
        # Step 5: Restrict the PSM data
        data.restrict_psm_data()

        self.assertEqual(len(data.main_data_restricted), 26)

        del data.digest.protein_to_peptide_dictionary["Q96BA7_HUMAN|Q96BA7"]
        del data.digest.protein_to_peptide_dictionary["B3KX72_HUMAN|B3KX72"]

        for pep in list(data.digest.peptide_to_protein_dictionary.keys()):
            if "Q96BA7_HUMAN|Q96BA7" in data.digest.peptide_to_protein_dictionary[pep]:
                data.digest.peptide_to_protein_dictionary[pep].remove("Q96BA7_HUMAN|Q96BA7")
            if "B3KX72_HUMAN|B3KX72" in data.digest.peptide_to_protein_dictionary[pep]:
                data.digest.peptide_to_protein_dictionary[pep].remove("B3KX72_HUMAN|B3KX72")

        # Try to run exclusion with proteins missing from the database digest but existing in the main psm data
        with self.assertRaises(KeyError):
            data.create_scoring_input()
            data.exclude_non_distinguishing_peptides()

        # Recover the proteins and check the digests
        data.recover_mapping()

        self.assertIn("Q96BA7_HUMAN|Q96BA7", list(data.digest.protein_to_peptide_dictionary.keys()))
        self.assertIn("B3KX72_HUMAN|B3KX72", list(data.digest.protein_to_peptide_dictionary.keys()))
        self.assertIn(
            "Q96BA7_HUMAN|Q96BA7", data.digest.peptide_to_protein_dictionary['LQAALDDEEAGGRPAMEPGNGSLDLGGDSAGR']
        )
        self.assertIn(
            "B3KX72_HUMAN|B3KX72", data.digest.peptide_to_protein_dictionary['LQAALDDEEAGGRPAMEPGNGSLDLGGDSAGR']
        )
        self.assertIn("B3KX72_HUMAN|B3KX72", data.digest.peptide_to_protein_dictionary['AEGGGGGGRPGAPAAGDGK'])

        data.create_scoring_input()
        data.exclude_non_distinguishing_peptides()
