import os
import tempfile
from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference
from pyproteininference import in_silico_digest
from pyproteininference.parameters import ProteinInferenceParameter

TEST_DATABASE = resource_filename("pyproteininference", "../tests/data/test_database.fasta")
TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_inclusion.yaml")
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('pyproteininference', '../tests/output/')
for sub_dir in ["leads", "all", "peptides", "psms", "psm_ids"]:
    if not os.path.exists(os.path.join(OUTPUT_DIR, sub_dir)):
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))


class TestAltProteinRead(TestCase):
    def test_alt_protein_reader(self):

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
            append_alt_from_db=True,
        )
        pep_and_prot_data.read_psms()

        possible_proteins = [x.possible_proteins for x in pep_and_prot_data.psms]

        true_possible_proteins = [
            ['RAF1_HUMAN|P04049'],
            ['RAF1_HUMAN|P04049'],
            ['ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056', 'RAF1_HUMAN|P04049'],
            ['ARAF_HUMAN|P10398', 'RAF1_HUMAN|P04049'],
            ['BRAF_HUMAN|P15056', 'RAF1_HUMAN|P04049'],
            ['ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056'],
            ['ARAF_HUMAN|P10398'],
            ['RAF1_HUMAN|P04049'],
            ['TCAF1_HUMAN|Q9Y4C2'],
            ['TCAF1_HUMAN|Q9Y4C2', 'TCAF1_HUMAN|Q9Y4C2-2'],
            ['HNRPU_HUMAN|Q00839', 'B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7'],
            ['HNRPU_HUMAN|Q00839', 'B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['##TCAF1_HUMAN|##Q9Y4C2'],
            [
                '##TCAF2_HUMAN|##A6NFQ2',
                '##TCAF2_HUMAN|##A6NFQ2-2',
                '##TCAF2_HUMAN|##A6NFQ2-3',
            ],
            ['RAF1_HUMAN|P04049'],
        ]

        for i in range(len(possible_proteins)):
            self.assertSetEqual(set(possible_proteins[i]), set(true_possible_proteins[i]))

        pep_and_prot_data_no_append = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
            append_alt_from_db=False,
        )
        pep_and_prot_data_no_append.read_psms()

        possible_proteins_no_append = [x.possible_proteins for x in pep_and_prot_data_no_append.psms]

        true_possible_proteins_no_append = [
            ['RAF1_HUMAN|P04049'],
            ['RAF1_HUMAN|P04049'],
            ['ARAF_HUMAN|P10398', 'RAF1_HUMAN|P04049', 'BRAF_HUMAN|P15056'],
            ['ARAF_HUMAN|P10398', 'RAF1_HUMAN|P04049'],
            ['RAF1_HUMAN|P04049', 'BRAF_HUMAN|P15056'],
            ['ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056'],
            ['ARAF_HUMAN|P10398'],
            ['RAF1_HUMAN|P04049'],
            ['TCAF1_HUMAN|Q9Y4C2'],
            ['TCAF1_HUMAN|Q9Y4C2'],
            ['HNRPU_HUMAN|Q00839', 'B3KX72_HUMAN|B3KX72'],
            ['HNRPU_HUMAN|Q00839', 'B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['RPOC_SHIF8|Q0SY12'],
            ['##TCAF1_HUMAN|##Q9Y4C2'],
            ['##TCAF2_HUMAN|##A6NFQ2'],
            ['RAF1_HUMAN|P04049'],
        ]

        for i in range(len(possible_proteins_no_append)):
            self.assertSetEqual(
                set(possible_proteins_no_append[i]),
                set(true_possible_proteins_no_append[i]),
            )
