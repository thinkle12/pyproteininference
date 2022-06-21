from unittest import TestCase

from pkg_resources import resource_filename

import pyproteininference
from pyproteininference import in_silico_digest
from pyproteininference.parameters import ProteinInferenceParameter

TARGET_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_target.txt")
DECOY_FILE = resource_filename("pyproteininference", "../tests/data/test_perc_data_decoy.txt")
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_inclusion.yaml")
COMBINED_FILE = resource_filename("pyproteininference", "../tests/data/test_data_many_alternative_proteins.txt")


class TestAltProteinFomInput(TestCase):
    def test_read_alt_proteins(self):

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        digest = in_silico_digest.PyteomicsDigest(
            database_path=None,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )

        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        reader = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
            append_alt_from_db=False,
        )
        reader.read_psms()

        expected_result = [
            ['RAF1_HUMAN|P04049'],
            ['RAF1_HUMAN|P04049'],
            ['RAF1_HUMAN|P04049', 'ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056'],
            ['RAF1_HUMAN|P04049', 'ARAF_HUMAN|P10398'],
            ['RAF1_HUMAN|P04049', 'BRAF_HUMAN|P15056'],
            ['ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056'],
            ['ARAF_HUMAN|P10398'],
            ['RAF1_HUMAN|P04049'],
            ['TCAF1_HUMAN|Q9Y4C2'],
            ['TCAF1_HUMAN|Q9Y4C2'],
            ['B3KX72_HUMAN|B3KX72', 'HNRPU_HUMAN|Q00839'],
            ['B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7', 'HNRPU_HUMAN|Q00839'],
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

        for i in range(len(reader.psms)):
            self.assertSetEqual(set(reader.psms[i].possible_proteins), set(expected_result[i]))

        self.assertEqual(len(reader.psms), 27)

    def test_read_alt_proteins_over_maximum(self):

        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        digest = in_silico_digest.PyteomicsDigest(
            database_path=None,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )

        # Try with generic reader

        reader = pyproteininference.reader.GenericReader(
            combined_files=COMBINED_FILE,
            parameter_file_object=protein_inference_parameters,
            digest=digest,
            append_alt_from_db=True,
        )
        reader.read_psms()

        self.assertEqual(reader.MAX_ALLOWED_ALTERNATIVE_PROTEINS, 50)

        protein_list1 = [
            'Protein1',
            'Protein10',
            'Protein11',
            'Protein12',
            'Protein13',
            'Protein14',
            'Protein15',
            'Protein16',
            'Protein17',
            'Protein18',
            'Protein19',
            'Protein2',
            'Protein20',
            'Protein21',
            'Protein22',
            'Protein23',
            'Protein24',
            'Protein25',
            'Protein26',
            'Protein27',
            'Protein28',
            'Protein29',
            'Protein3',
            'Protein30',
            'Protein31',
            'Protein32',
            'Protein33',
            'Protein34',
            'Protein35',
            'Protein36',
            'Protein37',
            'Protein38',
            'Protein39',
            'Protein4',
            'Protein40',
            'Protein41',
            'Protein42',
            'Protein43',
            'Protein44',
            'Protein45',
            'Protein46',
            'Protein47',
            'Protein48',
            'Protein49',
            'Protein5',
            'Protein50',
            'Protein51',
            'Protein52',
            'Protein53',
            'Protein54',
        ]

        protein_list2 = [
            '##Protein1',
            '##Protein10',
            '##Protein11',
            '##Protein12',
            '##Protein13',
            '##Protein14',
            '##Protein15',
            '##Protein16',
            '##Protein17',
            '##Protein18',
            '##Protein19',
            '##Protein2',
            '##Protein20',
            '##Protein21',
            '##Protein22',
            '##Protein23',
            '##Protein24',
            '##Protein25',
            '##Protein26',
            '##Protein27',
            '##Protein28',
            '##Protein29',
            '##Protein3',
            '##Protein30',
            '##Protein31',
            '##Protein32',
            '##Protein33',
            '##Protein34',
            '##Protein35',
            '##Protein36',
            '##Protein37',
            '##Protein38',
            '##Protein39',
            '##Protein4',
            '##Protein40',
            '##Protein41',
            '##Protein42',
            '##Protein43',
            '##Protein44',
            '##Protein45',
            '##Protein46',
            '##Protein47',
            '##Protein48',
            '##Protein49',
            '##Protein5',
            '##Protein50',
            '##Protein51',
            '##Protein52',
            '##Protein53',
            '##Protein54',
        ]

        self.assertListEqual(
            reader.psms[0].possible_proteins,
            protein_list1,
        )
        self.assertListEqual(
            reader.psms[1].possible_proteins,
            protein_list2,
        )

        percdigest = in_silico_digest.PyteomicsDigest(
            database_path=None,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )

        # Try with percolator reader
        percreader = pyproteininference.reader.PercolatorReader(
            combined_files=COMBINED_FILE,
            parameter_file_object=protein_inference_parameters,
            digest=percdigest,
            append_alt_from_db=True,
        )

        percreader.read_psms()

        self.assertEqual(percreader.MAX_ALLOWED_ALTERNATIVE_PROTEINS, 50)

        self.assertListEqual(
            percreader.psms[0].possible_proteins,
            protein_list1,
        )
        self.assertListEqual(
            percreader.psms[1].possible_proteins,
            protein_list2,
        )
