import csv
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


class TestExportTypes(TestCase):
    def test_export_types(self):

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters = pyproteininference.parameters.ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        self.assertEqual(protein_inference_parameters.digest_type, "trypsin")
        self.assertEqual(protein_inference_parameters.export, "q_value")
        self.assertEqual(protein_inference_parameters.fdr, 0.01)
        self.assertEqual(protein_inference_parameters.missed_cleavages, 3)
        self.assertEqual(protein_inference_parameters.picker, True)
        self.assertEqual(protein_inference_parameters.restrict_pep, 0.9)
        self.assertEqual(protein_inference_parameters.restrict_peptide_length, 7)
        self.assertEqual(protein_inference_parameters.restrict_q, 0.9)
        self.assertEqual(protein_inference_parameters.protein_score, "multiplicative_log")
        self.assertEqual(protein_inference_parameters.psm_score, "posterior_error_prob")
        self.assertEqual(protein_inference_parameters.psm_score_type, "multiplicative")
        self.assertEqual(protein_inference_parameters.decoy_symbol, "##")
        self.assertEqual(protein_inference_parameters.isoform_symbol, "-")
        self.assertEqual(protein_inference_parameters.reviewed_identifier_symbol, "sp|")
        self.assertEqual(protein_inference_parameters.inference_type, "parsimony")
        self.assertEqual(protein_inference_parameters.tag, "test_parsimony")
        self.assertEqual(protein_inference_parameters.grouping_type, "shared_peptides")
        self.assertEqual(protein_inference_parameters.max_identifiers_peptide_centric, 5)
        self.assertEqual(protein_inference_parameters.lp_solver, "pulp")

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
        score.score_psms(score_method=protein_inference_parameters.protein_score)

        # STEP 9: Run protein picker on the data
        # STEP 9: Run protein picker on the data
        # STEP 9: Run protein picker on the data
        if protein_inference_parameters.picker:
            data.protein_picker()
        else:
            pass

        # STEP 10: Apply Inference
        # STEP 10: Apply Inference
        # STEP 10: Apply Inference
        inference_type = protein_inference_parameters.inference_type

        if inference_type == pyproteininference.inference.Inference.PARSIMONY:
            group = pyproteininference.inference.Parsimony(data=data, digest=digest)
            group.infer_proteins()

        if inference_type == pyproteininference.inference.Inference.INCLUSION:
            group = pyproteininference.inference.Inclusion(data=data, digest=digest)
            group.infer_proteins()

        if inference_type == pyproteininference.inference.Inference.EXCLUSION:
            group = pyproteininference.inference.Exclusion(data=data, digest=digest)
            group.infer_proteins()

        # STEP 11: Run FDR and Q value Calculations
        # STEP 11: Run FDR and Q value Calculations
        # STEP 11: Run FDR and Q value Calculations
        data.calculate_q_values()

        export_type = "peptides"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have peptides space separated
        self.assertListEqual(
            output[1],
            [
                "RPOC_SHIF8|Q0SY12",
                "82.89306334778564",
                "0.0",
                "12",
                "Reviewed",
                "1",
                "CGVEVTQTK EGLNVLQY#FISTHGAR FATSDLNDLYR IALASPDMIR IPQESGGTK LIPAGTGYAYHQDR "
                "MGAEAIQALLK NTLLHEQWCDLLEENSVDAVK RVDYSGR VADLFEAR VIDIWAAANDR VTAEDVLKPGTADILVPR",
            ],
        )

        export_type = "psms"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have psms space separated
        self.assertListEqual(
            output[1],
            [
                "RPOC_SHIF8|Q0SY12",
                "82.89306334778564",
                "0.0",
                "12",
                "Reviewed",
                "1",
                "CGVEVTQTK EGLNVLQY#FISTHGAR FATSDLNDLYR IALASPDMIR IPQESGGTK LIPAGTGYAYHQDR "
                "MGAEAIQALLK NTLLHEQWCDLLEENSVDAVK RVDYSGR VADLFEAR VIDIWAAANDR VTAEDVLKPGTADILVPR",
            ],
        )

        export_type = "psm_ids"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have psm_ids space separated
        self.assertListEqual(
            output[1],
            [
                "RPOC_SHIF8|Q0SY12",
                "82.89306334778564",
                "0.0",
                "12",
                "Reviewed",
                "1",
                "13 14 15 16 17 18 19 20 21 22 23 24",
            ],
        )

        export_type = "q_value"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have peptides comma separated
        self.assertListEqual(
            output[1],
            [
                "RPOC_SHIF8|Q0SY12",
                "82.89306334778564",
                "0.0",
                "12",
                "Reviewed",
                "1",
                "CGVEVTQTK",
                "EGLNVLQY#FISTHGAR",
                "FATSDLNDLYR",
                "IALASPDMIR",
                "IPQESGGTK",
                "LIPAGTGYAYHQDR",
                "MGAEAIQALLK",
                "NTLLHEQWCDLLEENSVDAVK",
                "RVDYSGR",
                "VADLFEAR",
                "VIDIWAAANDR",
                "VTAEDVLKPGTADILVPR",
            ],
        )

        export_type = "q_value_all"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have peptides comma separated
        self.assertListEqual(
            output[1],
            [
                "RPOC_SHIF8|Q0SY12",
                "82.89306334778564",
                "0.0",
                "12",
                "Reviewed",
                "1",
                "CGVEVTQTK",
                "EGLNVLQY#FISTHGAR",
                "FATSDLNDLYR",
                "IALASPDMIR",
                "IPQESGGTK",
                "LIPAGTGYAYHQDR",
                "MGAEAIQALLK",
                "NTLLHEQWCDLLEENSVDAVK",
                "RVDYSGR",
                "VADLFEAR",
                "VIDIWAAANDR",
                "VTAEDVLKPGTADILVPR",
            ],
        )

        export_type = "q_value_comma_sep"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have proteins comma separated
        self.assertListEqual(
            output[2],
            [
                "RAF1_HUMAN|P04049",
                "70.7434325345954",
                "0.0",
                "6",
                "Reviewed",
                "2",
                "ARAF_HUMAN|P10398",
                "BRAF_HUMAN|P15056",
            ],
        )

        export_type = "leads"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have peptides comma separated with no q values
        self.assertListEqual(
            output[2],
            [
                "RAF1_HUMAN|P04049",
                "70.7434325345954",
                "6",
                "Reviewed",
                "{2, 3}",
                "CQTCGYKFHEHCSTK",
                "FQMFQLIDIAR",
                "QTAQGMDYLHAK",
                "SASEPSLHR",
                "VFLPNKQR",
                "WHGDVAVKILK",
            ],
        )

        export_type = "all"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have peptides comma separated with no q values
        self.assertListEqual(
            output[-1],
            [
                "Q96BA7_HUMAN|Q96BA7",
                "6.907755278982137",
                "1",
                "Unreviewed",
                "{5}",
                "LQAALDDEEAGGRPAMEPGNGSLDLGGDSAGR",
            ],
        )

        export_type = "comma_sep"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should have proteins comma separated with no q values
        self.assertListEqual(
            output[-1],
            [
                "HNRPU_HUMAN|Q00839",
                "15.316094065486292",
                "2",
                "Reviewed",
                "{5}",
                "B3KX72_HUMAN|B3KX72",
                "Q96BA7_HUMAN|Q96BA7",
            ],
        )

        export_type = "long"
        export = pyproteininference.export.Export(data=data)
        export.export_to_csv(directory=OUTPUT_DIR, export_type=export_type)

        output = []
        with open(export.filepath, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                output.append(row)

        # Output should be one line for every peptide
        self.assertListEqual(
            output[1],
            ['RPOC_SHIF8|Q0SY12', '82.89306334778564', '0.0', '12', 'Reviewed', '1', 'CGVEVTQTK'],
        )

        self.assertListEqual(
            output[2],
            ['RPOC_SHIF8|Q0SY12', '82.89306334778564', '0.0', '12', 'Reviewed', '1', 'EGLNVLQY#FISTHGAR'],
        )

        self.assertListEqual(
            output[-1],
            [
                '##TCAF2_HUMAN|##A6NFQ2',
                '2.4079456086518722',
                '0.3333333333333333',
                '1',
                'Reviewed',
                '6',
                'MEPTPVPFCGAK',
            ],
        )
