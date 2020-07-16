#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""
import csv
import tempfile
from unittest import TestCase
from pkg_resources import resource_filename

import protein_inference
from protein_inference import in_silico_digest
from protein_inference.parameters import ProteinInferenceParameter
import os
import logging

TEST_DATABASE = resource_filename(
    "protein_inference", "../tests/data/test_database.fasta"
)
TARGET_FILE = resource_filename(
    "protein_inference", "../tests/data/test_perc_data_target.txt"
)
DECOY_FILE = resource_filename(
    "protein_inference", "../tests/data/test_perc_data_decoy.txt"
)
PARAMETER_FILE = resource_filename(
    "protein_inference", "../tests/data/test_params_exclusion.yaml"
)
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')
for sub_dir in ["leads", "all", "peptides", "psms", "psm_ids"]:
    if not os.path.exists(os.path.join(OUTPUT_DIR, sub_dir)):
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))

LEAD_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/leads/test_exclusion_q_value_leads_ml_posterior_error_prob.csv",
)
ALL_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/all/test_exclusion_q_value_all_ml_posterior_error_prob.csv",
)
PEPTIDE_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/peptides/test_exclusion_q_value_leads_peptides_ml_posterior_error_prob.csv",
)
PSM_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/psms/test_exclusion_q_value_leads_psms_ml_posterior_error_prob.csv",
)
PSM_ID_OUTPUT_FILE = resource_filename(
    "protein_inference",
    "../tests/output/psm_ids/test_exclusion_q_value_leads_psm_ids_ml_posterior_error_prob.csv",
)

IDENTIFIER_INDEX = 0
SCORE_INDEX = 1
Q_VALUE_INDEX = 2
GROUP_ID_INDEX = 5
PEPTIDES_INDEX = 6

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("protein_inference.tests.test_003_exclusion_pipeline.py")


class TestLoadExclusionWorkflow(TestCase):
    def test_workflow_exclusion(self):

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        # TODO check to make sure proper parameters are being loaded
        self.assertEqual(protein_inference_parameters.digest_type, "trypsin")
        self.assertEqual(protein_inference_parameters.export, "q_value")
        self.assertEqual(protein_inference_parameters.fdr, 0.01)
        self.assertEqual(protein_inference_parameters.glpk_path, None)
        self.assertEqual(protein_inference_parameters.missed_cleavages, 3)
        self.assertEqual(protein_inference_parameters.picker, True)
        self.assertEqual(protein_inference_parameters.restrict_pep, 0.9)
        self.assertEqual(protein_inference_parameters.restrict_peptide_length, 7)
        self.assertEqual(protein_inference_parameters.restrict_q, 0.9)
        self.assertEqual(
            protein_inference_parameters.score_method, "multiplicative_log"
        )
        self.assertEqual(protein_inference_parameters.psm_score, "posterior_error_prob")
        self.assertEqual(protein_inference_parameters.score_type, "multiplicative")
        self.assertEqual(protein_inference_parameters.decoy_symbol, "##")
        self.assertEqual(protein_inference_parameters.isoform_symbol, "-")
        self.assertEqual(protein_inference_parameters.reviewed_identifier_symbol, "sp|")
        self.assertEqual(protein_inference_parameters.inference_type, "exclusion")
        self.assertEqual(protein_inference_parameters.tag, "test_exclusion")
        self.assertEqual(protein_inference_parameters.grouping_type, None)
        self.assertEqual(
            protein_inference_parameters.max_identifiers_peptide_centric, 5
        )
        self.assertEqual(protein_inference_parameters.lp_solver, "None")

        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = in_silico_digest.InSilicoDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters.digest_type,
            missed_cleavages=protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters.restrict_peptide_length,
            id_splitting=True,
        )
        digest.digest_fasta_database()

        # TODO dont test digest... test this in a separate test

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data = protein_inference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters,
            digest_class=digest,
        )
        pep_and_prot_data.read_psms()

        self.assertEqual(len(pep_and_prot_data.psms), 27)

        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data = protein_inference.datastore.DataStore(
            pep_and_prot_data, digest_class=digest
        )

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data.restrict_psm_data()

        self.assertEqual(len(data.main_data_restricted), 26)

        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        data.create_scoring_input()

        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters.inference_type == protein_inference.inference.Inference.EXCLUSION:
            # This gets ran if we run exclusion...
            data.exclude_non_distinguishing_peptides()

        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        score = protein_inference.scoring.Score(data_class=data)
        score.score_psms(score_method=protein_inference_parameters.score_method)

        ### STEP 9: Run protein picker on the data
        ### STEP 9: Run protein picker on the data
        ### STEP 9: Run protein picker on the data
        if protein_inference_parameters.picker:
            data.protein_picker()
        else:
            pass

        ### STEP 10: Apply Inference
        ### STEP 10: Apply Inference
        ### STEP 10: Apply Inference
        inference_type = protein_inference_parameters.inference_type

        # For parsimony... Run GLPK setup, runner, grouper...
        if inference_type == protein_inference.inference.Inference.PARSIMONY:
            group = protein_inference.inference.Parsimony(
                data_class=data, digest_class=digest
            )
            group.infer_proteins()

        if inference_type == protein_inference.inference.Inference.INCLUSION:
            group = protein_inference.inference.Inclusion(
                data_class=data, digest_class=digest
            )
            group.infer_proteins()

        if inference_type == protein_inference.inference.Inference.EXCLUSION:
            group = protein_inference.inference.Exclusion(
                data_class=data, digest_class=digest
            )
            group.infer_proteins()

        ### STEP 11: Run FDR and Q value Calculations
        ### STEP 11: Run FDR and Q value Calculations
        ### STEP 11: Run FDR and Q value Calculations
        data.set_based_fdr()
        data.calculate_q_values()

        # Print the len of restricted data... which is how many protein groups pass FDR threshold
        logger.info(
            "Number of Proteins passing an FDR of "
            + str(protein_inference_parameters.fdr)
            + " = "
            + str(len(data.fdr_restricted_grouped_scored_proteins))
        )

        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        export_type = protein_inference_parameters.export
        export = protein_inference.export.Export(data_class=data)
        export.export_to_csv(
            directory=os.path.join(OUTPUT_DIR, "leads"), export_type=export_type
        )

        logger.info("Protein Inference Finished")

        lead_output = []
        with open(LEAD_OUTPUT_FILE, "r") as lead_output_file:
            reader = csv.reader(lead_output_file, delimiter=",")
            for row in reader:
                lead_output.append(row)

        del lead_output[0]

        protein_groups = data.protein_group_objects

        for i in range(len(protein_groups)):
            lead_protein = protein_groups[i].proteins[0]
            self.assertEqual(lead_protein.identifier, lead_output[i][IDENTIFIER_INDEX])
            self.assertAlmostEqual(
                lead_protein.score, float(lead_output[i][SCORE_INDEX])
            )
            self.assertEqual(
                protein_groups[i].q_value, float(lead_output[i][Q_VALUE_INDEX])
            )
            self.assertEqual(
                protein_groups[i].number_id, int(lead_output[i][GROUP_ID_INDEX])
            )
            self.assertEqual(
                lead_protein.peptides, set(lead_output[i][PEPTIDES_INDEX:])
            )

        export.export_to_csv(
            directory=os.path.join(OUTPUT_DIR, "all"), export_type="q_value_all"
        )

        all_output = []
        with open(ALL_OUTPUT_FILE, "r") as all_output_file:
            reader = csv.reader(all_output_file, delimiter=",")
            for row in reader:
                all_output.append(row)

        del all_output[0]

        all_output_new = []
        with open(export.filepath, "r") as all_output_file_new:
            reader = csv.reader(all_output_file_new, delimiter=",")
            for row in reader:
                all_output_new.append(row)

        del all_output_new[0]

        for i in range(len(all_output)):
            self.assertEqual(
                all_output_new[i][IDENTIFIER_INDEX], all_output[i][IDENTIFIER_INDEX]
            )
            self.assertAlmostEqual(
                float(all_output_new[i][SCORE_INDEX]), float(all_output[i][SCORE_INDEX])
            )
            self.assertEqual(
                float(all_output_new[i][Q_VALUE_INDEX]),
                float(all_output[i][Q_VALUE_INDEX]),
            )
            self.assertEqual(
                int(all_output_new[i][GROUP_ID_INDEX]),
                int(all_output[i][GROUP_ID_INDEX]),
            )
            self.assertEqual(
                set(all_output_new[i][PEPTIDES_INDEX:]),
                set(all_output[i][PEPTIDES_INDEX:]),
            )

        export.export_to_csv(
            directory=os.path.join(OUTPUT_DIR, "peptides"), export_type="peptides"
        )

        peptide_output = []
        with open(PEPTIDE_OUTPUT_FILE, "r") as peptide_output_file:
            reader = csv.reader(peptide_output_file, delimiter=",")
            for row in reader:
                peptide_output.append(row)

        del peptide_output[0]

        peptide_output_new = []
        with open(export.filepath, "r") as peptide_output_file_new:
            reader = csv.reader(peptide_output_file_new, delimiter=",")
            for row in reader:
                peptide_output_new.append(row)

        del peptide_output_new[0]

        for i in range(len(peptide_output)):
            self.assertEqual(
                peptide_output_new[i][IDENTIFIER_INDEX],
                peptide_output[i][IDENTIFIER_INDEX],
            )
            self.assertAlmostEqual(
                float(peptide_output_new[i][SCORE_INDEX]),
                float(peptide_output[i][SCORE_INDEX]),
            )
            self.assertEqual(
                float(peptide_output_new[i][Q_VALUE_INDEX]),
                float(peptide_output[i][Q_VALUE_INDEX]),
            )
            self.assertEqual(
                int(peptide_output_new[i][GROUP_ID_INDEX]),
                int(peptide_output[i][GROUP_ID_INDEX]),
            )
            self.assertEqual(
                set(peptide_output_new[i][PEPTIDES_INDEX:]),
                set(peptide_output[i][PEPTIDES_INDEX:]),
            )

        export.export_to_csv(
            directory=os.path.join(OUTPUT_DIR, "psms"), export_type="psms"
        )

        psm_output = []
        with open(PSM_OUTPUT_FILE, "r") as psm_output_file:
            reader = csv.reader(psm_output_file, delimiter=",")
            for row in reader:
                psm_output.append(row)

        del psm_output[0]

        psm_output_new = []
        with open(export.filepath, "r") as psm_output_file_new:
            reader = csv.reader(psm_output_file_new, delimiter=",")
            for row in reader:
                psm_output_new.append(row)

        del psm_output_new[0]

        for i in range(len(psm_output)):
            self.assertEqual(
                psm_output_new[i][IDENTIFIER_INDEX], psm_output[i][IDENTIFIER_INDEX]
            )
            self.assertAlmostEqual(
                float(psm_output_new[i][SCORE_INDEX]), float(psm_output[i][SCORE_INDEX])
            )
            self.assertEqual(
                float(psm_output_new[i][Q_VALUE_INDEX]),
                float(psm_output[i][Q_VALUE_INDEX]),
            )
            self.assertEqual(
                int(psm_output_new[i][GROUP_ID_INDEX]),
                int(psm_output[i][GROUP_ID_INDEX]),
            )
            self.assertEqual(
                set(psm_output_new[i][PEPTIDES_INDEX:]),
                set(psm_output[i][PEPTIDES_INDEX:]),
            )

        export.export_to_csv(
            directory=os.path.join(OUTPUT_DIR, "psm_ids"), export_type="psm_ids"
        )

        psm_id_output = []
        with open(PSM_ID_OUTPUT_FILE, "r") as psm_id_output_file:
            reader = csv.reader(psm_id_output_file, delimiter=",")
            for row in reader:
                psm_id_output.append(row)

        del psm_id_output[0]

        psm_id_output_new = []
        with open(export.filepath, "r") as psm_id_output_file_new:
            reader = csv.reader(psm_id_output_file_new, delimiter=",")
            for row in reader:
                psm_id_output_new.append(row)

        del psm_id_output_new[0]

        for i in range(len(psm_id_output)):
            self.assertEqual(
                psm_id_output_new[i][IDENTIFIER_INDEX],
                psm_id_output[i][IDENTIFIER_INDEX],
            )
            self.assertAlmostEqual(
                float(psm_id_output_new[i][SCORE_INDEX]),
                float(psm_id_output[i][SCORE_INDEX]),
            )
            self.assertEqual(
                float(psm_id_output_new[i][Q_VALUE_INDEX]),
                float(psm_id_output[i][Q_VALUE_INDEX]),
            )
            self.assertEqual(
                int(psm_id_output_new[i][GROUP_ID_INDEX]),
                int(psm_id_output[i][GROUP_ID_INDEX]),
            )
            self.assertEqual(
                set(psm_id_output_new[i][PEPTIDES_INDEX:]),
                set(psm_id_output[i][PEPTIDES_INDEX:]),
            )

        pipeline = protein_inference.pipeline.ProteinInferencePipeline(
            parameter_file=PARAMETER_FILE,
            database_file=TEST_DATABASE,
            target_files=TARGET_FILE,
            decoy_files=DECOY_FILE,
            combined_files=None,
            output_directory=os.path.join(OUTPUT_DIR, "leads"),
            id_splitting=True,
        )

        pipeline.execute()

        pipeline_protein_groups = pipeline.data.protein_group_objects

        for i in range(len(pipeline_protein_groups)):
            lead_protein_pipeline = pipeline_protein_groups[i].proteins[0]
            lead_protein = protein_groups[i].proteins[0]
            self.assertEqual(lead_protein_pipeline.identifier, lead_protein.identifier)
            self.assertAlmostEqual(lead_protein_pipeline.score, lead_protein.score)
            self.assertEqual(
                pipeline_protein_groups[i].q_value, protein_groups[i].q_value
            )
            self.assertEqual(
                pipeline_protein_groups[i].number_id, protein_groups[i].number_id
            )
            self.assertEqual(lead_protein_pipeline.peptides, lead_protein.peptides)
