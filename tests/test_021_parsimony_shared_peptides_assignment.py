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
PARAMETER_FILE = resource_filename("pyproteininference", "../tests/data/test_params_parsimony_pulp.yaml")
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('pyproteininference', '../tests/output/')
for sub_dir in ["leads", "all", "peptides", "psms", "psm_ids"]:
    if not os.path.exists(os.path.join(OUTPUT_DIR, sub_dir)):
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))


IDENTIFIER_INDEX = 0
SCORE_INDEX = 1
Q_VALUE_INDEX = 2
GROUP_ID_INDEX = 5
PEPTIDES_INDEX = 6


class TestLoadParsimonyPulpWorkflowSharedPeptideReassignment(TestCase):
    def test_workflow_parsimony_pulp_shared_peptides_best(self):

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        # Set the shared_peptides param to "best"
        protein_inference_parameters.shared_peptides = "best"

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

        lead_protein_psms = [x[0].psms for x in data.grouped_scored_proteins]

        # Get the PSM Identifiers from our lead proteins
        psm_id_list = []
        for protein in lead_protein_psms:
            for psm in protein:
                psm_id_list.append(psm.identifier)

        # Make sure the there are no duplicate psms in our lead proteins
        self.assertEqual(len(psm_id_list), len(set(psm_id_list)))

    def test_workflow_parsimony_pulp_shared_peptides_all(self):

        # Do the same analysis but put shared peptides everywhere

        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        # STEP 1: Load parameter file #
        protein_inference_parameters_all = ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        # Set the shared_peptides param to "best"
        protein_inference_parameters_all.shared_peptides = "all"

        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        # STEP 2: Start with running an In Silico Digestion #
        digest_all = in_silico_digest.PyteomicsDigest(
            database_path=TEST_DATABASE,
            digest_type=protein_inference_parameters_all.digest_type,
            missed_cleavages=protein_inference_parameters_all.missed_cleavages,
            reviewed_identifier_symbol=protein_inference_parameters_all.reviewed_identifier_symbol,
            max_peptide_length=protein_inference_parameters_all.restrict_peptide_length,
            id_splitting=True,
        )
        digest_all.digest_fasta_database()

        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        # STEP 3: Read PSM Data #
        pep_and_prot_data_all = pyproteininference.reader.GenericReader(
            target_file=TARGET_FILE,
            decoy_file=DECOY_FILE,
            parameter_file_object=protein_inference_parameters_all,
            digest=digest_all,
            append_alt_from_db=False,
        )
        pep_and_prot_data_all.read_psms()

        self.assertEqual(len(pep_and_prot_data_all.psms), 27)

        # STEP 4: Initiate the datastore class #
        # STEP 4: Initiate the datastore class #
        # STEP 4: Initiate the datastore class #
        data_all = pyproteininference.datastore.DataStore(pep_and_prot_data_all, digest=digest_all)

        # Step 5: Restrict the PSM data
        # Step 5: Restrict the PSM data
        # Step 5: Restrict the PSM data
        data_all.restrict_psm_data()

        self.assertEqual(len(data_all.main_data_restricted), 26)

        # Step 6: Generate protein scoring input
        # Step 6: Generate protein scoring input
        # Step 6: Generate protein scoring input
        data_all.create_scoring_input()

        # Step 7: Remove non unique peptides if running exclusion
        # Step 7: Remove non unique peptides if running exclusion
        # Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters_all.inference_type == pyproteininference.inference.Inference.EXCLUSION:
            # This gets ran if we run exclusion...
            data_all.exclude_non_distinguishing_peptides()

        # STEP 8: Score our PSMs given a score method
        # STEP 8: Score our PSMs given a score method
        # STEP 8: Score our PSMs given a score method
        score_all = pyproteininference.scoring.Score(data=data_all)
        score_all.score_psms(score_method=protein_inference_parameters_all.protein_score)

        # STEP 9: Run protein picker on the data
        # STEP 9: Run protein picker on the data
        # STEP 9: Run protein picker on the data
        if protein_inference_parameters_all.picker:
            data_all.protein_picker()
        else:
            pass

        # STEP 10: Apply Inference
        # STEP 10: Apply Inference
        # STEP 10: Apply Inference
        inference_type = protein_inference_parameters_all.inference_type

        if inference_type == pyproteininference.inference.Inference.PARSIMONY:
            group = pyproteininference.inference.Parsimony(data=data_all, digest=digest_all)
            group.infer_proteins()

        if inference_type == pyproteininference.inference.Inference.INCLUSION:
            group = pyproteininference.inference.Inclusion(data=data_all, digest=digest_all)
            group.infer_proteins()

        if inference_type == pyproteininference.inference.Inference.EXCLUSION:
            group = pyproteininference.inference.Exclusion(data=data_all, digest=digest_all)
            group.infer_proteins()

        lead_protein_psms = [x[0].psms for x in data_all.grouped_scored_proteins]

        # Get the PSM Identifiers from our lead proteins
        psm_id_list = []
        for protein in lead_protein_psms:
            for psm in protein:
                psm_id_list.append(psm.identifier)

        # Make sure the there are duplicate psms in our lead protein objects
        self.assertNotEqual(len(psm_id_list), len(set(psm_id_list)))
