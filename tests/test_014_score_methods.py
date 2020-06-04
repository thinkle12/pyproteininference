import tempfile
from unittest import TestCase
from pkg_resources import resource_filename

import protein_inference

TEST_DATABASE = resource_filename('protein_inference', '../tests/data/test_database.fasta')
PARAMETER_FILE = resource_filename('protein_inference', '../tests/data/test_params_inclusion.yaml')
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')

TARGET_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_target.txt')
DECOY_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_decoy.txt')

temp_dir = tempfile.gettempdir()

class TestScoreMethods(TestCase):

    def test_scoring_methods(self):
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = protein_inference.parameters.ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE)

        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = protein_inference.in_silico_digest.InSilicoDigest(database_path=TEST_DATABASE,
                                                                   parameter_file_object=protein_inference_parameters,
                                                                   id_splitting=True)
        digest.digest_fasta_database()

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data = protein_inference.reader.GenericReader(target_file=TARGET_FILE,
                                                                   decoy_file=DECOY_FILE,
                                                                   combined_files=None,
                                                                   parameter_file_object=protein_inference_parameters,
                                                                   digest_class=digest)
        pep_and_prot_data.read_psms()

        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data = protein_inference.datastore.DataStore(pep_and_prot_data, digest_class=digest)

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters)

        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        data.create_scoring_input(score_input=protein_inference_parameters.score)

        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters.inference_type == "exclusion":
            # This gets ran if we run exclusion...
            data.exclude_non_distinguishing_peptides(digest_class=digest)

        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        score = protein_inference.scoring.Score(data_class=data)
        score.score_psms(score_method='best_peptide_per_protein')
        self.assertEqual(data.scored_proteins[0].score,3.5e-06)

        score.score_psms(score_method='iterative_downweighted_log')
        self.assertEqual(data.scored_proteins[0].score,64.16418132258529)

        score.score_psms(score_method='multiplicative_log')
        self.assertEqual(data.scored_proteins[0].score,82.89306334778564)

        score.score_psms(score_method='downweighted_multiplicative_log')
        self.assertEqual(data.scored_proteins[0].score,38.42870136439401)

        score.score_psms(score_method='downweighted_version2')
        self.assertEqual(data.scored_proteins[0].score,29.723342219506584)

        score.score_psms(score_method='top_two_combined')
        self.assertEqual(data.scored_proteins[0].score,25.523842819277114)

        score.score_psms(score_method='geometric_mean')
        self.assertEqual(data.scored_proteins[0].score,11.790572089099232)