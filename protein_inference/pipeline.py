import os
import sys
import logging
import protein_inference

class ProteinInferencePipeline(object):

    def __init__(self, parameter_file, database_file, target_files=None, decoy_files=None, combined_files=None, target_directory=None, decoy_directory=None, combined_directory=None, output_directory=None):


        self.logger = logging.getLogger("protein_inference.pipeline.ProteinInferencePipeline")

        # set up our logger
        logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        self.parameter_file = parameter_file
        self.database_file = database_file
        self.target_files = target_files
        self.decoy_files = decoy_files
        self.combined_files = combined_files
        self.target_directory = target_directory
        self.decoy_directory = decoy_directory
        self.combined_directory = combined_directory
        self.output_directory = output_directory
        self.data = None
        self.digest = None

        self._validate_input()

        self._set_output_directory()

    def execute(self):

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = protein_inference.parameters.ProteinInferenceParameter(yaml_param_filepath=self.parameter_file)

        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = protein_inference.in_silico_digest.InSilicoDigest(database_path=self.database_file,
                                                 parameter_file_object=protein_inference_parameters, id_splitting=True)
        digest.digest_fasta_database()

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data = protein_inference.reader.GenericReader(target_file=self.target_files,
                                                                   decoy_file=self.decoy_files,
                                                                   combined_files=self.combined_files,
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
        if inference_type == 'parsimony':
            group = protein_inference.inference.Parsimony(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "inclusion":
            group = protein_inference.inference.Inclusion(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "exclusion":
            group = protein_inference.inference.Exclusion(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "none":
            group = protein_inference.inference.FirstProtein(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "peptide_centric":
            group = protein_inference.inference.PeptideCentric(data_class=data, digest_class=digest)
            group.infer_proteins()

        ### STEP 11: Run FDR and Q value Calculations
        ### STEP 11: Run FDR and Q value Calculations
        ### STEP 11: Run FDR and Q value Calculations
        data.set_based_fdr(false_discovery_rate=float(protein_inference_parameters.fdr))
        data.calculate_q_values()

        # Print the len of restricted data... which is how many protein groups pass FDR threshold
        self.logger.info('Number of Proteins passing an FDR of ' + str(protein_inference_parameters.fdr) + ' = ' + str(
            len(data.fdr_restricted_grouped_scored_proteins)))

        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        export_type = protein_inference_parameters.export
        export = protein_inference.export.Export(data_class=data)
        export.export_to_csv(directory=self.output_directory, export_type=export_type)

        self.data = data
        self.digest = digest

        self.logger.info('Protein Inference Finished')


    def _validate_input(self):

        if self.target_files and self.decoy_files and not self.combined_files and not self.target_directory and not self.decoy_directory and not self.combined_directory:
            self.logger.info("Validating input as target_files and decoy_files")
        elif self.combined_files and not self.target_files and not self.decoy_files and not self.decoy_directory and not self.target_directory and not self.combined_directory:
            self.logger.info("Validating input as combined_files")
        elif self.target_directory and self.decoy_directory and not self.target_files and not self.decoy_files and not self.combined_directory and not self.combined_files:
            self.logger.info("Validating input as target_directory and decoy_directory")
            self._transform_directory_to_files()
        elif self.combined_directory and not self.combined_files and not self.decoy_files and not self.decoy_directory and not self.target_files and not self.target_directory:
            self.logger.info("Validating input as combined_directory")
            self._transform_directory_to_files()
        else:
            raise ValueError("To run Protein inference please supply either: "
                             "(1) either one or multiple target_files and decoy_files, "
                             "(2) either one or multiple combined_files that include target and decoy data"
                             "(3) a directory that contains target files (target_directory) as well as a directory that contains decoy files (decoy_directory)"
                             "(4) a directory that contains combined target/decoy files (combined_directory)")


    def _transform_directory_to_files(self):

        if self.target_directory and self.decoy_directory:
            self.logger.info("Transforming target_directory and decoy_directory into files")
            target_files = os.listdir(self.target_directory)
            target_files_full = [os.path.join(self.target_directory,x) for x in target_files if x.endswith(".txt") or x.endswith(".tsv")]

            decoy_files = os.listdir(self.decoy_directory)
            decoy_files_full = [os.path.join(self.decoy_directory,x) for x in decoy_files if x.endswith(".txt") or x.endswith(".tsv")]

            self.target_files = target_files_full
            self.decoy_files = decoy_files_full

        elif self.combined_directory:
            self.logger.info("Transforming combined_directory into files")
            combined_files = os.listdir(self.combined_directory)
            combined_files_full = [os.path.join(self.combined_directory,x) for x in combined_files if x.endswith(".txt") or x.endswith(".tsv")]
            self.combined_files = combined_files_full


    def _set_output_directory(self):
        if not self.output_directory:
            self.output_directory = os.getcwd()
        else:
            pass