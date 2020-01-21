import logging
import protein_inference

class ProteinInferencePipeline(object):

    def __init__(self, parameter_file, database_file, target_files=None, decoy_files=None, files=None, output_directory=None):
        self.parameter_file = parameter_file
        self.database_file = database_file
        self.target_files = target_files
        self.decoy_files = decoy_files
        # Add an option to just have files... which is target/decoy already combined...
        self.files = files
        self.output_directory = output_directory

    def execute(self):

        logging.basicConfig(level=logging.INFO)

        logger = logging.getLogger("protein_inference.scripts.Command_Line_PI_Runner_Yaml")

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = protein_inference.ProteinInferenceParameter(yaml_param_filepath=self.parameter_file)

        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = protein_inference.in_silico_digest.InSilicoDigest(database_path=self.database_file,
                                                 parameter_file_object=protein_inference_parameters, id_splitting=True)
        digest.digest_fasta_database()

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        # TODO for the Reader functions add an option to have a combined target/decoy input file instead of having target/decoy files separate
        pep_and_prot_data = protein_inference.reader.GenericReader(target_file=self.target_files,
                                                                      decoy_file=self.decoy_files,
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
        logger.info('Number of Proteins passing an FDR of' + str(protein_inference_parameters.fdr) + ' = ' + str(
            len(data.fdr_restricted_grouped_scored_proteins)))

        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        ### STEP 12: Export to CSV
        export_type = protein_inference_parameters.export
        export = protein_inference.export.Export(data_class=data)
        export.export_to_csv(directory=self.output_directory, export_type=export_type)

        logger.info('Protein Inference Finished')


    def _validate_input(self):
        # TODO write an internal function to validate the input
        pass



