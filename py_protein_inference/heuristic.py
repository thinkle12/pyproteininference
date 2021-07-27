import os
import sys
import numpy
import logging
import py_protein_inference
from py_protein_inference.pipeline import ProteinInferencePipeline
from py_protein_inference.inference import Inference
import matplotlib.pyplot as plt

class HeuristicPipeline(ProteinInferencePipeline):
    """
    This is the Protein Inference Heuristic classs which houses the logic to run the Protein Inference Heuristic method
     to determine the best inference method for the given data
    Logic is executed in the :py:meth:`py_protein_inference.heuristic.HeuristicPipeline.execute` method

    Attributes:
        parameter_file (str): Path to Protein Inference Yaml Parameter File
        database_file (str): Path to Fasta database used in proteomics search
        target_files (str/list): Path to Target Psm File (Or a list of files)
        decoy_files (str/list): Path to Decoy Psm File (Or a list of files)
        combined_files (str/list): Path to Combined Psm File (Or a list of files)
        target_directory (str): Path to Directory containing Target Psm Files
        decoy_directory (str): Path to Directory containing Decoy Psm Files
        combined_directory (str): Path to Directory containing Combined Psm Files
        output_directory (str): Path to Directory where output will be written
        id_splitting (bool): True/False on whether to split protein IDs in the digest. Leave as False unless you know what you are doing
        append_alt_from_db (bool): True/False on whether to append alternative proteins from the DB digestion in Reader class
        roc_plot_filepath (str): Filepath to be written to by Heuristic Plotting method. This is optional and a default filename will be created in output_directory if this is left as None
        fdr_max (float): The Maximum FDR to display on the ROC Plot generated to compare inference methods
        logger (logging.logger): Logger object for logging
        inference_method_list: (list) List of inference methods used in heuristic determination
        datastore_dict: (dict) Dictionary of :py:class:`py_protein_inference.datastore.DataStore` objects generated in hueristic determination with the inference method as the key of each entry
        selected_method: (str) String representation of the selected inference method based on the heuristic
        heuristic: (float) Heuristic Value as determined from the data
        selected_datastore: (:py:class:`py_protein_inference.datastore.DataStore`) The DataStore object as selected by the heuristic

    """

    def __init__(
            self,
            parameter_file=None,
            database_file=None,
            target_files=None,
            decoy_files=None,
            combined_files=None,
            target_directory=None,
            decoy_directory=None,
            combined_directory=None,
            output_directory=None,
            id_splitting=False,
            append_alt_from_db=True,
            roc_plot_filepath=None,
            fdr_max=0.2,
    ):
        """

        Args:
            parameter_file (str): Path to Protein Inference Yaml Parameter File
            database_file (str): Path to Fasta database used in proteomics search
            target_files (str/list): Path to Target Psm File (Or a list of files)
            decoy_files (str/list): Path to Decoy Psm File (Or a list of files)
            combined_files (str/list): Path to Combined Psm File (Or a list of files)
            target_directory (str): Path to Directory containing Target Psm Files
            decoy_directory (str): Path to Directory containing Decoy Psm Files
            combined_directory (str): Path to Directory containing Combined Psm Files
            output_directory (str): Path to Directory where output will be written
            id_splitting (bool): True/False on whether to split protein IDs in the digest. Leave as False unless you know what you are doing
            append_alt_from_db (bool): True/False on whether to append alternative proteins from the DB digestion in Reader class
            roc_plot_filepath (str): Filepath to be written to by Heuristic Plotting method. This is optional and a default filename will be created in output_directory if this is left as None
            fdr_max (float): The Maximum FDR to display on the ROC Plot generated to compare inference methods    
                
        Returns:
            object:

        Example:
            >>> heuristic = py_protein_inference.heuristic.HeuristicPipeline(
            >>>     parameter_file=yaml_params,
            >>>     database_file=database,
            >>>     target_files=target,
            >>>     decoy_files=decoy,
            >>>     combined_files=combined_files,
            >>>     target_directory=target_directory,
            >>>     decoy_directory=decoy_directory,
            >>>     combined_directory=combined_directory,
            >>>     output_directory=dir_name,
            >>>     append_alt_from_db=append_alt,
            >>>     roc_plot_filepath=roc_plot_filepath,
            >>>     fdr_max=0.2,
            >>> )
        """

        self.logger = logging.getLogger(
            "py_protein_inference.heuristic.HeuristicPipeline"
        )

        # set up our logger
        logging.basicConfig(
            stream=sys.stderr,
            level=logging.INFO,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        )

        self.parameter_file = parameter_file
        self.database_file = database_file
        self.target_files = target_files
        self.decoy_files = decoy_files
        self.combined_files = combined_files
        self.target_directory = target_directory
        self.decoy_directory = decoy_directory
        self.combined_directory = combined_directory
        self.output_directory = output_directory
        self.id_splitting = id_splitting
        self.append_alt_from_db = append_alt_from_db
        if not roc_plot_filepath:
            self.roc_plot_filepath = os.path.join(self.output_directory, "roc_plot.pdf")
        else:
            self.roc_plot_filepath = roc_plot_filepath
        if not fdr_max:
            self.fdr_max = 0.2
        else:
            self.fdr_max = fdr_max

        self.inference_method_list = [Inference.INCLUSION, Inference.EXCLUSION, Inference.PARSIMONY, Inference.PEPTIDE_CENTRIC]
        self.datastore_dict = {}
        self.heuristic = None
        self.selected_method = None
        self.selected_datastore = None

        self._validate_input()

        self._set_output_directory()

        self._log_append_alt_from_db()

    def execute(self, fdr_threshold=0.01, skip_plot=False):
        """
        This method is the main driver of the heuristic method
        This method calls other classes and methods that make up the heuristic pipeline
        This includes but is not limited to:

        1. Loops over the main inference methods: Inclusion, Exclusion, Parsimony, and Peptide Centric
        2. Determines the optimal inference method based on the input data as well as the database file
        3. Outputs the optimal results

        Args:
            fdr_threshold (float): The Qvalue/FDR threshold the heuristic method uses to base calculations from
            skip_plot (bool): True/False on whether to skip outputting ROC Plot

        Returns:
            None

        Example:
            >>> heuristic = py_protein_inference.heuristic.HeuristicPipeline(
            >>>     parameter_file=yaml_params,
            >>>     database_file=database,
            >>>     target_files=target,
            >>>     decoy_files=decoy,
            >>>     combined_files=combined_files,
            >>>     target_directory=target_directory,
            >>>     decoy_directory=decoy_directory,
            >>>     combined_directory=combined_directory,
            >>>     output_directory=dir_name,
            >>>     append_alt_from_db=append_alt,
            >>>     roc_plot_filepath=roc_plot_filepath,
            >>>     fdr_max=0.2,
            >>> )
            >>> heuristic.execute(fdr_threshold=0.01)

        """

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        py_protein_inference_parameters = py_protein_inference.parameters.ProteinInferenceParameter(
            yaml_param_filepath=self.parameter_file
        )

        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = py_protein_inference.in_silico_digest.InSilicoDigest(
            database_path=self.database_file,
            digest_type=py_protein_inference_parameters.digest_type,
            missed_cleavages=py_protein_inference_parameters.missed_cleavages,
            reviewed_identifier_symbol=py_protein_inference_parameters.reviewed_identifier_symbol,
            max_peptide_length=py_protein_inference_parameters.restrict_peptide_length,
            id_splitting=self.id_splitting,
        )
        if self.database_file:
            self.logger.info(
                "Running In Silico Database Digest on file {}".format(
                    self.database_file
                )
            )
            digest.digest_fasta_database()
        else:
            self.logger.warning(
                "No Database File provided, Skipping database digest and only taking protein-peptide mapping from the input files."
            )

        for inference_method in self.inference_method_list:

            py_protein_inference_parameters = py_protein_inference.parameters.ProteinInferenceParameter(
                yaml_param_filepath=self.parameter_file
            )

            self.logger.info("Overriding inference type {}".format(py_protein_inference_parameters.inference_type))

            py_protein_inference_parameters.inference_type = inference_method
            py_protein_inference_parameters.fdr = fdr_threshold

            self.logger.info("New inference type {}".format(py_protein_inference_parameters.inference_type))
            self.logger.info("FDR Threshold Set to {}".format(py_protein_inference_parameters.fdr))


            ### STEP 3: Read PSM Data ###
            ### STEP 3: Read PSM Data ###
            ### STEP 3: Read PSM Data ###
            reader = py_protein_inference.reader.GenericReader(
                target_file=self.target_files,
                decoy_file=self.decoy_files,
                combined_files=self.combined_files,
                parameter_file_object=py_protein_inference_parameters,
                digest_class=digest,
                append_alt_from_db=self.append_alt_from_db,
            )
            reader.read_psms()

            ### STEP 4: Initiate the datastore class ###
            ### STEP 4: Initiate the datastore class ###
            ### STEP 4: Initiate the datastore class ###
            data = py_protein_inference.datastore.DataStore(
                reader_class=reader, digest_class=digest
            )

            ### Step 5: Restrict the PSM data
            ### Step 5: Restrict the PSM data
            ### Step 5: Restrict the PSM data
            data.restrict_psm_data()

            ### Step 6: Generate protein scoring input
            ### Step 6: Generate protein scoring input
            ### Step 6: Generate protein scoring input
            data.create_scoring_input()

            ### Step 7: Remove non unique peptides if running exclusion
            ### Step 7: Remove non unique peptides if running exclusion
            ### Step 7: Remove non unique peptides if running exclusion
            if py_protein_inference_parameters.inference_type == Inference.EXCLUSION:
                # This gets ran if we run exclusion...
                data.exclude_non_distinguishing_peptides()

            ### STEP 8: Score our PSMs given a score method
            ### STEP 8: Score our PSMs given a score method
            ### STEP 8: Score our PSMs given a score method
            score = py_protein_inference.scoring.Score(data_class=data)
            score.score_psms(score_method=py_protein_inference_parameters.protein_score)

            ### STEP 9: Run protein picker on the data
            ### STEP 9: Run protein picker on the data
            ### STEP 9: Run protein picker on the data
            if py_protein_inference_parameters.picker:
                data.protein_picker()
            else:
                pass

            ### STEP 10: Apply Inference
            ### STEP 10: Apply Inference
            ### STEP 10: Apply Inference
            py_protein_inference.inference.Inference.run_inference(data_class=data, digest_class=digest)

            ### STEP 11: Q value Calculations
            ### STEP 11: Q value Calculations
            ### STEP 11: Q value Calculations
            data.calculate_q_values()

            self.datastore_dict[inference_method] = data

        self.selected_method = self.determine_optimal_inference_method()
        self.selected_datastore = self.datastore_dict[self.selected_method]

        export = py_protein_inference.export.Export(data_class=self.selected_datastore)
        export.export_to_csv(directory=self.output_directory, export_type=py_protein_inference_parameters.export)

        if not skip_plot:
            self.generate_roc_plot(fdr_max=self.fdr_max, pdf_filename=self.roc_plot_filepath)

        else:
            self.logger.info("skip_plot is set to True. Not creating ROC Plot.")

    def determine_optimal_inference_method(self, empirical_threshold=.2):

        # Get the number of passing proteins
        filtered_protein_objects = {x:self.datastore_dict[x].get_protein_objects(fdr_restricted=True) for x in self.datastore_dict.keys()}
        number_passing_proteins = {x:len(filtered_protein_objects[x]) for x in filtered_protein_objects.keys()}

        self.logger.info("Number of Passing Proteins per Inference Method")
        self.logger.info(number_passing_proteins)

        # Calculate how similar the number of passing proteins is for each method
        similarity_dict = {}
        for key in number_passing_proteins.keys():
            cur_value = number_passing_proteins[key]
            other_values = [x for x in number_passing_proteins.values() if x != cur_value]
            similarity_dict[key] = cur_value/(numpy.mean(other_values))

        # Simple transformation for getting max below
        diff_dict = {x:abs(1-similarity_dict[x]) for x in number_passing_proteins.keys()}

        self.logger.info("Initial Heuristic Scores")
        self.logger.info(diff_dict)

        # Remove the most dissimilar method, which is the max
        key_to_delete = max(diff_dict, key=lambda k: diff_dict[k])
        self.logger.info("Removing {} with score {}".format(key_to_delete, diff_dict[key_to_delete]))
        del diff_dict[key_to_delete]
        del number_passing_proteins[key_to_delete]

        # Redo above on the restricted set of 3 methods
        similarity_dict2 = {}
        for key in number_passing_proteins.keys():
            cur_value = number_passing_proteins[key]
            other_values = [x for x in number_passing_proteins.values() if x != cur_value]
            similarity_dict2[key] = cur_value/(numpy.mean(other_values))

        diff_dict2 = {x:abs(1-similarity_dict2[x]) for x in number_passing_proteins.keys()}

        self.logger.info("Final Heuristic Scores")
        self.logger.info(diff_dict2)

        # Remove Inclusion and Exclusion if they are poor. IE if their heuristic scores are above an empirical threshold value
        # .2 was determined as a proper threshold in testing different databases (Uniprot, Swissprot, Swissprot no isoforms)
        if Inference.EXCLUSION in diff_dict2.keys():
            if diff_dict2[Inference.EXCLUSION] <= empirical_threshold:
                self.logger.info("Keeping {} with score {}".format(Inference.EXCLUSION, diff_dict2[Inference.EXCLUSION]))

            else:
                # If not the remove it
                self.logger.info("Removing {} with score {}".format(Inference.EXCLUSION, diff_dict2[Inference.EXCLUSION]))
                del diff_dict2[Inference.EXCLUSION]
                del number_passing_proteins[Inference.EXCLUSION]

        if Inference.INCLUSION in diff_dict2.keys():
            if diff_dict2[Inference.INCLUSION] <= empirical_threshold:
                self.logger.info("Keeping {} with score {}".format(Inference.INCLUSION, diff_dict2[Inference.INCLUSION]))

            else:
                # If not then remove it
                self.logger.info("Removing {} with score {}".format(Inference.INCLUSION, diff_dict2[Inference.INCLUSION]))
                del diff_dict2[Inference.INCLUSION]
                del number_passing_proteins[Inference.INCLUSION]

        remaining_inference_methods = list(diff_dict2.keys())

        ## At this point we have 3, 2 or 1 inference types remaining... So lets do branching if statement for all possible combinations
        ## Each combination will have different rules based on empirical knowledgee
        if len(remaining_inference_methods)==3:
            if set([Inference.PARSIMONY, Inference.EXCLUSION, Inference.INCLUSION]) == set(remaining_inference_methods):
                # If inclusion is over double parsimony remove it
                if number_passing_proteins[Inference.INCLUSION] / 2 > number_passing_proteins[Inference.PARSIMONY]:
                    self.logger.info("Removing {} with score {}".format(Inference.INCLUSION, diff_dict2[Inference.INCLUSION]))
                    del diff_dict2[Inference.INCLUSION]
                    del number_passing_proteins[Inference.INCLUSION]

                # If exclusion is less than half of parsimony remove it...
                if number_passing_proteins[Inference.EXCLUSION] * 2 < number_passing_proteins[Inference.PARSIMONY]:
                    self.logger.info("Removing {} with score {}".format(Inference.EXCLUSION, diff_dict2[Inference.EXCLUSION]))
                    del diff_dict2[Inference.EXCLUSION]
                    del number_passing_proteins[Inference.EXCLUSION]

                if len(diff_dict2.keys())==3:
                    # if neither are removed select Exclusion. Since all methods are close to parsimony take exclusion because it wouldnt have removed too many hits
                    selected_method = Inference.EXCLUSION
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method


            elif set([Inference.PARSIMONY, Inference.PEPTIDE_CENTRIC, Inference.INCLUSION]) == set(remaining_inference_methods):
                # If inclusion is double of parsimony remove it...
                if number_passing_proteins[Inference.INCLUSION] / 2 > number_passing_proteins[Inference.PARSIMONY]:
                    self.logger.info("Removing {} with score {}".format(Inference.INCLUSION, diff_dict2[Inference.INCLUSION]))
                    del diff_dict2[Inference.INCLUSION]
                    del number_passing_proteins[Inference.INCLUSION]

                # Check to see if inclusion is still present
                if Inference.INCLUSION in diff_dict2.keys():
                    # If inclusion is still present and double peptide-centric remove it...
                    if number_passing_proteins[Inference.INCLUSION] / 2 > number_passing_proteins[Inference.PEPTIDE_CENTRIC]:
                        self.logger.info("Removing {} with score {}".format(Inference.INCLUSION, diff_dict2[Inference.INCLUSION]))
                        del diff_dict2[Inference.INCLUSION]
                        del number_passing_proteins[Inference.INCLUSION]

                if len(diff_dict2.keys())==3:
                    # if Inclusion is not removed select Inclusion.
                    selected_method = Inference.INCLUSION
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method


            elif set([Inference.PARSIMONY, Inference.EXCLUSION, Inference.PEPTIDE_CENTRIC]) == set(remaining_inference_methods):
                # If exclusion is less than half of parsimony remove it...
                if number_passing_proteins[Inference.EXCLUSION] * 2 < number_passing_proteins[Inference.PARSIMONY]:
                    self.logger.info("Removing {} with score {}".format(Inference.EXCLUSION, diff_dict2[Inference.EXCLUSION]))
                    del diff_dict2[Inference.EXCLUSION]
                    del number_passing_proteins[Inference.EXCLUSION]

                # Check to see if exclusion is still present
                if Inference.EXCLUSION in diff_dict2.keys():
                    # If exclusion is still present and less than half of peptide-centric remove it...
                    if diff_dict2[Inference.EXCLUSION] * 2 < diff_dict2[Inference.PEPTIDE_CENTRIC]:
                        self.logger.info("Removing {} with score {}".format(Inference.EXCLUSION, diff_dict2[Inference.EXCLUSION]))
                        del diff_dict2[Inference.EXCLUSION]
                        del number_passing_proteins[Inference.EXCLUSION]

                if len(diff_dict2.keys())==3:
                    # if Exclusion is not removed select Exclusion. Since it is close to parsimony and peptide-centric take exclusion because it wouldnt have removed too many hits
                    selected_method = Inference.EXCLUSION
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method


            elif set([Inference.PEPTIDE_CENTRIC, Inference.EXCLUSION, Inference.INCLUSION]) == set(remaining_inference_methods):
                # If inclusion is over double peptide-centric remove it
                if number_passing_proteins[Inference.INCLUSION] / 2 > number_passing_proteins[Inference.PEPTIDE_CENTRIC]:
                    self.logger.info("Removing {} with score {}".format(Inference.INCLUSION, diff_dict2[Inference.INCLUSION]))
                    del diff_dict2[Inference.INCLUSION]
                    del number_passing_proteins[Inference.INCLUSION]

                # If exclusion is less than half of peptide-centric remove it...
                if diff_dict2[Inference.EXCLUSION] * 2 < diff_dict2[Inference.PEPTIDE_CENTRIC]:
                    self.logger.info("Removing {} with score {}".format(Inference.EXCLUSION, diff_dict2[Inference.EXCLUSION]))
                    del diff_dict2[Inference.EXCLUSION]
                    del number_passing_proteins[Inference.EXCLUSION]

                if len(diff_dict2.keys())==3:
                    # if neither are removed select Exclusion. Since both are close to peptide-centric take exclusion because it wouldnt have removed too many hits
                    selected_method = Inference.EXCLUSION
                    # If we have one remaining just return it
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method

            else:
                pass

        remaining_inference_methods = list(diff_dict2.keys())

        if len(remaining_inference_methods)==2:
            if set([Inference.PEPTIDE_CENTRIC, Inference.EXCLUSION]) == set(remaining_inference_methods):
                # Take peptide centric
                selected_method = Inference.PEPTIDE_CENTRIC
                self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                return selected_method


            elif set([Inference.PEPTIDE_CENTRIC, Inference.INCLUSION]) == set(remaining_inference_methods):
                # Take peptide centric
                selected_method = Inference.PEPTIDE_CENTRIC
                self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                return selected_method

            elif set([Inference.PEPTIDE_CENTRIC, Inference.PARSIMONY]) == set(remaining_inference_methods):

                # First If peptide centric is less than parsimony return parsimony
                if number_passing_proteins[Inference.PEPTIDE_CENTRIC] < number_passing_proteins[Inference.PARSIMONY]:
                    selected_method = Inference.PARSIMONY
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method

                # check to see if parsimony and peptide centric are close... If half of peptide centric is greater than parsimony pick parsimony
                if number_passing_proteins[Inference.PEPTIDE_CENTRIC] / 2 > number_passing_proteins[Inference.PARSIMONY]:
                    selected_method = Inference.PARSIMONY
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method
                else:
                    # If not (Meaning the values are close... return peptide-centric)
                    selected_method = Inference.PEPTIDE_CENTRIC
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method

            elif set([Inference.EXCLUSION, Inference.INCLUSION]) == set(remaining_inference_methods):
                # This situation should never occur
                selected_method = Inference.INCLUSION
                self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                return selected_method
                pass

            elif set([Inference.EXCLUSION, Inference.PARSIMONY]) == set(remaining_inference_methods):

                # First If exclusion is greater than parsimony return exclusion
                # This is highly unlikely to happen
                if number_passing_proteins[Inference.EXCLUSION] > number_passing_proteins[Inference.PARSIMONY]:
                    selected_method = Inference.EXCLUSION
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method

                # check to see if parsimony and exclusion are close... If half of parsimony is greater than exclusion pick parsimony
                # This means that exclusion is removing a lot of peptides
                if number_passing_proteins[Inference.PARSIMONY] / 2 > number_passing_proteins[Inference.EXCLUSION]:
                    selected_method = Inference.PARSIMONY
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method
                else:
                    # If not (Meaning the values are close... return exclusion)
                    selected_method = Inference.EXCLUSION
                    self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                    return selected_method

            elif set([Inference.PARSIMONY, Inference.INCLUSION]) == set(remaining_inference_methods):
                # take parsimony
                selected_method = Inference.PARSIMONY
                self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
                return selected_method

        remaining_inference_methods = list(diff_dict2.keys())

        if len(remaining_inference_methods)==1:
            selected_method = list(diff_dict2.keys())[0]
            # If we have one remaining just return it
            self.logger.info("Inference {} Selected with score {}".format(selected_method, diff_dict2[selected_method]))
            return selected_method

        if len(remaining_inference_methods)==0:
            raise ValueError("Not able to determine optimal Inference Method for your dataset")

    def generate_roc_plot(self, fdr_max=0.2, pdf_filename=None, target_fdr=None):
        f = plt.figure()
        for inference_method in self.datastore_dict.keys():
            fdr_vs_target_hits = self.datastore_dict[inference_method].generate_fdr_vs_target_hits(fdr_max=fdr_max)
            fdrs = [x[0] for x in fdr_vs_target_hits]
            target_hits = [x[1] for x in fdr_vs_target_hits]
            plt.plot(fdrs, target_hits, '-', label = inference_method)
            target_fdr = self.datastore_dict[inference_method].parameter_file_object.fdr
            if inference_method == self.selected_method:
                best_value = min(fdrs, key=lambda x:abs(x-target_fdr))
                best_index = fdrs.index(best_value)
                best_target_hit_value = target_hits[best_index]

        plt.axvline(target_fdr, color="black", linestyle='--', alpha=.75, label="Target FDR")
        plt.legend()
        plt.xlabel('Decoy FDR')
        plt.ylabel('Target Protein Hits')
        plt.xlim([-.01, fdr_max])
        plt.legend(loc='lower right')
        plt.title("FDR vs Target Protein Hits per Inference Method")
        if pdf_filename:
            f.savefig(pdf_filename)
        plt.show()
        plt.close()

