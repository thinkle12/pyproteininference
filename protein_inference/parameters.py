import yaml
from logging import getLogger

from protein_inference.in_silico_digest import InSilicoDigest
from protein_inference.export import Export
from protein_inference.scoring import Score
from protein_inference.inference import Inference


class ProteinInferenceParameter(object):
    """
    Class that handles data retrieval, storage, and validation of Protein Inference Parameters

    Attributes:
        yaml_param_filepath (str): path to properly formatted parameter file specific to Protein Inference
        digest_type (str): String that determines that type of in silico digestion for :py:class:`protein_inference.in_silico_digest.Digest`. Typically "trypsin"
        export (str): String to indicate the export type for :py:class:`protein_inference.export.Export`. Typically this is "psms", "peptides", or "psm_ids"
        fdr (float): Float to indicate FDR filtering
        glpk_path (str): Path to local installation of glpsol if inference_type="parsimony" and lp_solver="glpk"
        missed_cleavages (int): Integer to determine the number of missed cleavages in the database digestion :py:class:`protein_inference.in_silico_digest.Digest`
        picker (bool): True/False on whether or not to run the protein picker algorithm :py:meth:protein_inference.datastore.DataStore.protein_picker`
        restrict_pep (float/None): Float to restrict the posterior error probability values by in the PSM input. Used in :py:meth:protein_inference.datastore.DataStore.restrict_psm_data`
        restrict_peptide_length (int/None): Float to restrict the peptide length values by in the PSM input. Used in :py:meth:protein_inference.datastore.DataStore.restrict_psm_data`
        restrict_q (float/None): Float to restrict the q values by in the PSM input. Used in :py:meth:protein_inference.datastore.DataStore.restrict_psm_data`
        restrict_custom (float/None): Float to restrict the custom values by in the PSM input. Used in :py:meth:protein_inference.datastore.DataStore.restrict_psm_data`. Filtering depends on score_type variable. If score_type is multiplicative then values that are less than restrict_custom are kept. If score_type is additive then values that are more than restrict_custom are kept.
        protein_score (str): String to determine the way in which Proteins are scored can be any of the SCORE_METHODS in :py:class:`protein_inference.scoring.Score`
        psm_score_type (str): String to determine the type of score that the PSM scores are (Additive or Multiplicative) can be any of the SCORE_TYPES in :py:class:`protein_inference.scoring.Score`
        decoy_symbol (str): String to denote decoy proteins from target proteins. IE "##"
        isoform_symbol (str): String to denote isoforms from regular proteins. IE "-". Can also be None
        reviewed_identifier_symbol (str): String to denote a "Reviewed" Protein. Typically this is: "sp|" if using Uniprot Fasta database
        inference_type (str): String to determine the inference procedure. Can be any value of INFERENCE_TYPES of :py:class:`protein_inference.inference.Inference` object
        tag (str): String to be added to output files
        psm_score (str): String that indicates the PSM input score. The value should match the string in the input data of the score you want to use for PSM score. This score will be used in scoring methods here: :py:class:`protein_inference.scoring.Score`
        grouping_type (str/None): String to determine the grouping procedure. Can be any value of GROUPING_TYPES of :py:class:`protein_inference.inference.Inference` object
        max_identifiers_peptide_centric (int): Maximum number of identifiers to assign to a group when running peptide_centric inference. Typically this is 10 or 5.
        lp_solver (str/None): The LP solver to use if inference_type="Parsimony". Can be any value in LP_SOLVERS in the :py:class:`protein_inference.inference.Inference` object
        logger (logger.logging): Logger object

    """

    def __init__(self, yaml_param_filepath, validate=True):
        """ Class to store Protein Inference parameter information as an object

        Args:
            yaml_param_filepath (str): path to properly formatted parameter file specific to Protein Inference
            validate (bool): True/False on whether to validate the parameter file of interest

        Returns:
            None

        Example:
            >>> protein_inference.parameters.ProteinInferenceParameter(
            >>>     yaml_param_filepath = "/path/to/protein_inference_params.yaml", validate=True
            >>> )


        """
        self.yaml_param_filepath = yaml_param_filepath
        self.digest_type = None
        self.export = None
        self.fdr = None
        self.glpk_path = None
        self.missed_cleavages = None
        self.picker = None
        self.restrict_pep = None
        self.restrict_peptide_length = None
        self.restrict_q = None
        self.restrict_custom = None
        self.protein_score = None
        self.psm_score_type = None
        self.decoy_symbol = None
        self.isoform_symbol = None
        self.reviewed_identifier_symbol = None
        self.inference_type = None
        self.tag = None
        self.psm_score = None
        self.grouping_type = None
        self.max_identifiers_peptide_centric = None
        self.lp_solver = None
        self.logger = getLogger(
            "protein_inference.parameters.ProteinInferenceParameter.validate_parameters"
        )

        self.convert_to_object()

        if validate:
            self.validate_parameters()

        self._fix_none_parameters()

    def convert_to_object(self):
        """
        Function that takes a Protein Inference parameter file and converts it into a ProteinInferenceParameter object
        by assigning all Attributes of the ProteinInferenceParameter object

        If no parameter filepath is supplied the parameter object will be loaded with default params

        This function gets ran in the initilization of the ProteinInferenceParameter object

        Args:
            None

        Returns:
            None

        """
        if self.yaml_param_filepath:
            with open(self.yaml_param_filepath, "r") as stream:
                yaml_params = yaml.load(stream, Loader=yaml.Loader)

            self.digest_type = yaml_params["parameters"]["digest"]["digest_type"]
            self.export = yaml_params["parameters"]["general"]["export"]
            self.fdr = yaml_params["parameters"]["general"]["fdr"]
            self.glpk_path = yaml_params["parameters"]["parsimony"]["glpk_path"]
            self.missed_cleavages = yaml_params["parameters"]["digest"][
                "missed_cleavages"
            ]
            self.picker = yaml_params["parameters"]["general"]["picker"]
            self.restrict_pep = yaml_params["parameters"]["data_restriction"][
                "pep_restriction"
            ]
            self.restrict_peptide_length = yaml_params["parameters"][
                "data_restriction"
            ]["peptide_length_restriction"]
            self.restrict_q = yaml_params["parameters"]["data_restriction"][
                "q_value_restriction"
            ]
            self.restrict_custom = yaml_params["parameters"]["data_restriction"][
                "custom_restriction"
            ]
            self.protein_score = yaml_params["parameters"]["score"]["protein_score"]
            self.psm_score_type = yaml_params["parameters"]["score"]["psm_score_type"]
            self.decoy_symbol = yaml_params["parameters"]["identifiers"]["decoy_symbol"]
            self.isoform_symbol = yaml_params["parameters"]["identifiers"][
                "isoform_symbol"
            ]
            self.reviewed_identifier_symbol = yaml_params["parameters"]["identifiers"][
                "reviewed_identifier_symbol"
            ]
            self.inference_type = yaml_params["parameters"]["inference"][
                "inference_type"
            ]
            self.tag = yaml_params["parameters"]["general"]["tag"]
            self.psm_score = yaml_params["parameters"]["score"]["psm_score"]
            self.grouping_type = yaml_params["parameters"]["inference"]["grouping_type"]
            self.max_identifiers_peptide_centric = yaml_params["parameters"][
                "peptide_centric"
            ]["max_identifiers"]
            self.lp_solver = yaml_params["parameters"]["parsimony"]["lp_solver"]
            try:
                # Do try except here to make old param files backwards compatible
                self.shared_peptides = yaml_params["parameters"]["parsimony"]["shared_peptides"]
            except KeyError:
                self.shared_peptides = Inference.ALL_SHARED_PEPTIDES

        else:
            self.logger.warning(
                "Yaml parameter file not found, parameters set to default"
            )
            self.digest_type = "trypsin"
            self.export = "q_value"
            self.fdr = 0.01
            self.glpk_path = "glpsol"
            self.missed_cleavages = 3
            self.picker = True
            self.restrict_pep = 0.9
            self.restrict_peptide_length = 7
            self.restrict_q = 0.005
            self.restrict_custom = None
            self.protein_score = "multiplicative_log"
            self.psm_score = "posterior_error_prob"
            self.decoy_symbol = "##"
            self.isoform_symbol = "-"
            self.reviewed_identifier_symbol = "sp|"
            self.inference_type = "peptide_centric"
            self.tag = "example_tag"
            self.psm_score_type = "multiplicative"
            self.grouping_type = "shared_peptides"
            self.max_identifiers_peptide_centric = 5
            self.lp_solver = "pulp"

    def validate_parameters(self):
        """
        Class method to validate all parameters

        Returns:
            None

        """
        # Run all of the parameter validations
        self._validate_digest_type()
        self._validate_export_type()
        self._validate_floats()
        self._validate_bools()
        self._validate_score_type()
        self._validate_score_method()
        self._validate_score_combination()
        self._validate_inference_type()
        self._validate_grouping_type()
        self._validate_max_id()
        self._validate_lp_solver()
        self._validate_identifiers()

    def _validate_digest_type(self):
        """
        Internal ProteinInferenceParameter method to validate the digest type
        """
        # Make sure we have a valid digest type
        if self.digest_type in InSilicoDigest.LIST_OF_DIGEST_TYPES:
            self.logger.info("Using digest type '{}'".format(self.digest_type))
        else:
            raise ValueError(
                "Digest Type '{}' not supported, please use one of the following enyzme digestions: '{}'".format(
                    self.digest_type, ", ".join(InSilicoDigest.LIST_OF_DIGEST_TYPES)
                )
            )

    def _validate_export_type(self):
        """
        Internal ProteinInferenceParameter method to validate the export type
        """
        # Make sure we have a valid export type
        if self.export in Export.EXPORT_TYPES:
            self.logger.info("Using Export type '{}'".format(self.export))
        else:
            raise ValueError(
                "Export Type '{}' not supported, please use one of the following export types: '{}'".format(
                    self.export, ", ".join(Export.EXPORT_TYPES)
                )
            )
        pass

    def _validate_floats(self):
        """
        Internal ProteinInferenceParameter method to validate floats
        """
        # Validate that FDR, cleavages, and restrict values are all floats and or ints if they need to be

        try:
            if 0 <= float(self.fdr) <= 1:
                self.logger.info("FDR Input {}".format(self.fdr))

        except ValueError:
            raise ValueError(
                "FDR must be a decimal between 0 and 1, FDR provided: {}".format(
                    self.fdr
                )
            )

        try:
            if 0 <= float(self.restrict_pep) <= 1:
                self.logger.info("PEP restriction {}".format(self.restrict_pep))

        except ValueError:
            if not self.restrict_pep or self.restrict_pep.lower() == "none":
                self.restrict_pep = None
                self.logger.info("Not restrict by PEP Value")
            else:
                raise ValueError(
                    "PEP restriction must be a decimal between 0 and 1, PEP restriction provided: {}".format(
                        self.restrict_pep
                    )
                )

        try:
            if 0 <= float(self.restrict_q) <= 1:
                self.logger.info("Q Value restriction {}".format(self.restrict_q))

        except ValueError:
            if not self.restrict_q or self.restrict_q.lower() == "none":
                self.restrict_q = None
                self.logger.info("Not restrict by Q Value")
            else:
                raise ValueError(
                    "Q Value restriction must be a decimal between 0 and 1, Q Value restriction provided: {}".format(
                        self.restrict_q
                    )
                )

        try:
            int(self.missed_cleavages)
            self.logger.info(
                "Missed Cleavages selected: {}".format(self.missed_cleavages)
            )
        except ValueError:
            raise ValueError(
                "Missed Cleavages must be an integer, Provided Missed Cleavages value: {}".format(
                    self.missed_cleavages
                )
            )

        try:
            int(self.restrict_peptide_length)
            self.logger.info(
                "Peptide Length Restriction: Len {}".format(
                    self.restrict_peptide_length
                ))
        except ValueError:
            if not self.restrict_peptide_length or self.restrict_peptide_length.lower() == "none":
                self.restrict_peptide_length = None
                self.logger.info("Not Restricting by Peptide Length")
            else:
                raise ValueError(
                    "Peptide Length Restriction must be an integer, Provided Peptide Length Restriction value: {}".format(
                        self.restrict_peptide_length
                    )
            )

        try:
            float(self.restrict_custom)
            self.logger.info("Custom restriction {}".format(self.restrict_custom))
        except ValueError:
            if not self.restrict_custom or self.restrict_custom.lower() == "none":
                self.restrict_custom = None
                self.logger.info("Not Restricting by Custom Value")
            else:
                raise ValueError(
                    "Custom restriction must be a number, Custom restriction provided: {}".format(
                        self.restrict_custom
                    )
                )

    def _validate_bools(self):
        """
        Internal ProteinInferenceParameter method to validate the bools
        """
        # Make sure picker is a bool
        if type(self.picker) == bool:
            if self.picker:
                self.logger.info("Parameters loaded to run Picker")
            else:
                self.logger.info("Parameters loaded to NOT run Picker")
        else:
            raise ValueError(
                "Picker Variable must be set to True or False, Picker Variable provided: {}".format(
                    self.picker
                )
            )

    def _validate_score_method(self):
        """
        Internal ProteinInferenceParameter method to validate the score method
        """
        # Make sure we have the score method defined in code to use...
        if self.protein_score in Score.SCORE_METHODS:
            self.logger.info("Using Score Method '{}'".format(self.protein_score))
        else:
            raise ValueError(
                "Score Method '{}' not supported, "
                "please use one of the following Score Methods: '{}'".format(
                    self.protein_score, ", ".join(Score.SCORE_METHODS)
                )
            )

    def _validate_score_type(self):
        """
        Internal ProteinInferenceParameter method to validate the score type
        """
        # Make sure score type is multiplicative or additive
        if self.psm_score_type in Score.SCORE_TYPES:
            self.logger.info("Using Score Type '{}'".format(self.psm_score_type))
        else:
            raise ValueError(
                "Score Type '{}' not supported, "
                "please use one of the following Score Types: '{}'".format(
                    self.psm_score_type, ", ".join(Score.SCORE_TYPES)
                )
            )

    def _validate_score_combination(self):
        """
        Internal ProteinInferenceParameter method to validate combination of score method and score type
        """
        # Check to see if combination of score (column), method(multiplicative log, additive), and score type (multiplicative/additive) is possible...
        # This will be super custom

        if self.psm_score_type == Score.ADDITIVE_SCORE_TYPE and self.protein_score != Score.ADDITIVE:
            raise ValueError(
                "If Score type is 'additive' (Higher PSM score is better) then you must use the 'additive' score method"
            )

        elif self.psm_score_type == Score.MULTIPLICATIVE_SCORE_TYPE and self.protein_score == Score.ADDITIVE:
            raise ValueError(
                "If Score type is 'multiplicative' (Lower PSM score is better) "
                "then you must NOT use the 'additive' score method please "
                "select one of the following score methods: {}".format(
                    ", ".join([x for x in Score.SCORE_METHODS if x != "additive"])
                )
            )

        else:
            self.logger.info(
                "Combination of Score Type: '{}' and Score Method: '{}' is Ok".format(
                    self.psm_score_type, self.protein_score
                )
            )

    def _validate_inference_type(self):
        """
        Internal ProteinInferenceParameter method to validate the inference type
        """
        # Check if its parsimony, exclusion, inclusion, none
        if self.inference_type in Inference.INFERENCE_TYPES:
            self.logger.info("Using inference type '{}'".format(self.inference_type))
        else:
            raise ValueError(
                "Inferece Type '{}' not supported, please use one of the following Inferece Types: '{}'".format(
                    self.inference_type, ", ".join(Inference.INFERENCE_TYPES)
                )
            )

    def _validate_grouping_type(self):
        """
        Internal ProteinInferenceParameter method to validate the grouping type
        """
        # Check if its parsimony, exclusion, inclusion, none
        if self.grouping_type in Inference.GROUPING_TYPES:
            self.logger.info("Using Grouping type '{}'".format(self.grouping_type))
        else:
            if self.grouping_type.lower() == "none" or not self.grouping_type:
                self.grouping_type = None
                self.logger.info("Using Grouping type: None")
            else:

                raise ValueError(
                    "Grouping Type '{}' not supported, please use one of the following Grouping Types: '{}'".format(
                        self.grouping_type, ", ".join(Inference.GROUPING_TYPES)
                    )
                )

    def _validate_max_id(self):
        """
        Internal ProteinInferenceParameter method to validate the max peptide centric id
        """
        # Check if max_identifiers_peptide_centric param is an INT
        if type(self.max_identifiers_peptide_centric) == int:
            self.logger.info(
                "Max Number of Indentifiers for Peptide Centric Inference: '{}'".format(
                    self.max_identifiers_peptide_centric
                )
            )
        else:
            raise ValueError(
                "Max Number of Indentifiers for Peptide Centric Inference must be an integer, provided value: {}".format(
                    self.max_identifiers_peptide_centric
                )
            )

    def _validate_lp_solver(self):
        """
        Internal ProteinInferenceParameter method to validate the lp solver
        """
        # Check if its pulp, glpk, or None
        if self.lp_solver in Inference.LP_SOLVERS:
            self.logger.info("Using LP Solver '{}'".format(self.lp_solver))
        else:
            if self.lp_solver.lower() == "none" or not self.lp_solver:
                self.lp_solver = None
                self.logger.info("Setting LP Solver to None")
            else:
                raise ValueError(
                    "LP Solver '{}' not supported, please use one of the following LP Solvers: '{}'".format(
                        self.lp_solver, ", ".join(Inference.LP_SOLVERS)
                    )
                )

    def _validate_parsimony_shared_peptides(self):
        """
        Internal ProteinInferenceParameter method to validate the shared peptides parameter
        """
        # Check if its all, best, or none
        if self.shared_peptides in Inference.SHARED_PEPTIDE_TYPES:
            self.logger.info("Using Shared Peptide types '{}'".format(self.shared_peptides))
        else:
            if self.shared_peptides.lower() == "none" or not self.shared_peptides:
                self.shared_peptides = None
                self.logger.info("Setting Shared Peptide type to None")
            else:
                raise ValueError(
                    "Shared Peptide types '{}' not supported, please use one of the following Shared Peptide types: '{}'".format(
                        self.shared_peptides, ", ".join(Inference.SHARED_PEPTIDE_TYPES)
                    )
                )

    def _validate_identifiers(self):
        """
        Internal ProteinInferenceParameter method to validate the decoy symbol, isoform symbol, and reviewed identifier symbol

        """
        if type(self.decoy_symbol) == str:
            self.logger.info(
                "Decoy Symbol set to: '{}'".format(
                    self.decoy_symbol
                )
            )
        else:
            raise ValueError(
                "Decoy Symbol must be a string, provided value: {}".format(
                    self.decoy_symbol
                )
            )


        if type(self.isoform_symbol) == str:
            self.logger.info(
                "Isoform Symbol set to: '{}'".format(
                    self.isoform_symbol
                )
            )
            if self.isoform_symbol.lower() == "none" or not self.isoform_symbol:
                self.isoform_symbol = None
                self.logger.info("Isoform Symbol set to None")
        else:
            if self.isoform_symbol:
                self.isoform_symbol = None
                self.logger.info("Isoform Symbol set to None")
            raise ValueError(
                "Isoform Symbol must be a string, provided value: {}".format(
                    self.isoform_symbol
                )
            )

        if type(self.reviewed_identifier_symbol) == str:
            self.logger.info(
                "Reviewed Identifier Symbol set to: '{}'".format(
                    self.reviewed_identifier_symbol
                )
            )
            if self.reviewed_identifier_symbol.lower() == "none" or not self.reviewed_identifier_symbol:
                self.reviewed_identifier_symbol = None
                self.logger.info("Reviewed Identifier Symbol set to None")
        else:
            if not self.reviewed_identifier_symbol:
                self.reviewed_identifier_symbol = None
                self.logger.info("Reviewed Identifier Symbol set to None")
            raise ValueError(
                "Reviewed Identifier Symbol must be a string, provided value: {}".format(
                    self.reviewed_identifier_symbol
                )
            )


    def override_q_restrict(self, data_class):
        """
        ProteinInferenceParameter method to override restrict_q if the input data does not contain q values.

        Args:
            data_class (protein_inference.datastore.DataStore): Data class

        """
        data_has_q = data_class.input_has_q()
        if data_has_q:
            pass
        else:
            if self.restrict_q:
                self.logger.warning(
                    "No Q values found in the input data, overriding parameters to not filter on Q value"
                )
                self.restrict_q = None

    def override_pep_restrict(self, data_class):
        """
        ProteinInferenceParameter method to override restrict_pep if the input data does not contain pep values.

        Args:
            data_class (protein_inference.datastore.DataStore): Data class

        """
        data_has_pep = data_class.input_has_pep()
        if data_has_pep:
            pass
        else:
            if self.restrict_pep:
                self.logger.warning(
                    "No Pep values found in the input data, overriding parameters to not filter on Pep value"
                )
                self.restrict_pep = None

    def override_custom_restrict(self, data_class):
        """
        ProteinInferenceParameter method to override restrict_custom if the input data does not contain custom score values.

        Args:
            data_class (protein_inference.datastore.DataStore): Data class

        """
        data_has_custom = data_class.input_has_custom()
        if data_has_custom:
            pass
        else:
            if self.restrict_custom:
                self.logger.warning(
                    "No Custom values found in the input data, overriding parameters to not filter on Custom value"
                )
                self.restrict_custom = None


    def fix_parameters_from_datastore(self,data_class):
        """
        ProteinInferenceParameter method to override restriction values in the parameter file if those scores do not exist in the input files

        Args:
            data_class (protein_inference.datastore.DataStore): Data class

        """

        self.override_q_restrict(data_class=data_class)
        self.override_pep_restrict(data_class=data_class)
        self.override_custom_restrict(data_class=data_class)

    def _fix_none_parameters(self):
        """
        Internal ProteinInferenceParameter method to fix parameters that have been defined as None
        These get read in as strings with YAML reader and need to be converted to None type
        """

        self._fix_grouping_type()
        self._fix_glpk_path()
        self._fix_lp_solver()
        self._fix_shared_peptides()


    def _fix_grouping_type(self):
        """
        Internal ProteinInferenceParameter method to override grouping type for None value
        """
        if self.grouping_type in ["None", "none", None]:
            self.grouping_type = None

    def _fix_glpk_path(self):
        """
        Internal ProteinInferenceParameter method to override glpk_path for None value
        """
        if self.glpk_path in ["None", "none", None]:
            self.glpk_path = None

    def _fix_lp_solver(self):
        """
        Internal ProteinInferenceParameter method to override lp_solver for None value
        """
        if self.lp_solver in ["None", "none", None]:
            self.lp_solver = None

    def _fix_shared_peptides(self):
        """
        Internal ProteinInferenceParameter method to override shared_peptides for None value
        """
        if self.shared_peptides in ["None", "none", None]:
            self.shared_peptides = None
