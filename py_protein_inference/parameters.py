import yaml
from logging import getLogger

from py_protein_inference.in_silico_digest import InSilicoDigest
from py_protein_inference.export import Export
from py_protein_inference.scoring import Score
from py_protein_inference.inference import Inference


class ProteinInferenceParameter(object):
    """
    Class that handles data retrieval, storage, and validation of Protein Inference Parameters

    Attributes:
        yaml_param_filepath (str): path to properly formatted parameter file specific to Protein Inference
        digest_type (str): String that determines that type of in silico digestion for :py:class:`py_protein_inference.in_silico_digest.Digest`. Typically "trypsin"
        export (str): String to indicate the export type for :py:class:`py_protein_inference.export.Export`. Typically this is "psms", "peptides", or "psm_ids"
        fdr (float): Float to indicate FDR filtering
        glpk_path (str): Path to local installation of glpsol if inference_type="parsimony" and lp_solver="glpk"
        missed_cleavages (int): Integer to determine the number of missed cleavages in the database digestion :py:class:`py_protein_inference.in_silico_digest.Digest`
        picker (bool): True/False on whether or not to run the protein picker algorithm :py:meth:py_protein_inference.datastore.DataStore.protein_picker`
        restrict_pep (float/None): Float to restrict the posterior error probability values by in the PSM input. Used in :py:meth:py_protein_inference.datastore.DataStore.restrict_psm_data`
        restrict_peptide_length (int/None): Float to restrict the peptide length values by in the PSM input. Used in :py:meth:py_protein_inference.datastore.DataStore.restrict_psm_data`
        restrict_q (float/None): Float to restrict the q values by in the PSM input. Used in :py:meth:py_protein_inference.datastore.DataStore.restrict_psm_data`
        restrict_custom (float/None): Float to restrict the custom values by in the PSM input. Used in :py:meth:py_protein_inference.datastore.DataStore.restrict_psm_data`. Filtering depends on score_type variable. If score_type is multiplicative then values that are less than restrict_custom are kept. If score_type is additive then values that are more than restrict_custom are kept.
        protein_score (str): String to determine the way in which Proteins are scored can be any of the SCORE_METHODS in :py:class:`py_protein_inference.scoring.Score`
        psm_score_type (str): String to determine the type of score that the PSM scores are (Additive or Multiplicative) can be any of the SCORE_TYPES in :py:class:`py_protein_inference.scoring.Score`
        decoy_symbol (str): String to denote decoy proteins from target proteins. IE "##"
        isoform_symbol (str): String to denote isoforms from regular proteins. IE "-". Can also be None
        reviewed_identifier_symbol (str): String to denote a "Reviewed" Protein. Typically this is: "sp|" if using Uniprot Fasta database
        inference_type (str): String to determine the inference procedure. Can be any value of INFERENCE_TYPES of :py:class:`py_protein_inference.inference.Inference` object
        tag (str): String to be added to output files
        psm_score (str): String that indicates the PSM input score. The value should match the string in the input data of the score you want to use for PSM score. This score will be used in scoring methods here: :py:class:`py_protein_inference.scoring.Score`
        grouping_type (str/None): String to determine the grouping procedure. Can be any value of GROUPING_TYPES of :py:class:`py_protein_inference.inference.Inference` object
        max_identifiers_peptide_centric (int): Maximum number of identifiers to assign to a group when running peptide_centric inference. Typically this is 10 or 5.
        lp_solver (str/None): The LP solver to use if inference_type="Parsimony". Can be any value in LP_SOLVERS in the :py:class:`py_protein_inference.inference.Inference` object
        logger (logger.logging): Logger object

    """

    PARENT_PARAMETER_KEY = "parameters"

    GENERAL_PARAMETER_KEY = "general"
    DATA_RESTRICTION_PARAMETER_KEY = "data_restriction"
    SCORE_PARAMETER_KEY = "score"
    IDENTIFIERS_PARAMETER_KEY = "identifiers"
    INFERENCE_PARAMETER_KEY = "inference"
    DIGEST_PARAMETER_KEY = "digest"
    PARSIMONY_PARAMETER_KEY = "parsimony"
    PEPTIDE_CENTRIC_PARAMETER_KEY = "peptide_centric"

    PARAMETER_MAIN_KEYS = {
        GENERAL_PARAMETER_KEY,
        DATA_RESTRICTION_PARAMETER_KEY,
        SCORE_PARAMETER_KEY,
        IDENTIFIERS_PARAMETER_KEY,
        INFERENCE_PARAMETER_KEY,
        DIGEST_PARAMETER_KEY,
        PARSIMONY_PARAMETER_KEY,
        PEPTIDE_CENTRIC_PARAMETER_KEY,
    }

    EXPORT_PARAMETER = "export"
    FDR_PARAMETER = "fdr"
    PICKER_PARAMETER = "picker"
    TAG_PARAMETER = "tag"

    GENERAL_PARAMETER_SUB_KEYS = {EXPORT_PARAMETER, FDR_PARAMETER, PICKER_PARAMETER, TAG_PARAMETER}

    PEP_RESTRICT_PARAMETER = "pep_restriction"
    PEPTIDE_LENGTH_RESTRICT_PARAMETER = "peptide_length_restriction"
    Q_VALUE_RESTRICT_PARAMETER = "q_value_restriction"
    CUSTOM_RESTRICT_PARAMETER = "custom_restriction"

    DATA_RESTRICTION_PARAMETER_SUB_KEYS = {
        PEP_RESTRICT_PARAMETER,
        PEPTIDE_LENGTH_RESTRICT_PARAMETER,
        Q_VALUE_RESTRICT_PARAMETER,
        CUSTOM_RESTRICT_PARAMETER,
    }

    PROTEIN_SCORE_PARAMETER = "protein_score"
    PSM_SCORE_PARAMETER = "psm_score"
    PSM_SCORE_TYPE_PARAMETER = "psm_score_type"

    SCORE_PARAMETER_SUB_KEYS = {PROTEIN_SCORE_PARAMETER, PSM_SCORE_PARAMETER, PSM_SCORE_TYPE_PARAMETER}

    DECOY_SYMBOL_PARAMETER = "decoy_symbol"
    ISOFORM_SYMBOL_PARAMETER = "isoform_symbol"
    REVIEWED_IDENTIFIER_PARAMETER = "reviewed_identifier_symbol"

    IDENTIFIER_SUB_KEYS = {DECOY_SYMBOL_PARAMETER, ISOFORM_SYMBOL_PARAMETER, REVIEWED_IDENTIFIER_PARAMETER}

    INFERENCE_TYPE_PARAMETER = "inference_type"
    GROUPING_TYPE_PARAMETER = "grouping_type"

    INFERENCE_SUB_KEYS = {INFERENCE_TYPE_PARAMETER, GROUPING_TYPE_PARAMETER}

    DIGEST_TYPE_PARAMETER = "digest_type"
    MISSED_CLEAV_PARAMETER = "missed_cleavages"

    DIGEST_SUB_KEYS = {DIGEST_TYPE_PARAMETER, MISSED_CLEAV_PARAMETER}

    LP_SOLVER_PARAMETER = "lp_solver"
    GLPK_PATH_PARAMETER = "glpk_path"
    SHARED_PEPTIDES_PARAMETER = "shared_peptides"

    PARSIMONY_SUB_KEYS = {LP_SOLVER_PARAMETER, GLPK_PATH_PARAMETER, SHARED_PEPTIDES_PARAMETER}

    MAX_IDENTIFIERS_PARAMETER = "max_identifiers"

    PEPTIDE_CENTRIC_SUB_KEYS = {MAX_IDENTIFIERS_PARAMETER}

    def __init__(self, yaml_param_filepath, validate=True):
        """Class to store Protein Inference parameter information as an object

        Args:
            yaml_param_filepath (str): path to properly formatted parameter file specific to Protein Inference
            validate (bool): True/False on whether to validate the parameter file of interest

        Returns:
            None

        Example:
            >>> py_protein_inference.parameters.ProteinInferenceParameter(
            >>>     yaml_param_filepath = "/path/to/py_protein_inference_params.yaml", validate=True
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
        self.logger = getLogger("py_protein_inference.parameters.ProteinInferenceParameter.validate_parameters")
        self.validate = validate

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

            if self.validate:
                self._validate_parameter_shape(yaml_params=yaml_params)

            self.digest_type = yaml_params[self.PARENT_PARAMETER_KEY][self.DIGEST_PARAMETER_KEY][
                self.DIGEST_TYPE_PARAMETER
            ]
            self.export = yaml_params[self.PARENT_PARAMETER_KEY][self.GENERAL_PARAMETER_KEY][self.EXPORT_PARAMETER]
            self.fdr = yaml_params[self.PARENT_PARAMETER_KEY][self.GENERAL_PARAMETER_KEY][self.FDR_PARAMETER]
            self.glpk_path = yaml_params[self.PARENT_PARAMETER_KEY][self.PARSIMONY_PARAMETER_KEY][
                self.GLPK_PATH_PARAMETER
            ]
            self.missed_cleavages = yaml_params[self.PARENT_PARAMETER_KEY][self.DIGEST_PARAMETER_KEY][
                self.MISSED_CLEAV_PARAMETER
            ]
            self.picker = yaml_params[self.PARENT_PARAMETER_KEY][self.GENERAL_PARAMETER_KEY][self.PICKER_PARAMETER]
            self.restrict_pep = yaml_params[self.PARENT_PARAMETER_KEY][self.DATA_RESTRICTION_PARAMETER_KEY][
                self.PEP_RESTRICT_PARAMETER
            ]
            self.restrict_peptide_length = yaml_params[self.PARENT_PARAMETER_KEY][self.DATA_RESTRICTION_PARAMETER_KEY][
                self.PEPTIDE_LENGTH_RESTRICT_PARAMETER
            ]
            self.restrict_q = yaml_params[self.PARENT_PARAMETER_KEY][self.DATA_RESTRICTION_PARAMETER_KEY][
                self.Q_VALUE_RESTRICT_PARAMETER
            ]
            self.restrict_custom = yaml_params[self.PARENT_PARAMETER_KEY][self.DATA_RESTRICTION_PARAMETER_KEY][
                self.CUSTOM_RESTRICT_PARAMETER
            ]
            self.protein_score = yaml_params[self.PARENT_PARAMETER_KEY][self.SCORE_PARAMETER_KEY][
                self.PROTEIN_SCORE_PARAMETER
            ]
            self.psm_score_type = yaml_params[self.PARENT_PARAMETER_KEY][self.SCORE_PARAMETER_KEY][
                self.PSM_SCORE_TYPE_PARAMETER
            ]
            self.decoy_symbol = yaml_params[self.PARENT_PARAMETER_KEY][self.IDENTIFIERS_PARAMETER_KEY][
                self.DECOY_SYMBOL_PARAMETER
            ]
            self.isoform_symbol = yaml_params[self.PARENT_PARAMETER_KEY][self.IDENTIFIERS_PARAMETER_KEY][
                self.ISOFORM_SYMBOL_PARAMETER
            ]
            self.reviewed_identifier_symbol = yaml_params[self.PARENT_PARAMETER_KEY][self.IDENTIFIERS_PARAMETER_KEY][
                self.REVIEWED_IDENTIFIER_PARAMETER
            ]
            self.inference_type = yaml_params[self.PARENT_PARAMETER_KEY][self.INFERENCE_PARAMETER_KEY][
                self.INFERENCE_TYPE_PARAMETER
            ]
            self.tag = yaml_params[self.PARENT_PARAMETER_KEY][self.GENERAL_PARAMETER_KEY][self.TAG_PARAMETER]
            self.psm_score = yaml_params[self.PARENT_PARAMETER_KEY][self.SCORE_PARAMETER_KEY][self.PSM_SCORE_PARAMETER]
            self.grouping_type = yaml_params[self.PARENT_PARAMETER_KEY][self.INFERENCE_PARAMETER_KEY][
                self.GROUPING_TYPE_PARAMETER
            ]
            self.max_identifiers_peptide_centric = yaml_params[self.PARENT_PARAMETER_KEY][
                self.PEPTIDE_CENTRIC_PARAMETER_KEY
            ][self.MAX_IDENTIFIERS_PARAMETER]
            self.lp_solver = yaml_params[self.PARENT_PARAMETER_KEY][self.PARSIMONY_PARAMETER_KEY][
                self.LP_SOLVER_PARAMETER
            ]
            try:
                # Do try except here to make old param files backwards compatible
                self.shared_peptides = yaml_params[self.PARENT_PARAMETER_KEY][self.PARSIMONY_PARAMETER_KEY][
                    self.SHARED_PEPTIDES_PARAMETER
                ]
            except KeyError:
                self.shared_peptides = Inference.ALL_SHARED_PEPTIDES

        else:
            self.logger.warning("Yaml parameter file not found, parameters set to default")
            self.digest_type = "trypsin"
            self.export = "q_value"
            self.fdr = 0.01
            self.glpk_path = "glpsol"
            self.missed_cleavages = 3
            self.picker = True
            self.restrict_pep = 0.9
            self.restrict_peptide_length = 7
            self.restrict_q = 0.005
            self.restrict_custom = "None"
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
            self.shared_peptides = "all"

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
        self._validate_parsimony_shared_peptides()

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
            raise ValueError("FDR must be a decimal between 0 and 1, FDR provided: {}".format(self.fdr))

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
            self.logger.info("Missed Cleavages selected: {}".format(self.missed_cleavages))
        except ValueError:
            raise ValueError(
                "Missed Cleavages must be an integer, Provided Missed Cleavages value: {}".format(self.missed_cleavages)
            )

        try:
            int(self.restrict_peptide_length)
            self.logger.info("Peptide Length Restriction: Len {}".format(self.restrict_peptide_length))
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
        except ValueError or TypeError:
            if not self.restrict_custom or self.restrict_custom.lower() == "none":
                self.restrict_custom = None
                self.logger.info("Not Restricting by Custom Value")
            else:
                raise ValueError(
                    "Custom restriction must be a number, Custom restriction provided: {}".format(self.restrict_custom)
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
                "Picker Variable must be set to True or False, Picker Variable provided: {}".format(self.picker)
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
                        self.grouping_type, Inference.GROUPING_TYPES
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
                        self.shared_peptides, Inference.SHARED_PEPTIDE_TYPES
                    )
                )

    def _validate_identifiers(self):
        """
        Internal ProteinInferenceParameter method to validate the decoy symbol, isoform symbol, and reviewed identifier symbol

        """
        if type(self.decoy_symbol) == str:
            self.logger.info("Decoy Symbol set to: '{}'".format(self.decoy_symbol))
        else:
            raise ValueError("Decoy Symbol must be a string, provided value: {}".format(self.decoy_symbol))

        if type(self.isoform_symbol) == str:
            self.logger.info("Isoform Symbol set to: '{}'".format(self.isoform_symbol))
            if self.isoform_symbol.lower() == "none" or not self.isoform_symbol:
                self.isoform_symbol = None
                self.logger.info("Isoform Symbol set to None")
        else:
            if self.isoform_symbol:
                self.isoform_symbol = None
                self.logger.info("Isoform Symbol set to None")
            raise ValueError("Isoform Symbol must be a string, provided value: {}".format(self.isoform_symbol))

        if type(self.reviewed_identifier_symbol) == str:
            self.logger.info("Reviewed Identifier Symbol set to: '{}'".format(self.reviewed_identifier_symbol))
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

    def _validate_parameter_shape(self, yaml_params):
        if self.PARENT_PARAMETER_KEY in yaml_params.keys():
            self.logger.info("Main Parameter Key is Present")
        else:
            raise ValueError(
                "Key {} needs to be defined as the outermost parameter group".format(self.PARENT_PARAMETER_KEY)
            )

        if self.PARAMETER_MAIN_KEYS.issubset(yaml_params[self.PARENT_PARAMETER_KEY]):
            self.logger.info("All Sub Parameter Keys Present")
        else:
            raise ValueError(
                "All of the following values: {}. Need to be Sub Parameters in the Yaml Parameter file".format(
                    ", ".join(self.PARAMETER_MAIN_KEYS),
                )
            )

        try:
            general_params = yaml_params[self.PARENT_PARAMETER_KEY][self.GENERAL_PARAMETER_KEY]
            for gkey in self.GENERAL_PARAMETER_SUB_KEYS:
                if gkey in general_params.keys():
                    pass
                else:
                    raise ValueError(
                        "General Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the general parameter field".format(
                            gkey
                        )
                    )

        except KeyError:
            raise ValueError("'general' sub Parameter not defined in the parameter file")

        try:
            data_res_params = yaml_params[self.PARENT_PARAMETER_KEY][self.DATA_RESTRICTION_PARAMETER_KEY]
            for drkey in self.DATA_RESTRICTION_PARAMETER_SUB_KEYS:
                if drkey in data_res_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Data Restriction Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the data_restriction parameter field".format(
                            drkey
                        )
                    )

        except KeyError:
            raise ValueError("'data_restriction' sub Parameter not defined in the parameter file")

        try:
            score_params = yaml_params[self.PARENT_PARAMETER_KEY][self.SCORE_PARAMETER_KEY]
            for skey in self.SCORE_PARAMETER_SUB_KEYS:
                if skey in score_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Score Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the score parameter field".format(
                            skey
                        )
                    )

        except KeyError:
            raise ValueError("'score' sub Parameter not defined in the parameter file")

        try:
            id_params = yaml_params[self.PARENT_PARAMETER_KEY][self.IDENTIFIERS_PARAMETER_KEY]
            for ikey in self.IDENTIFIER_SUB_KEYS:
                if ikey in id_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Identifiers Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the identifiers parameter field".format(
                            ikey
                        )
                    )

        except KeyError:
            raise ValueError("'identifiers' sub Parameter not defined in the parameter file")

        try:
            inf_params = yaml_params[self.PARENT_PARAMETER_KEY][self.INFERENCE_PARAMETER_KEY]
            for infkey in self.INFERENCE_SUB_KEYS:
                if infkey in inf_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Inference Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the inference parameter field".format(
                            infkey
                        )
                    )

        except KeyError:
            raise ValueError("'inference' sub Parameter not defined in the parameter file")

        try:
            digest_params = yaml_params[self.PARENT_PARAMETER_KEY][self.DIGEST_PARAMETER_KEY]
            for dkey in self.DIGEST_SUB_KEYS:
                if dkey in digest_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Digest Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the digest parameter field".format(
                            dkey
                        )
                    )

        except KeyError:
            raise ValueError("'digest' sub Parameter not defined in the parameter file")

        try:
            parsimony_params = yaml_params[self.PARENT_PARAMETER_KEY][self.PARSIMONY_PARAMETER_KEY]
            for pkey in self.PARSIMONY_SUB_KEYS:
                if pkey in parsimony_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Parsimony Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the parsimony parameter field".format(
                            pkey
                        )
                    )

        except KeyError:
            raise ValueError("'parsimony' sub Parameter not defined in the parameter file")

        try:
            pep_cen_params = yaml_params[self.PARENT_PARAMETER_KEY][self.PEPTIDE_CENTRIC_PARAMETER_KEY]
            for pckey in self.PEPTIDE_CENTRIC_SUB_KEYS:
                if pckey in pep_cen_params.keys():
                    pass
                else:
                    raise ValueError(
                        "Peptide Centric Sub Parameter '{}' is not found in the parameter file. Please add it as a sub parameter of the peptide_centric parameter field".format(
                            pckey
                        )
                    )

        except KeyError:
            raise ValueError("'peptide_centric' sub Parameter not defined in the parameter file")

    def override_q_restrict(self, data):
        """
        ProteinInferenceParameter method to override restrict_q if the input data does not contain q values.

        Args:
            data (py_protein_inference.datastore.DataStore): Data Object

        """
        data_has_q = data.input_has_q()
        if data_has_q:
            pass
        else:
            if self.restrict_q:
                self.logger.warning(
                    "No Q values found in the input data, overriding parameters to not filter on Q value"
                )
                self.restrict_q = None

    def override_pep_restrict(self, data):
        """
        ProteinInferenceParameter method to override restrict_pep if the input data does not contain pep values.

        Args:
            data (py_protein_inference.datastore.DataStore): Data Object

        """
        data_has_pep = data.input_has_pep()
        if data_has_pep:
            pass
        else:
            if self.restrict_pep:
                self.logger.warning(
                    "No Pep values found in the input data, overriding parameters to not filter on Pep value"
                )
                self.restrict_pep = None

    def override_custom_restrict(self, data):
        """
        ProteinInferenceParameter method to override restrict_custom if the input data does not contain custom score values.

        Args:
            data (py_protein_inference.datastore.DataStore): Data Object

        """
        data_has_custom = data.input_has_custom()
        if data_has_custom:
            pass
        else:
            if self.restrict_custom:
                self.logger.warning(
                    "No Custom values found in the input data, overriding parameters to not filter on Custom value"
                )
                self.restrict_custom = None

    def fix_parameters_from_datastore(self, data):
        """
        ProteinInferenceParameter method to override restriction values in the parameter file if those scores do not exist in the input files

        Args:
            data (py_protein_inference.datastore.DataStore): Data Object

        """

        self.override_q_restrict(data=data)
        self.override_pep_restrict(data=data)
        self.override_custom_restrict(data=data)

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
