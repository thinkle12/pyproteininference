import yaml
from logging import getLogger

from protein_inference.in_silico_digest import InSilicoDigest
from protein_inference.export import Export
from protein_inference.scoring import Score
from protein_inference.inference import Inference


class ProteinInferenceParameter(object):
    """
    Class that handles Percolator Parameters

    """

    def __init__(self, yaml_param_filepath, validate=True):
        """ Class to store percolator parameter information as an object

        Args:
            yaml_param_filepath (str): path to properly formatted parameter file specific to percoaltor

        Returns:
            None

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
        self.score_method = None
        self.score_type = None
        self.decoy_symbol = None
        self.isoform_symbol = None
        self.reviewed_identifier_symbol = None
        self.inference_type = None
        self.tag = None
        self.score = None
        self.max_identifiers_peptide_centric = None
        self.logger = getLogger('protein_inference.parameters.ProteinInferenceParameter.validate_parameters')

        self.convert_to_object()

        if validate:
            self.validate_parameters()



    def convert_to_object(self):
        """Function that takes a Percolator parameter file and converts it into a PercolatorParameter object

        Args:
            None

        Returns:
            None

        """
        if self.yaml_param_filepath:
            with open(self.yaml_param_filepath, 'r') as stream:
                yaml_params = yaml.load(stream)

            self.digest_type = yaml_params['parameters']['digest']['digest_type']
            self.export = yaml_params['parameters']['general']['export']
            self.fdr = yaml_params['parameters']['general']['fdr']
            self.glpk_path = yaml_params['parameters']['parsimony']['glpk_path']
            self.missed_cleavages = yaml_params['parameters']['digest']['missed_cleavages']
            self.picker = yaml_params['parameters']['general']['picker']
            self.restrict_pep = yaml_params['parameters']['data_restriction']['pep_restriction']
            self.restrict_peptide_length = yaml_params['parameters']['data_restriction']['peptide_length_restriction']
            self.restrict_q = yaml_params['parameters']['data_restriction']['q_value_restriction']
            self.restrict_custom = yaml_params['parameters']['data_restriction']['custom_restriction']
            self.score_method = yaml_params['parameters']['score']['score_method']
            self.score_type = yaml_params['parameters']['score']['score_type']
            self.decoy_symbol = yaml_params['parameters']['identifiers']['decoy_symbol']
            self.isoform_symbol = yaml_params['parameters']['identifiers']['isoform_symbol']
            self.reviewed_identifier_symbol = yaml_params['parameters']['identifiers']['reviewed_identifier_symbol']
            self.inference_type = yaml_params['parameters']['inference']['inference_type']
            self.tag = yaml_params['parameters']['general']["tag"]
            self.score = yaml_params["parameters"]['score']["score"]
            self.grouping_type = yaml_params["parameters"]['inference']["grouping_type"]
            self.max_identifiers_peptide_centric = yaml_params["parameters"]['peptide_centric']["max_identifiers"]
            self.lp_solver = yaml_params["parameters"]['parsimony']["lp_solver"]

        else:
            self.logger.info("Yaml parameter file not found, parameters set to default")
            self.digest_type = "trypsin"
            self.export = "q_value"
            self.fdr = 0.01
            self.glpk_path = "glpsol"
            self.missed_cleavages = 3
            self.picker = True
            self.restrict_pep = 0.9
            self.restrict_peptide_length = 7
            self.restrict_q = 0.005
            self.score_method = "multiplicative_log"
            self.score = "posterior_error_prob"
            self.decoy_symbol = "##"
            self.isoform_symbol = "-"
            self.reviewed_identifier_symbol = "sp|"
            self.inference_type = "peptide_centric"
            self.tag = "Test"
            self.score_type = "multiplicative"
            self.grouping_type = "shared_peptides"
            self.max_identifiers_peptide_centric = 5
            self.lp_solver = "pulp"

    def validate_parameters(self):
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


    def _validate_digest_type(self):
        # Make sure we have a valid digest type
        if self.digest_type in InSilicoDigest.LIST_OF_DIGEST_TYPES:
            self.logger.info("Using digest type '{}'".format(self.digest_type))
        else:
            raise ValueError("Digest Type '{}' not supported, please use one of the following enyzme digestions: '{}'".format(self.digest_type, ", ".join(InSilicoDigest.LIST_OF_DIGEST_TYPES)))

    def _validate_export_type(self):
        # Make sure we have a valid export type
        if self.export in Export.EXPORT_TYPES:
            self.logger.info("Using Export type '{}'".format(self.export))
        else:
            raise ValueError("Export Type '{}' not supported, please use one of the following export types: '{}'".format(self.export, ", ".join(Export.EXPORT_TYPES)))
        pass

    def _validate_floats(self):
        # Validate that FDR, cleavages, and restrict values are all floats and or ints if they need to be

        try:
            if 0<=float(self.fdr)<=1:
                self.logger.info("FDR Input {}".format(self.fdr))
            else:
                raise ValueError("FDR must be a decimal between 0 and 1, FDR provided: {}".format(self.fdr))
        except ValueError:
            raise ValueError("FDR must be a decimal between 0 and 1, FDR provided: {}".format(self.fdr))

        try:
            if 0<=float(self.restrict_pep)<=1:
                self.logger.info("PEP restriction {}".format(self.restrict_pep))
            else:
                raise ValueError("PEP restriction must be a decimal between 0 and 1, PEP restriction provided: {}".format(self.restrict_pep))
        except ValueError:
            if not self.restrict_pep or self.restrict_pep=="None":
                self.restrict_pep=None
            else:
                raise ValueError("PEP restriction must be a decimal between 0 and 1, PEP restriction provided: {}".format(
                    self.restrict_pep))

        try:
            if 0<=float(self.restrict_q)<=1:
                self.logger.info("Q Value restriction {}".format(self.restrict_q))

            else:
                raise ValueError("Q Value restriction must be a decimal between 0 and 1, Q Value restriction provided: {}".format(self.restrict_q))
        except ValueError:
            if not self.restrict_q or self.restrict_q=="None":
                self.restrict_q=None
            else:
                raise ValueError(
                    "Q Value restriction must be a decimal between 0 and 1, Q Value restriction provided: {}".format(
                        self.restrict_q))

        try:
            int(self.missed_cleavages)
        except ValueError:
            raise ValueError("Missed Cleavages must be an integer, Provided Missed Cleavages value: {}".format(
                self.missed_cleavages))
        if type(self.missed_cleavages)==int:
            self.logger.info("Missed Cleavages selected: {}".format(self.missed_cleavages))

        else:
            raise ValueError("Missed Cleavages must be an integer, Provided Missed Cleavages value: {}".format(self.missed_cleavages))

        if self.restrict_peptide_length:
            try:
                int(self.restrict_peptide_length)
            except ValueError:
                raise ValueError("Missed Cleavages must be an integer, Provided Missed Cleavages value: {}".format(
                    self.restrict_peptide_length))
            if type(self.restrict_peptide_length)==int:
                self.logger.info("Peptide Length Restriction: Len {}".format(self.restrict_peptide_length))
            else:
                raise ValueError("Missed Cleavages must be an integer, Provided Missed Cleavages value: {}".format(self.restrict_peptide_length))
        else:
            self.logger.info("Not Restricting by Peptide Length")

    def _validate_bools(self):
        # Make sure picker is a bool
        if type(self.picker)==bool:
            if self.picker:
                self.logger.info("Parameters loaded to run Picker")
            else:
                self.logger.info("Parameters loaded to NOT run Picker")
        else:
            raise ValueError("Picker Variable must be set to True or False, Picker Variable provided: {}".format(self.picker))


    def _validate_score_method(self):
        # Make sure we have the score method defined in code to use...
        if self.score_method in Score.SCORE_METHODS:
            self.logger.info("Using Score Method '{}'".format(self.score_method))
        else:
            raise ValueError("Score Method '{}' not supported, "
                             "please use one of the following Score Methods: '{}'".format(self.score_method, ", ".join(Score.SCORE_METHODS)))
        pass


    def _validate_score_type(self):
        # Make sure score type is multiplicative or additive
        if self.score_type in Score.SCORE_TYPES:
            self.logger.info("Using Score Type '{}'".format(self.score_type))
        else:
            raise ValueError("Score Type '{}' not supported, "
                             "please use one of the following Score Types: '{}'".format(self.score_type,
                                                                                      ", ".join(Score.SCORE_TYPES)))
        pass

    def _validate_score_combination(self):
        # Check to see if combination of score (column), method(multiplicative log, additive), and score type (multiplicative/additive) is possible...
        # This will be super custom

        if self.score_type=="additive" and self.score_method!= "additive":
            raise ValueError("If Score type is 'additive' (Higher PSM score is better) then you must use the 'additive' score method")

        elif self.score_type=="multiplicative" and self.score_method=="additive":
            raise ValueError("If Score type is 'multiplicative' (Lower PSM score is better) "
                             "then you must NOT use the 'additive' score method please "
                             "select one of the following score methods: {}".format(", ".join([x for x in Score.SCORE_METHODS if x!="additive"])))

        else:
            self.logger.info("Combination of Score Type: '{}' and Score Method: '{}' is Ok".format(self.score_type, self.score_method))

    def _validate_inference_type(self):
        # Check if its parsimony, exclusion, inclusion, none
        if self.inference_type in Inference.INFERENCE_TYPES:
            self.logger.info("Using inference type '{}'".format(self.inference_type))
        else:
            raise ValueError(
                "Inferece Type '{}' not supported, please use one of the following Inferece Types: '{}'".format(
                    self.inference_type, ", ".join(Inference.INFERENCE_TYPES)))

    def _validate_grouping_type(self):
        # Check if its parsimony, exclusion, inclusion, none
        if self.grouping_type in Inference.GROUPING_TYPES:
            self.logger.info("Using Grouping type '{}'".format(self.grouping_type))
        else:
            raise ValueError(
                "Grouping Type '{}' not supported, please use one of the following Grouping Types: '{}'".format(
                    self.grouping_type, ", ".join(Inference.GROUPING_TYPES)))

    def _validate_max_id(self):
        # Check if max_identifiers_peptide_centric param is an INT
        if type(self.max_identifiers_peptide_centric)==int:
            self.logger.info("Max Number of Indentifiers for Peptide Centric Inference: '{}'".format(self.max_identifiers_peptide_centric))
        else:
            raise ValueError("Max Number of Indentifiers for Peptide Centric Inference must be an integer, provided value: {}".format(self.max_identifiers_peptide_centric))

    def _validate_lp_solver(self):
        # Check if its pulp or glpk
        if self.lp_solver in Inference.LP_SOLVERS:
            self.logger.info("Using LP Solver '{}'".format(self.lp_solver))
        else:
            raise ValueError(
                "LP Solver '{}' not supported, please use one of the following LP Solvers: '{}'".format(
                    self.lp_solver, ", ".join(Inference.LP_SOLVERS)))

    def _override_inference_none(self):
        if self.inference_type in ["None", "none", None]:
            self.inference_type="none"