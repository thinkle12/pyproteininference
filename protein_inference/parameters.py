import yaml

from proteologic.data.domain.mapper import AbstractPhysicalObject



class ProteinInferenceParameter(object):
    """
    Class that handles Percolator Parameters

    """

    def __init__(self, yaml_param_filepath):
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

        self.convert_to_object()



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

            self.digest_type = yaml_params['Parameters']['Digest_Type']
            self.export = yaml_params['Parameters']['Export']
            self.fdr = yaml_params['Parameters']['FDR']
            self.glpk_path = yaml_params['Parameters']['GLPK_Path']
            self.missed_cleavages = yaml_params['Parameters']['Missed_Cleavages']
            self.picker = yaml_params['Parameters']['Picker']
            self.restrict_pep = yaml_params['Parameters']['Restrict_Pep']
            self.restrict_peptide_length = yaml_params['Parameters']['Restrict_Peptide_Length']
            self.restrict_q = yaml_params['Parameters']['Restrict_Q']
            self.score_method = yaml_params['Parameters']['Score_Method']
            self.score_type = yaml_params['Parameters']['Score_Type']
            self.decoy_symbol = yaml_params['Parameters']['Decoy_Symbol']
            self.isoform_symbol = yaml_params['Parameters']['Isoform_Symbol']
            self.reviewed_identifier_symbol = yaml_params['Parameters']['Reviewed_Identifier_Symbol']
            self.inference_type = yaml_params['Parameters']['Inference_Type']
            self.tag = yaml_params['Parameters']["Tag"]

        else:
            print("Yaml parameter file not found, parameters set to default")
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
            self.score_type = "pep_value"
            self.decoy_symbol = "##"
            self.isoform_symbol = "-"
            self.reviewed_identifier_symbol = "sp|"
            self.inference_type = "parsimony"
            self.tag = "Test"

