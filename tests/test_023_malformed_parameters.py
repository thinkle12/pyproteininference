#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""
import tempfile
from unittest import TestCase
from pkg_resources import resource_filename

from py_protein_inference.parameters import ProteinInferenceParameter
import yaml
import logging
import os

PARAMETER_FILE = resource_filename(
    "py_protein_inference", "../tests/data/test_params_inclusion.yaml"
)
OUTPUT_DIR = tempfile.gettempdir()

NEW_PARAMETER_FILE = os.path.join(OUTPUT_DIR, "malformed_protein_inference_params.yaml")

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("py_protein_inference.tests.test_023_malformed_parameters")


class TestMalformedParameters(TestCase):

    def test_malformed_parameters(self):

        # TODO maybe add a shape catcher... need certain headers and subheaders in params....

        # Read params...
        # Make sure everything is as it should be...
        protein_inference_parameters = ProteinInferenceParameter(
            yaml_param_filepath=PARAMETER_FILE
        )

        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["digest"]["digest_type"] = "Unimplemented Digest Type"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )

        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["general"]["export"] = None

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["general"]["fdr"] = "string"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["digest"]["missed_cleavages"] = "not an integer"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["general"]["picker"] = 456

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["data_restriction"]["pep_restriction"] = "string"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["data_restriction"]["peptide_length_restriction"] = "string"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["data_restriction"]["q_value_restriction"] = "string"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["data_restriction"]["custom_restriction"] = "string"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["score"]["protein_score"] = "invalid protein score"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["score"]["psm_score_type"] = "invalid psm score type"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["identifiers"]["decoy_symbol"] = None

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["inference"]["inference_type"] = "improper inference type"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["inference"]["grouping_type"] = "improper grouping type"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["peptide_centric"]["max_identifiers"] = "string"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["parsimony"]["lp_solver"] = "improper lp solver"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )


        # Read params with pyyaml
        with open(PARAMETER_FILE, "r") as stream:
            yaml_params = yaml.load(stream, Loader=yaml.Loader)

        # Edit the params to something incorrect...
        edited_yaml_params = yaml_params
        edited_yaml_params["parameters"]["parsimony"]["shared_peptides"] = "improper shared peptides"

        # Rewrite the params...
        with open(NEW_PARAMETER_FILE, 'w') as file:
             new_params = yaml.dump(edited_yaml_params, file)

        # Re-read the params with ProteinInferenceParameter...
        with self.assertRaises(ValueError):
            protein_inference_parameters = ProteinInferenceParameter(
                yaml_param_filepath=NEW_PARAMETER_FILE
            )

    def test_empty_parameters(self):
        protein_inference_parameters = ProteinInferenceParameter(
            yaml_param_filepath=None
        )

        self.assertEqual(protein_inference_parameters.digest_type, "trypsin")
        self.assertEqual(protein_inference_parameters.export, "q_value")
        self.assertEqual(protein_inference_parameters.fdr, 0.01)
        self.assertEqual(protein_inference_parameters.glpk_path, "glpsol")
        self.assertEqual(protein_inference_parameters.missed_cleavages, 3)
        self.assertEqual(protein_inference_parameters.picker, True)
        self.assertEqual(protein_inference_parameters.restrict_pep, 0.9)
        self.assertEqual(protein_inference_parameters.restrict_peptide_length, 7)
        self.assertEqual(protein_inference_parameters.restrict_q, 0.005)
        self.assertEqual(
            protein_inference_parameters.protein_score, "multiplicative_log"
        )
        self.assertEqual(protein_inference_parameters.psm_score, "posterior_error_prob")
        self.assertEqual(protein_inference_parameters.psm_score_type, "multiplicative")
        self.assertEqual(protein_inference_parameters.decoy_symbol, "##")
        self.assertEqual(protein_inference_parameters.isoform_symbol, "-")
        self.assertEqual(protein_inference_parameters.reviewed_identifier_symbol, "sp|")
        self.assertEqual(protein_inference_parameters.inference_type, "peptide_centric")
        self.assertEqual(protein_inference_parameters.tag, "example_tag")
        self.assertEqual(protein_inference_parameters.grouping_type, "shared_peptides")
        self.assertEqual(
            protein_inference_parameters.max_identifiers_peptide_centric, 5
        )
        self.assertEqual(protein_inference_parameters.lp_solver, "pulp")
        self.assertEqual(protein_inference_parameters.restrict_custom, None)
        self.assertEqual(protein_inference_parameters.shared_peptides, "all")




