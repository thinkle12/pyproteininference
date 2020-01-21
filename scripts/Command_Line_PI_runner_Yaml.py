#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""

import protein_inference
import argparse


parser = argparse.ArgumentParser(description='Protein Inference')
parser.add_argument("-t","--target", dest="target", required=True,
                    help="Input target psm output from percolator", metavar="FILE")
parser.add_argument("-d","--decoy", dest="decoy", required=True,
                    help="Input decoy psm output from percolator", metavar="FILE")
parser.add_argument("-o","--output", dest="dir_name", required=False,
                    help="protein_inference Result Directory to write to - Name of file will be determined by parameters selected and searchID", metavar="FILE")
parser.add_argument("-db","--database", dest="database", required=True,
                    help="Provide the database used in the MS search", metavar="FILE")
parser.add_argument("-ym","--yaml_params", dest="yaml_params", required=False,
                    help="Provide a Protein Inference Yaml Parameter File... If none given, default parameters will be ran", metavar="FILE")
args = parser.parse_args()


pipeline = protein_inference.pipeline.ProteinInferencePipeline(parameter_file=args.yaml_params,
                                                               database_file=args.database,
                                                               target_files=args.target,
                                                               decoy_files=args.decoy,
                                                               files=None,
                                                               output_directory=args.dir_name)
pipeline.execute()