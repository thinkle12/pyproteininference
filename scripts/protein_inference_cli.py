#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""

import protein_inference
import argparse


parser = argparse.ArgumentParser(description='Protein Inference')
parser.add_argument("-t","--target", dest="target", required=False,
                    help="Input target psm output from percolator. Can either input one file or a list of files", metavar="FILE", nargs='+', type=str)
parser.add_argument("-d","--decoy", dest="decoy", required=False,
                    help="Input decoy psm output from percolator. Can either input one file or a list of files", metavar="FILE", nargs='+', type=str)
parser.add_argument("-f","--combined_files", dest="combined_files", required=False,
                    help="Input combined psm output from percolator. This should contain Target and Decoy PSMS. Can either input one file or a list of files", metavar="FILE", nargs='+', type=str)
parser.add_argument("-o","--output", dest="dir_name", required=False,
                    help="protein_inference Result Directory to write to - Name of file will be determined by parameters selected and parameter tag", metavar="DIR")
parser.add_argument("-a","--target_directory", dest="target_directory", required=False,
                    help="Directory that contains either .txt or .tsv input target psm data. Make sure the directory ONLY contains result files", metavar="DIR")
parser.add_argument("-b","--decoy_directory", dest="decoy_directory", required=False,
                    help="Directory that contains either .txt or .tsv input decoy psm data. Make sure the directory ONLY contains result files", metavar="DIR")
parser.add_argument("-c","--combined_directory", dest="combined_directory", required=False,
                    help="Directory that contains either .txt or .tsv input data with targets/decoys combined. Make sure the directory ONLY contains result files", metavar="DIR")
parser.add_argument("-db","--database", dest="database", required=True,
                    help="Provide the fasta formatted database used in the MS search", metavar="FILE")
parser.add_argument("-y","--yaml_params", dest="yaml_params", required=True,
                    help="Provide a Protein Inference Yaml Parameter File", metavar="FILE")
args = parser.parse_args()


pipeline = protein_inference.pipeline.ProteinInferencePipeline(parameter_file=args.yaml_params,
                                                               database_file=args.database,
                                                               target_files=args.target,
                                                               decoy_files=args.decoy,
                                                               combined_files=args.combined_files,
                                                               target_directory=args.target_directory,
                                                               decoy_directory=args.decoy_directory,
                                                               combined_directory=args.combined_directory,
                                                               output_directory=args.dir_name)
pipeline.execute()