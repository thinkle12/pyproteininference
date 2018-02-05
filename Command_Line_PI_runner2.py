#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:16:17 2017

@author: hinklet
"""


import ProteinInference
import argparse
from Digest import insilicodigest
import os


parser = argparse.ArgumentParser(description='Protein Inference')
parser.add_argument("-t","--target", dest="target", required=True,
                    help="Input target psm output from percolator", metavar="FILE")
parser.add_argument("-d","--decoy", dest="decoy", required=True,
                    help="Input decoy psm output from percolator", metavar="FILE")
parser.add_argument("-o","--output", dest="write_file", required=True,
                    help="ProteinInference Result File to Write", metavar="FILE")
parser.add_argument("-db","--database", dest="database", required=True,
                    help="Provide the database used in the MS search", metavar="FILE")
parser.add_argument("-mc","--missed_cleavages", dest="missed_cleavages", required=True,
                    help="Input a number of missed cleavages allowed - should be the same as what was done in the search", metavar="INT")
parser.add_argument("-dt","--digestion_type", dest="digestion_type", required=True,
                    help="Input the digestion type - example: trypsin", metavar="STR")
parser.add_argument("-pl","--peptide_legth", dest="peptide_length", required=False,
                    help="Provide the length of peptide to restrict by... 7 is a good minimum", metavar="INT")
parser.add_argument("-pep","--posterior_error_probability", dest="posterior_error_probability", required=False,
                    help="Provide the posterior error probability to restrict by...", metavar="FLOAT")
parser.add_argument("-qv","--q_value", dest="q_value", required=False,
                    help="Provide the q value to restrict by...", metavar="FLOAT")
parser.add_argument("-sm","--score_method", dest="score_method", required=True,
                    help="Select one of the following 'best_peptide_per_protein','iterative_downweighted_log','multiplicative_log','downweighted_multiplicative_log','downweighted_version2','top_two_combined','geometric_mean'", metavar="TYPE")
parser.add_argument("-st","--score_type", dest="score_type", required=True,
                    help="Select either 'pep_value' or 'q_value'", metavar="TYPE")
parser.add_argument("-ppk","--picker", dest="picker", action="store_true")
parser.add_argument("-gp","--group", dest="group", required=True,
                    help="select 'simple_subsetting', 'glpk', or 'multi_subsetting'", metavar="TYPE")
parser.add_argument("-fdr","--fdrcalc", dest="fdrcalc", required=True,
                    help="Input a number for fdr calculation ie .01 for a 1 percent fdr", metavar="FLOAT")
parser.add_argument("-roc","--roc_curve_filename", dest="roc_curve", required=False,
                    help="Provide the a filename with .pdf extension for roc curve output", metavar="FILE")
parser.add_argument("-ex","--export_type", dest="export_type", required=False,
                    help="select 'all', 'leads', 'comma_sep', 'q_value': q_value output not yet supported", metavar="LIST")
args = parser.parse_args()


tag = args.target.split('_percolator_')[0].split('/')[-1]

#Initiate the reader...
#Input for now is a target percolator output and a decoy percolator output
pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file=args.target,
                                                          decoy_file=args.decoy)

#Execeute the reader instance, this loads the data into the reader class
pep_and_prot_data.execute()
#Next create a data store which is a class that stores all data for all steps of the PI process
#Each method and each class calls from this data class to gather information for analyses
data = ProteinInference.datastore.DataStore(pep_and_prot_data)

#Here restrict the data to having peptides with length 7 or greater
if args.posterior_error_probability:
    pep_restrict = float(args.posterior_error_probability)
else:
    pep_restrict = None
if args.q_value:
    q_restrict = float(args.q_value)
else:
    q_restrict = None

if args.peptide_length:
    pl_restrict = int(args.peptide_length)
else:
    pl_restrict = None


print 'restricting data'
restrict = ProteinInference.datastore.RestrictMainData(data,peptide_length=pl_restrict,posterior_error_prob_threshold=pep_restrict,q_value_threshold=q_restrict)
restrict.execute()

#Here generate the pre score data using 'PEP' or 'Q' values
if args.score_type=='pep_value':
    score_setup = ProteinInference.datastore.PreScorePepValue(data)
if args.score_type=='q_value':
    score_setup = ProteinInference.datastore.PreScoreQValue(data)

#Execute score setup...
score_setup.execute()

#Here select scoring
if args.score_method=='best_peptide_per_protein':
    score = ProteinInference.scoring.BestPeptidePerProtein(data_class=data)
if args.score_method=='iterative_downweighted_log':
    score = ProteinInference.scoring.IterativeDownweightedLog(data_class=data)
if args.score_method=='multiplicative_log':
    score = ProteinInference.scoring.MultiplicativeLog(data_class=data)
if args.score_method=='downweighted_multiplicative_log':
    score = ProteinInference.scoring.DownweightedMultiplicativeLog(data_class=data)
if args.score_method=='downweighted_version2':
    score = ProteinInference.scoring.DownweightedVersion2(data_class=data)
if args.score_method=='top_two_combined':
    score = ProteinInference.scoring.TopTwoCombined(data_class=data)
if args.score_method=='geometric_mean':
    score = ProteinInference.scoring.GeometricMeanLog(data_class=data)

#Execute scoring...
score.execute()

#Run protein picker on the data
if args.picker:
    picker = ProteinInference.picker.StandardPicker(data_class=data)
    picker.execute()
else:
    pass

#Do in silico digest....
digest = insilicodigest.InSilicoDigest(database_path=args.database, num_miss_cleavs=int(args.missed_cleavages), digest_type=args.digestion_type)
digest.execute()

try:
    os.mkdir('glpkinout/')
except OSError:
    pass

#Run simple group subsetting
if args.group=='simple_subsetting':
    group = ProteinInference.grouping.SimpleSubsetting(data_class=data)
    group.execute()

#Run GLPK setup, runner, grouper...
if args.group=='glpk':
    glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data,glpkin_filename='glpkinout/glpkin_'+tag+'.mod')
    glpksetup.execute()
    glpkrun = ProteinInference.grouping.GlpkRunner(path_to_glpsol = '/gne/research/apps/protchem/glpk/bin/glpsol',glpkin='glpkinout/glpkin_'+tag+'.mod',glpkout='glpkinout/glpkout_'+tag+'.sol',file_override=False)
    glpkrun.execute()
    group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft', glpksolution_filename='glpkinout/glpkout_'+tag+'.sol')
    group.execute()

if args.group=='multi_subsetting':
    group = ProteinInference.grouping.MultiSubsetting(data_class=data)
    group.execute()


#Next run fdrcalc on the data....
fdr = ProteinInference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=float(args.fdrcalc))
fdr.execute()
#Finally we have our output restricted data...
restricted = data.fdr_restricted_grouped_scored_proteins
#Print the len of restricted data... which is how many protein groups pass FDR threshold
print 'Number of Proteins passing an FDR of'+str(args.fdrcalc)+' = '+str(len(restricted))



#Write the output to a csv...
if 'leads' in args.export_type:
    export = ProteinInference.export.CsvOutLeads(data_class=data,filename_out='leads_'+args.write_file)
if 'all' in args.export_type:
    export = ProteinInference.export.CsvOutAll(data_class=data,filename_out='all_'+args.write_file)
if 'comma_sep' in args.export_type:
    export = ProteinInference.export.CsvOutCommaSep(data_class=data, filename_out='comma_sep_'+args.write_file)
if 'q_value' in args.export_type:
    export = ProteinInference.export.CsvOutLeadsQValues(data_class=data,filename_out=args.write_file)
export.execute()

# if args.roc_curve:
#     roc = ProteinInference.benchmark.RocPlot(data_class=data)
#     roc.execute(pdf=args.roc_curve)


