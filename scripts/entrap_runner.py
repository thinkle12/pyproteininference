#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:04:08 2017

@author: hinklet
"""
# import matplotlib
#
# matplotlib.use('Agg')
import os
from matplotlib.backends.backend_pdf import PdfPages
import numpy

qvalue_restriction = [.2,.1,.05]
pepvalue_restriction = [.9,.8,]

list_of_searchids = ['AB']

list_of_databases = ['complete_and_reversed_AB.fasta']
list_of_true_db = ['prest_pool_ab.fasta']

# other_db = ['all_benchmark_entrapment_data/entrapment_data/prest_pool_b.fasta']
other_db = [None]


pick_nopick = [True, False]

entrap_list = [None]
entrap_list2 = [None]

scoring_methods = ['ml','bppp','ttc','dwml','idwl','dw2','gm']
# scoring_methods = ['bppp']
import ProteinInference
import time

restricted_sets = []

large_list_of_results = [['scoring_method','score_type','picked','full_vert_dist_mean','restricted_vert_dist_mean','full_perp_dist_mean','restricted_perp_dist_mean','q_value_restriction','pep_value_restriction','number_of_proteins_passing','pearson']]
with PdfPages('plots/'+list_of_searchids[0]+'_entrap_plots_pep.pdf') as pdf:
    for k in range(len(scoring_methods)):
        for i in range(len(list_of_searchids)):
            for pnp in pick_nopick:
                for qvr in qvalue_restriction:
                    for pvr in pepvalue_restriction:
                        start_time = time.time()
                        #Initiate the reader...
                        #Input for now is a target percolator output and a decoy percolator output
                        pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file='all_benchmark_entrapment_data/combined_perc_output_files/'+str(list_of_searchids[i])+'_percolator_target_psm.txt',
                                                                                  decoy_file='all_benchmark_entrapment_data/combined_perc_output_files/'+str(list_of_searchids[i])+'_percolator_decoy_psm.txt')
                        # pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file='data/134285_k562_target_psm.txt',
                        #                                                           decoy_file='data/134285_k562_decoy_psm.txt')

                        #Execeute the reader instance, this loads the data into the reader class
                        pep_and_prot_data.execute()

                        #Next create a data store which is a class that stores all data for all steps of the PI process
                        #Each method and each class calls from this data class to gather information for analyses
                        data = ProteinInference.datastore.DataStore(pep_and_prot_data)

                        #Here restrict the data to having peptides with length 7 or greater and a pep of less than .9
                        restrict = ProteinInference.datastore.RestrictMainData(data, peptide_length=7, posterior_error_prob_threshold=pvr,q_value_threshold=qvr)
                        restrict.execute()

                        #Here generate the pre score data using 'PEP' values
                        score_setup = ProteinInference.datastore.PreScorePepValue(data)
                        score_setup.execute()

                        #Here we do scoring
                        if scoring_methods[k]=='bppp':
                            score = ProteinInference.scoring.BestPeptidePerProtein(data_class=data)
                        if scoring_methods[k] == 'ttc':
                            score = ProteinInference.scoring.TopTwoCombined(data_class=data)
                        if scoring_methods[k] == 'ml':
                            score = ProteinInference.scoring.MultiplicativeLog(data_class=data)
                        if scoring_methods[k] == 'dwml':
                            score = ProteinInference.scoring.DownweightedMultiplicativeLog(data_class=data)
                        if scoring_methods[k] == 'dw2':
                            score = ProteinInference.scoring.DownweightedVersion2(data_class=data)
                        if scoring_methods[k] == 'fm':
                            score = ProteinInference.scoring.FishersMethod(data_class=data)
                        if scoring_methods[k] == 'gm':
                            score = ProteinInference.scoring.GeometricMeanLog(data_class=data)
                        if scoring_methods[k] == 'idwl':
                            score = ProteinInference.scoring.IterativeDownweightedLog(data_class=data)
                        score.execute()
                        #This variable becomes the scored proteins as a real variable
                        scored_prots = data.scored_proteins
                        #Run protein picker on the data
                        if pnp:
                            picker = ProteinInference.picker.StandardPicker(data_class=data)
                            picker.execute()


                        #Run simple group subsetting
                        #group = ProteinInference.grouping.Simple_Subsetting(data_class=data)
                        #group.execute()

                        #Do in silico trypsin digestion
                        from Digest import insilicodigest
                        digest = insilicodigest.InSilicoDigest(database_path='all_benchmark_entrapment_data/entrapment_data/'+list_of_databases[i], num_miss_cleavs=2, digest_type='trypsin')
                        digest.execute()
                        #
                        #Run GLPK to generate the minimal list of proteins that account for the peptides
                        #Running GLPK consists of 3 classes, setup, runner, and grouper which need to be run in succession
                        glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data,glpkin_filename='glpkinout/glpkout_'+list_of_searchids[i]+'.mod')
                        glpksetup.execute()
                        runglpk = ProteinInference.grouping.GlpkRunner(path_to_glpsol = 'glpsol',glpkin = 'glpkinout/glpkout_'+list_of_searchids[i]+'.mod',glpkout = 'glpkinout/glpkout_'+list_of_searchids[i]+'.sol',file_override = False)
                        runglpk.execute()
                        group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft', glpksolution_filename='glpkinout/glpkout_'+list_of_searchids[i]+'.sol')
                        group.execute()

                        #Next run fdrcalc on the data....
                        fdr = ProteinInference.fdrcalc.EntrapFdr(data_class=data,entrapment_database='all_benchmark_entrapment_data/entrapment_data/prest_1000_random.fasta',other_database=other_db[0],false_discovery_rate=.05)
                        fdr.execute()



                        #Finally we have our output restricted data...
                        # restricted = data.fdr_restricted_grouped_scored_proteins
                        #
                        # restricted_sets.append(len(restricted))
                        #print it as well as the len...
                        #print restricted

                        #Write the output to a csv...
                        #output = ProteinInference.export.CsvOutAll(data_class=data, filename_out='output/hard_sp_override.csv')
                        #output.execute()

                        # output_leads = ProteinInference.export.CsvOutLeads(data_class=data, filename_out='output/protein_lead_export_'+list_of_searchids[i]+scoring_methods[i]+'.csv')
                        # output_leads.execute()

                        #qval_out_leads = ProteinInference.export.CsvOutLeadsQValues(data_class=data, filename_out='output/soft_override_leads_qvalue_restrict_pep09.csv')
                        #qval_out_leads.execute()

                        # roc = ProteinInference.benchmark.RocPlot(data_class=data)
                        # roc.execute(pdf='output/k562_134285_dwml.pdf')

                        entrap = ProteinInference.entrapment.GeneratePlot(data_class = data, entrapment_db= 'all_benchmark_entrapment_data/entrapment_data/prest_1000_random.fasta', true_db='all_benchmark_entrapment_data/entrapment_data/'+list_of_true_db[i], search_id=list_of_searchids[i], pdf=pdf, picked=pnp, qvr=qvr, pvr=pvr, other_database=other_db[0])
                        entrap.execute()


                        #youden_data = data.max_youdens_data
                        large_list_of_results.append([scoring_methods[k]])
                        large_list_of_results[-1].append(data.score_type)
                        if pnp:
                            large_list_of_results[-1].append('picked')
                        else:
                            large_list_of_results[-1].append('not_picked')
                        large_list_of_results[-1].append(entrap.full_vertical_distance_mean)
                        large_list_of_results[-1].append(entrap.vertical_distance_mean)
                        large_list_of_results[-1].append(entrap.full_perp_dist_mean)
                        large_list_of_results[-1].append(entrap.perp_dist_mean)
                        large_list_of_results[-1].append(qvr)
                        large_list_of_results[-1].append(pvr)
                        large_list_of_results[-1].append(str(len(entrap.prots_that_pass)))
                        large_list_of_results[-1].append(str(entrap.pr))


        print 'score type = '+str(data.score_type)
        print 'score method = '+str(data.score_method)
        end_time = time.time()-start_time
        print str(end_time)+' seconds'
        print str(end_time/60)+' minutes'
        print str(end_time/(60*60))+' hours'
        # print restricted_sets

nplist = numpy.array(large_list_of_results)

import csv

with open('AB_entrap_results_pep.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(large_list_of_results)