#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:04:08 2017

@author: hinklet
"""
import matplotlib
#
matplotlib.use('Agg')
import os
from matplotlib.backends.backend_pdf import PdfPages

# list_of_searchids = ['145773','145774','145775']
# list_of_searchids = ['147169','147170','147171','147172','147173','147174','147175','147176','147177']
# list_of_searchids = ['147175','147176','147177','145773','145774','145775']
# list_of_searchids = ['145773','145774','145775','145779','145780','145781','145782','145783','145784']
# list_of_databases = ['complete_and_reversed_AB.fasta','complete_and_reversed_AB.fasta','complete_and_reversed_AB.fasta','complete_and_reversed_B.fasta','complete_and_reversed_B.fasta','complete_and_reversed_B.fasta','complete_and_reversed_A.fasta','complete_and_reversed_A.fasta','complete_and_reversed_A.fasta']
# list_of_databases = ['random_entrap_B_test.fasta','random_entrap_B_test.fasta','random_entrap_B_test.fasta','random_entrap_A_test.fasta','random_entrap_A_test.fasta','random_entrap_A_test.fasta','random_entrap_AB_test.fasta','random_entrap_AB_test.fasta','random_entrap_AB_test.fasta']
# list_of_true_db = ['prest_pool_b.fasta','prest_pool_b.fasta','prest_pool_b.fasta','prest_pool_a.fasta','prest_pool_a.fasta','prest_pool_a.fasta','prest_pool_ab.fasta','prest_pool_ab.fasta','prest_pool_ab.fasta']
# list_of_true_db = ['prest_pool_ab.fasta','prest_pool_ab.fasta','prest_pool_ab.fasta','prest_pool_b.fasta','prest_pool_b.fasta','prest_pool_b.fasta','prest_pool_a.fasta','prest_pool_a.fasta','prest_pool_a.fasta']
# scoring_methods = ['bppp','ttc','ml','dwml','dw2','fm','gm','mdwl']
# list_of_databases = os.listdir('complete_random_entrap_dbs')


runs_to_use = 'AB_rep1'


list_of_search_filenames = os.listdir('/gne/research/data/protchem/gfy/working/benchmark_comet_directories/perc_bsub/perc_output_files_new/')

list_of_search_filenames = [x for x in list_of_search_filenames if runs_to_use in x and 'target' in x ]
list_of_true_db = ['prest_pool_ab.fasta' for x in range(len(list_of_search_filenames))]

list_of_searchids = [x.split('/')[-1].split('_percolator_')[0] for x in list_of_search_filenames]

list_of_db_numbers = [x.split('.')[0].split('_')[3] for x in list_of_search_filenames]

list_of_databases = ['complete_random_entrap_dbs_new/complete_random_entrap_db_AB_'+x+'_new.fasta' for x in list_of_db_numbers]

entrap_list = [None]
entrap_list2 = [None]

scoring_methods = ['ml']
import ProteinInference
import time

with PdfPages('plots/'+scoring_methods[0] + '_' + runs_to_use + '.pdf') as pdf:
    for k in range(len(scoring_methods)):
        for i in range(len(list_of_search_filenames)):
            if list_of_search_filenames[i]!='AB_rep1_db_32_new_percolator_target_psm.txt':
                start_time = time.time()
                #Initiate the reader...
                #Input for now is a target percolator output and a decoy percolator output
                pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file='/gne/research/data/protchem/gfy/working/benchmark_comet_directories/perc_bsub/perc_output_files_new/'+list_of_search_filenames[i],
                                                                          decoy_file='/gne/research/data/protchem/gfy/working/benchmark_comet_directories/perc_bsub/perc_output_files_new/'+'_decoy_'.join(list_of_search_filenames[i].split('_target_')))
                # pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file='data/134285_k562_target_psm.txt',
                #                                                           decoy_file='data/134285_k562_decoy_psm.txt')

                #Execeute the reader instance, this loads the data into the reader class
                pep_and_prot_data.execute()

                #Next create a data store which is a class that stores all data for all steps of the PI process
                #Each method and each class calls from this data class to gather information for analyses
                data = ProteinInference.datastore.DataStore(pep_and_prot_data)

                #Here restrict the data to having peptides with length 7 or greater and a pep of less than .9
                restrict = ProteinInference.datastore.RestrictMainData(data, peptide_length=7, posterior_error_prob_threshold=.9)
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
                picker = ProteinInference.picker.StandardPicker(data_class=data)
                picker.execute()


                #Run simple group subsetting
                #group = ProteinInference.grouping.Simple_Subsetting(data_class=data)
                #group.execute()

                #Do in silico trypsin digestion
                from Digest import insilicodigest
                digest = insilicodigest.InSilicoDigest(database_path=list_of_databases[i], num_miss_cleavs=2, digest_type='trypsin')
                digest.execute()
                #
                #Run GLPK to generate the minimal list of proteins that account for the peptides
                #Running GLPK consists of 3 classes, setup, runner, and grouper which need to be run in succession
                glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data,glpkin_filename='glpkinout/glpkout_'+list_of_searchids[i]+'.mod')
                glpksetup.execute()
                runglpk = ProteinInference.grouping.GlpkRunner(path_to_glpsol = '/gne/research/apps/protchem/glpk/bin/glpsol',glpkin = 'glpkinout/glpkout_'+list_of_searchids[i]+'.mod',glpkout = 'glpkinout/glpkout_'+list_of_searchids[i]+'.sol',file_override = False)
                runglpk.execute()
                group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft', glpksolution_filename='glpkinout/glpkout_'+list_of_searchids[i]+'.sol')
                group.execute()

                # #Next run fdrcalc on the data....
                # fdr = ProteinInference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=.05)
                # fdr.execute()



                #Finally we have our output restricted data...
                # restricted = data.fdr_restricted_grouped_scored_proteins
                #print it as well as the len...
                #print restricted

                #Write the output to a csv...
                #output = ProteinInference.export.CsvOutAll(data_class=data, filename_out='output/hard_sp_override.csv')
                #output.execute()

                output_leads = ProteinInference.export.CsvOutLeads(data_class=data, filename_out='output/protein_lead_export_'+list_of_searchids[i]+'_'+scoring_methods[i]+'.csv')
                output_leads.execute()

                # qval_out_leads = ProteinInference.export.CsvOutLeadsQValues(data_class=data, filename_out='output/soft_override_leads_qvalue_restrict_pep09.csv')
                # qval_out_leads.execute()

                # roc = ProteinInference.benchmark.RocPlot(data_class=data)
                # roc.execute(pdf='output/k562_134285_dwml.pdf')

                entrap = ProteinInference.entrapment.GeneratePlot(data_class = data, entrapment_db= 'random_entrap_dbs_new/random_entrap_db_new_'+str(list_of_db_numbers[i])+'.fasta', true_db='entrapment_data/'+list_of_true_db[i], search_id=list_of_searchids[i],pdf=pdf)
                entrap.execute()


                #youden_data = data.max_youdens_data


    print 'score type = '+str(data.score_type)
    print 'score method = '+str(data.score_method)
    end_time = time.time()-start_time
    print str(end_time)+' seconds'
    print str(end_time/60)+' minutes'
    print str(end_time/(60*60))+' hours'
