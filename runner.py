#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:04:08 2017

@author: hinklet
"""


import ProteinInference
import time
from matplotlib.backends.backend_pdf import PdfPages
import numpy

qvalue_restriction = [.1]
pepvalue_restriction = [.9]

list_of_searchids = ['134316']



pick_nopick = [True]


scoring_methods = ['dwml']

large_list_of_results = [['SearchID','Scoring_Methods','Scoring_Type','Picker','Qvr','Pvr','Proteins_Passing_1_Percent_FDR']]

with PdfPages('plots/'+list_of_searchids[0]+'_plot_ml.pdf') as pdf:
    for k in range(len(scoring_methods)):
        for i in range(len(list_of_searchids)):
            for pnp in pick_nopick:
                for qvr in qvalue_restriction:
                    for pvr in pepvalue_restriction:
                        start_time = time.time()
                        #Initiate the reader...
                        #Input for now is a target percolator output and a decoy percolator output
                        pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file='data/'+list_of_searchids[i]+'_percolator_target_psm.txt',
                                                                                  decoy_file='data/'+list_of_searchids[i]+'_percolator_decoy_psm.txt')

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

                        score = ProteinInference.scoring.DownweightedMultiplicativeLog(data_class=data)
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

                        digest = insilicodigest.InSilicoDigest(database_path='data/UniprotKBConcat1708_HUMAN.fasta',
                                                               num_miss_cleavs=2,
                                                               digest_type='trypsin')
                        digest.execute()

                        #
                        #Run GLPK to generate the minimal list of proteins that account for the peptides
                        #Running GLPK consists of 3 classes, setup, runner, and grouper which need to be run in succession
                        glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data,glpkin_filename='glpkinout/glpkout09_k562_134316.mod')
                        glpksetup.execute()
                        runglpk = ProteinInference.grouping.GlpkRunner(path_to_glpsol = 'glpsol',glpkin = 'glpkinout/glpkout09_k562_134316.mod',glpkout = 'glpkinout/glpkout09_k562_134316.sol',file_override = False)
                        runglpk.execute()
                        group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest, swissprot_override='soft', glpksolution_filename='glpkinout/glpkout09_k562_134316.sol')
                        group.execute()

                        #Next run fdrcalc on the data....
                        fdr = ProteinInference.fdrcalc.SetBasedFdr(data_class=data,false_discovery_rate=.01)
                        fdr.execute()



                        #Finally we have our output restricted data...
                        restricted = data.fdr_restricted_grouped_scored_proteins
                        #print it as well as the len...
                        #print restricted

                        #Write the output to a csv...
                        #output = ProteinInference.export.CsvOutAll(data_class=data, filename_out='output/hard_sp_override.csv')
                        #output.execute()

                        output_leads = ProteinInference.export.CsvOutLeads(data_class=data, filename_out='output/protein_lead_export_ml_k562_134316.csv')
                        output_leads.execute()

                        # qval_out_leads = ProteinInference.export.CsvOutLeadsQValues(data_class=data, filename_out='output/soft_override_leads_qvalue_restrict_pep09.csv')
                        # qval_out_leads.execute()

                        roc = ProteinInference.benchmark.RocPlot(data_class=data)
                        roc.execute(pdf=pdf)

                        # entrap = ProteinInference.entrapment.GeneratePlot(data_class = data, entrapment_db= 'entrapment_data/prest_1000_random.fasta', true_db='entrapment_data/'+list_of_true_db[i], search_id=list_of_searchids[i])
                        # entrap.execute()

                        #youden_data = data.max_youdens_data

                        large_list_of_results.append([list_of_searchids[i]])
                        large_list_of_results[-1].append(scoring_methods[k])
                        large_list_of_results[-1].append('PepValue')
                        if pnp:
                            large_list_of_results[-1].append('Yes Picker')
                        else:
                            large_list_of_results[-1].append('No Picker')
                        large_list_of_results[-1].append(str(qvr))
                        large_list_of_results[-1].append(str(pvr))
                        large_list_of_results[-1].append(str(len(restricted)))


                        print 'score type = '+str(data.score_type)
                        print 'score method = '+str(data.score_method)
                        end_time = time.time()-start_time
                        print str(end_time)+' seconds'
                        print str(end_time/60)+' minutes'
                        print str(end_time/(60*60))+' hours'

# import csv
#
# with open('test.csv', "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(large_list_of_results)
