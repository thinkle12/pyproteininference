#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:04:08 2017

@author: hinklet
"""
import matplotlib
from Digest import insilicodigest
import yaml
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

# Put a dictionary mapper in here...
# This will map the runs to use... to the db letter and the rep number...
# Make if more automated....

runs_to_use = 'A_rep1'

list_of_search_filenames = os.listdir('/gne/research/data/protchem/gfy/working/benchmark_comet_directories/perc_bsub/perc_output_files_new/')

list_of_search_filenames = [x for x in list_of_search_filenames if runs_to_use in x and 'target' in x ]
list_of_true_db = ['prest_pool_a.fasta' for x in range(len(list_of_search_filenames))]

list_of_searchids = [x.split('/')[-1].split('_percolator_')[0] for x in list_of_search_filenames]

list_of_db_numbers = [x.split('.')[0].split('_')[3] for x in list_of_search_filenames]

list_of_databases = ['complete_random_entrap_dbs_new/complete_random_entrap_db_A_'+x+'_new.fasta' for x in list_of_db_numbers]

entrap_list = [None]
entrap_list2 = [None]


dir_name = "/Users/hinklet/PI_output_benchmark/"
yaml_params = "/Users/hinklet/PythonPackages/ProteinInference/parameters/Protein_Inference_Params.yaml"
scoring_methods = ['idwl']
import ProteinInference
import time

with PdfPages('plots/'+scoring_methods[0] + '_' + runs_to_use + '.pdf') as pdf:
    for i in range(len(list_of_search_filenames)):
        target = '/gne/research/data/protchem/gfy/working/benchmark_comet_directories/perc_bsub/perc_output_files_new/'+list_of_search_filenames[i]
        decoy ='/gne/research/data/protchem/gfy/working/benchmark_comet_directories/perc_bsub/perc_output_files_new/'+'_decoy_'.join(list_of_search_filenames[i].split('_target_'))
        tag = list_of_searchids[i]
        start_time = time.time()
        with open(yaml_params, 'r') as stream:
            yaml_parameteres_for_digest = yaml.load(stream)

        # Do in silico digest....
        digest = insilicodigest.InSilicoDigest(database_path=list_of_databases[i],
                                               num_miss_cleavs=int(yaml_parameteres_for_digest['Parameters'][
                                                                       'Missed_Cleavages']),
                                               digest_type=yaml_parameteres_for_digest['Parameters'][
                                                   'Digest_Type'])
        digest.execute()

        # Initiate the reader...
        # Input for now is a target percolator output and a decoy percolator output
        pep_and_prot_data = ProteinInference.reader.PercolatorRead(target_file=target,
                                                                   decoy_file=decoy,
                                                                   yaml_param_file=yaml_params,
                                                                   digest_class=digest)

        # Execeute the reader instance, this loads the data into the reader class
        pep_and_prot_data.execute()
        # Next create a data store which is a class that stores all data for all steps of the PI process
        # Each method and each class calls from this data class to gather information and store data for analyses
        data = ProteinInference.datastore.DataStore(pep_and_prot_data)

        # Here restrict the data to having peptides with length 7 or greater
        if data.yaml_params['Parameters']['Restrict_Pep']:
            pep_restrict = float(data.yaml_params['Parameters']['Restrict_Pep'])
        else:
            pep_restrict = None
        if data.yaml_params['Parameters']['Restrict_Q']:
            q_restrict = float(data.yaml_params['Parameters']['Restrict_Q'])
        else:
            q_restrict = None

        if data.yaml_params['Parameters']['Restrict_Peptide_Length']:
            pl_restrict = int(data.yaml_params['Parameters']['Restrict_Peptide_Length'])
        else:
            pl_restrict = None

        print 'restricting data'
        restrict = ProteinInference.datastore.RestrictMainData(data, peptide_length=pl_restrict,
                                                               posterior_error_prob_threshold=pep_restrict,
                                                               q_value_threshold=q_restrict)
        restrict.execute()

        # Here generate the pre score data using 'PEP' or 'Q' values
        if data.yaml_params['Parameters']['Score_Type'] == 'pep_value':
            score_setup = ProteinInference.datastore.PreScorePepValue(data)
        if data.yaml_params['Parameters']['Score_Type'] == 'q_value':
            score_setup = ProteinInference.datastore.PreScoreQValue(data)

        # Execute score setup...
        score_setup.execute()

        score_method = data.yaml_params['Parameters']['Score_Method']

        # Here select scoring
        if score_method == 'best_peptide_per_protein':
            score = ProteinInference.scoring.BestPeptidePerProtein(data_class=data)
        if score_method == 'iterative_downweighted_log':
            score = ProteinInference.scoring.IterativeDownweightedLog(data_class=data)
        if score_method == 'multiplicative_log':
            score = ProteinInference.scoring.MultiplicativeLog(data_class=data)
        if score_method == 'downweighted_multiplicative_log':
            score = ProteinInference.scoring.DownweightedMultiplicativeLog(data_class=data)
        if score_method == 'downweighted_version2':
            score = ProteinInference.scoring.DownweightedVersion2(data_class=data)
        if score_method == 'top_two_combined':
            score = ProteinInference.scoring.TopTwoCombined(data_class=data)
        if score_method == 'geometric_mean':
            score = ProteinInference.scoring.GeometricMeanLog(data_class=data)

        # Execute scoring...
        score.execute()

        # Run protein picker on the data
        if data.yaml_params['Parameters']['Picker']:
            picker = ProteinInference.picker.StandardPicker(data_class=data)
            picker.execute()
        else:
            pass

        try:
            os.mkdir('glpkinout/')
        except OSError:
            pass

        grouping_type = data.yaml_params['Parameters']['Group']

        # Run simple group subsetting
        if grouping_type == 'simple_subsetting':
            group = ProteinInference.grouping.SimpleSubsetting(data_class=data)
            group.execute()

        server_glpk_path = '/gne/research/apps/protchem/glpk/bin/glpsol'

        # Run GLPK setup, runner, grouper...
        if grouping_type == 'glpk':
            glpksetup = ProteinInference.grouping.GlpkSetup(data_class=data,
                                                            glpkin_filename='glpkinout/glpkin_' + tag + '.mod')
            glpksetup.execute()
            glpkrun = ProteinInference.grouping.GlpkRunner(
                path_to_glpsol=data.yaml_params['Parameters']['GLPK_Path'],
                glpkin='glpkinout/glpkin_' + tag + '.mod', glpkout='glpkinout/glpkout_' + tag + '.sol',
                file_override=False)
            glpkrun.execute()
            group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest,
                                                          swissprot_override='soft',
                                                          glpksolution_filename='glpkinout/glpkout_' + tag + '.sol')
            group.execute()

        if grouping_type == 'multi_subsetting':
            group = ProteinInference.grouping.MultiSubsetting(data_class=data)
            group.execute()

        # Next run fdrcalc on the data....
        fdr = ProteinInference.fdrcalc.SetBasedFdr(data_class=data, false_discovery_rate=float(
            data.yaml_params['Parameters']['FDR']))
        fdr.execute()
        q = ProteinInference.fdrcalc.QValueCalculation(data_class=data)
        q.execute()
        # Finally we have our output restricted data...
        restricted = data.fdr_restricted_grouped_scored_proteins
        # Print the len of restricted data... which is how many protein groups pass FDR threshold
        print 'Number of Proteins passing an FDR of' + str(data.yaml_params['Parameters']['FDR']) + ' = ' + str(
            len(restricted))

        export_type = data.yaml_params['Parameters']['Export']

        # Write the output to a csv...
        if 'q_value' in export_type:
            export = ProteinInference.export.CsvOutLeadsQValues(data_class=data,
                                                         filename_out=dir_name + tag + '_' + 'leads' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export.execute()

        print 'Protein Inference Finished'



        entrap = ProteinInference.entrapment.GeneratePlot(data_class = data, entrapment_db= 'random_entrap_dbs_new/random_entrap_db_new_'+str(list_of_db_numbers[i])+'.fasta',
                                                          true_db='entrapment_data/'+list_of_true_db[i], search_id=list_of_searchids[i],pdf=pdf,
                                                          picked=data.yaml_params['Parameters']['Picker'],qvr=data.yaml_params['Parameters']['Restrict_Q'],
                                                          pvr=data.yaml_params['Parameters']['Restrict_Pep'])
        entrap.execute()


        #youden_data = data.max_youdens_data


    print 'score type = '+str(data.score_type)
    print 'score method = '+str(data.score_method)
    end_time = time.time()-start_time
    print str(end_time)+' seconds'
    print str(end_time/60)+' minutes'
    print str(end_time/(60*60))+' hours'
