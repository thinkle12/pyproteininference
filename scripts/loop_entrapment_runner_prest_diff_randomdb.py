#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:04:08 2017

@author: hinklet
"""
import matplotlib
from digest import insilicodigest
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

param_tags = "RNUP_NoPicker_78rdb_standard_new"
RemoveNonUniquePeptides = True


db_path = "/Users/hinklet/PythonPackages/py_protein_inference/old/all_benchmark_entrapment_data/entrapment_data/"



combined_targets = [['181049_percolator_target_psm_all_psms-3403b50e-f906-4806-a6b7-4eff7daac98f.txt',
                     '181050_percolator_target_psm_all_psms-dbdc690a-88c2-4662-a200-e6403efdd932.txt',
                     '181051_percolator_target_psm_all_psms-83d5a11e-5b09-46f9-be04-cc3293ab43dc.txt'],
                    ['181052_percolator_target_psm_all_psms-b3dacac3-c786-429a-965e-9e978b6e295f.txt',
                     '181053_percolator_target_psm_all_psms-545ad8f2-f1bd-43f3-b082-eca80c57650f.txt',
                     '181054_percolator_target_psm_all_psms-6614140f-a997-44d5-b39a-c84d10d0d33b.txt'],
                    ['181055_percolator_target_psm_all_psms-a1c1a34d-fc80-40c9-ae16-35e9e202576a.txt',
                     '181056_percolator_target_psm_all_psms-b4b471a2-0fb3-4b52-b687-a559e118a0ec.txt',
                     '181057_percolator_target_psm_all_psms-9269c0e2-26b8-41d8-b5f5-7005f4378309.txt']
                    ]




combined_decoys = [['181049_percolator_decoy_psm_all_psms-3403b50e-f906-4806-a6b7-4eff7daac98f.txt',
                     '181050_percolator_decoy_psm_all_psms-dbdc690a-88c2-4662-a200-e6403efdd932.txt',
                     '181051_percolator_decoy_psm_all_psms-83d5a11e-5b09-46f9-be04-cc3293ab43dc.txt'],
                   ['181052_percolator_decoy_psm_all_psms-b3dacac3-c786-429a-965e-9e978b6e295f.txt',
                     '181053_percolator_decoy_psm_all_psms-545ad8f2-f1bd-43f3-b082-eca80c57650f.txt',
                     '181054_percolator_decoy_psm_all_psms-6614140f-a997-44d5-b39a-c84d10d0d33b.txt'],
                   ['181055_percolator_decoy_psm_all_psms-a1c1a34d-fc80-40c9-ae16-35e9e202576a.txt',
                     '181056_percolator_decoy_psm_all_psms-b4b471a2-0fb3-4b52-b687-a559e118a0ec.txt',
                     '181057_percolator_decoy_psm_all_psms-9269c0e2-26b8-41d8-b5f5-7005f4378309.txt']
                   ]

# mapper2 = {'181055':"prest_1000_random.fasta",
#           '181056':"prest_1000_random.fasta",
#           '181057':"prest_1000_random.fasta",
#           '147863':"random_entrap_B_test.fasta",
#           '147862':"random_entrap_B_test.fasta",
#           '147861':"random_entrap_B_test.fasta",
#           '147860':"random_entrap_A_test.fasta",
#           '147859':"random_entrap_A_test.fasta",
#           '147858':"random_entrap_A_test.fasta"}

mapper3 = {'Rep123_A':'random_entrap_db_new_78_withB.fasta',
           'Rep123_B':'random_entrap_db_new_78_withA.fasta',
           'Rep123_AB':'random_entrap_db_new_78.fasta'
}




list_of_tags = ['Rep123_B','Rep123_A','Rep123_AB']

list_of_true_db = ['complete_random_entrap_db_AB_78_new.fasta' for x in range(3)]

list_of_db_numbers = ["B", "A", "AB"]



list_of_databases = [db_path+'prest_pool_'+x.lower()+".fasta" for x in list_of_db_numbers]

entrap_list = [None]
entrap_list2 = [None]

dir_name = '/Users/hinklet/prest_random_db_benchmark/protein_inference_files/'
yaml_params = "/Users/hinklet/PythonPackages/py_protein_inference/parameters/Protein_Inference_Params_prest.yaml"
import py_protein_inference
import time

with open(yaml_params, 'r') as stream:
    yaml_parameteres_for_digest = yaml.load(stream)


# Make this only a 3x for loop not a 9x for loop... waste of time...
with PdfPages('/Users/hinklet/prest_random_db_benchmark/plots/'+yaml_parameteres_for_digest['Parameters']['Score_Method'] + '_' + yaml_parameteres_for_digest['Parameters']['Score_Type'] +'_' +param_tags + '.pdf') as pdf:
    for i in range(len(combined_targets)):
        target = ['/Users/hinklet/prest_random_db_benchmark/percolator_out/' + x for x in combined_targets[i]]
        decoy =['/Users/hinklet/prest_random_db_benchmark/percolator_out/' + x for x in combined_decoys[i]]
        tag = list_of_tags[i]
        entrapdb = mapper3[tag]
        start_time = time.time()
        with open(yaml_params, 'r') as stream:
            yaml_parameteres_for_digest = yaml.load(stream)

        # Do in silico digest....
        digest = insilicodigest.InSilicoDigest(database_path=db_path+list_of_true_db[i],
                                               num_miss_cleavs=int(yaml_parameteres_for_digest['Parameters'][
                                                                       'Missed_Cleavages']),
                                               digest_type=yaml_parameteres_for_digest['Parameters'][
                                                   'Digest_Type'],
                                               id_splitting=False)
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

        if RemoveNonUniquePeptides:
            rnup = ProteinInference.datastore.RemoveNonUniquePeptides(data)
            rnup.execute()

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
                                                            glpkin_filename='glpkinout/glpkin_' + tag +'_'+param_tags+ '.mod')
            glpksetup.execute()
            glpkrun = ProteinInference.grouping.GlpkRunner(
                path_to_glpsol=data.yaml_params['Parameters']['GLPK_Path'],
                glpkin='glpkinout/glpkin_' + tag +'_'+param_tags+ '.mod', glpkout='glpkinout/glpkout_' + tag +'_'+param_tags+ '.sol',
                file_override=False)
            glpkrun.execute()
            group = ProteinInference.grouping.GlpkGrouper(data_class=data, digest_class=digest,
                                                          swissprot_override='soft',
                                                          glpksolution_filename='glpkinout/glpkout_' + tag +'_'+param_tags+ '.sol')
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
            export = ProteinInference.export.CsvOutCommaSepQValues(data_class=data,
                                                         filename_out=dir_name + param_tags + "_" + tag + '_' + 'CommaSep' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export.execute()

            export2 = ProteinInference.export.CsvOutLeadsQValues(data_class=data,
                                                                   filename_out=dir_name + param_tags + "_" + tag + '_' + 'Leads' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export2.execute()

            export3 = ProteinInference.export.CsvOutAllQValues(data_class=data,
                                                                 filename_out=dir_name + param_tags + "_" + tag + '_' + 'All' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export3.execute()

            # export4 = ProteinInference.export.CsvOutAll(data_class=data,
            #                                                    filename_out=dir_name + param_tags + "_" + tag + '_' + 'All' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            # export4.execute()


        print 'Protein Inference Finished'



        entrap = ProteinInference.entrapment.GeneratePlot(data_class = data, true_db=list_of_databases[i], search_id=tag,pdf=pdf,
                                                          picked=data.yaml_params['Parameters']['Picker'],qvr=data.yaml_params['Parameters']['Restrict_Q'],pvr=data.yaml_params['Parameters']['Restrict_Pep'])
        entrap.execute()


        #youden_data = data.max_youdens_data


print 'score type = '+str(data.score_type)
print 'score method = '+str(data.score_method)
end_time = time.time()-start_time
print str(end_time)+' seconds'
print str(end_time/60)+' minutes'
print str(end_time/(60*60))+' hours'
