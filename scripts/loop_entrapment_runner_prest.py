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

param_tags = "KEEP_NoPicker_standard_new"
RemoveNonUniquePeptides = False


db_path = "/Users/hinklet/PythonPackages/protein_inference/old/all_benchmark_entrapment_data/entrapment_data/"


combined_targets = [['147858_percolator_target_psm_all_psms-ae0ab4ec-78ab-4092-a505-1f5434b55a70.txt',
 '147859_percolator_target_psm_all_psms-1869d69c-b11a-410e-9f8e-43c133752c7e.txt',
 '147860_percolator_target_psm_all_psms-c7b08d19-689b-423e-96ff-1815d6257e91.txt'],
 ['147861_percolator_target_psm_all_psms-2b382860-425c-4672-b315-2bff0a04021e.txt',
 '147862_percolator_target_psm_all_psms-c1198f31-71aa-4857-b19c-26c82206711d.txt',
 '147863_percolator_target_psm_all_psms-8b2d061b-53b2-4475-b827-4f8b13410aa6.txt'],
 ['147864_percolator_target_psm_all_psms-37f39262-194b-4cbe-af62-c790d121c3c5.txt',
 '147865_percolator_target_psm_all_psms-44a039c5-c3e7-4b44-8db7-914da9904df6.txt',
 '147866_percolator_target_psm_all_psms-fc39906c-da4d-4bee-8ad7-09d5f2a268d4.txt']]

combined_decoys = [['147858_percolator_decoy_psm_all_psms-ae0ab4ec-78ab-4092-a505-1f5434b55a70.txt',
 '147859_percolator_decoy_psm_all_psms-1869d69c-b11a-410e-9f8e-43c133752c7e.txt',
 '147860_percolator_decoy_psm_all_psms-c7b08d19-689b-423e-96ff-1815d6257e91.txt'],
 ['147861_percolator_decoy_psm_all_psms-2b382860-425c-4672-b315-2bff0a04021e.txt',
 '147862_percolator_decoy_psm_all_psms-c1198f31-71aa-4857-b19c-26c82206711d.txt',
 '147863_percolator_decoy_psm_all_psms-8b2d061b-53b2-4475-b827-4f8b13410aa6.txt'],
 ['147864_percolator_decoy_psm_all_psms-37f39262-194b-4cbe-af62-c790d121c3c5.txt',
 '147865_percolator_decoy_psm_all_psms-44a039c5-c3e7-4b44-8db7-914da9904df6.txt',
 '147866_percolator_decoy_psm_all_psms-fc39906c-da4d-4bee-8ad7-09d5f2a268d4.txt']]


mapper3 = {'Rep123_A':"1000_random_plus_b.fasta",
           'Rep123_B':"1000_random_plus_a.fasta",
           'Rep123_AB':"prest_1000_random.fasta"
}

list_of_tags = ['Rep123_B','Rep123_A','Rep123_AB']

list_of_search_filenames_all = os.listdir('/Users/hinklet/new_benchmark_percolator/')

list_of_search_filenames = [x for x in list_of_search_filenames_all if 'combined' in x ]
list_of_target_filenames = [x for x in list_of_search_filenames_all if 'target' in x ]
list_of_decoy_filenames = [x for x in list_of_search_filenames_all if 'decoy' in x ]
list_of_true_db = ['complete_and_reversed_AB.fasta' for x in range(len(list_of_search_filenames))]

list_of_searchids = [x.split('/')[-1].split('_percolator_')[0] for x in list_of_target_filenames]

#list_of_db_numbers = [mapper[x] for x in list_of_searchids]

list_of_db_numbers = ["B", "A", "AB"]

list_of_databases = [db_path+"prest_pool_"+x.lower()+".fasta" for x in list_of_db_numbers]

entrap_list = [None]
entrap_list2 = [None]

dir_name = "/Users/hinklet/PI_output_benchmark/"
yaml_params = "/Users/hinklet/PythonPackages/protein_inference/parameters/Protein_Inference_Params_prest.yaml"
import protein_inference
import time

with open(yaml_params, 'r') as stream:
    yaml_parameteres_for_digest = yaml.load(stream)


# Make this only a 3x for loop not a 9x for loop... waste of time...
with PdfPages('plots/'+yaml_parameteres_for_digest['Parameters']['Score_Method'] + '_' + yaml_parameteres_for_digest['Parameters']['Score_Type'] +'_' +param_tags + '.pdf') as pdf:
    for i in range(len(combined_targets)):
        target = ['/Users/hinklet/new_benchmark_percolator/' + x for x in combined_targets[i]]
        decoy =['/Users/hinklet/new_benchmark_percolator/' + x for x in combined_decoys[i]]
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
        pep_and_prot_data = protein_inference.reader.PercolatorRead(target_file=target,
                                                                   decoy_file=decoy,
                                                                   yaml_param_file=yaml_params,
                                                                   digest_class=digest)

        # Execeute the reader instance, this loads the data into the reader class
        pep_and_prot_data.execute()
        # Next create a data store which is a class that stores all data for all steps of the PI process
        # Each method and each class calls from this data class to gather information and store data for analyses
        data = protein_inference.datastore.DataStore(pep_and_prot_data)

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
        restrict = protein_inference.datastore.RestrictMainData(data, peptide_length=pl_restrict,
                                                               posterior_error_prob_threshold=pep_restrict,
                                                               q_value_threshold=q_restrict)
        restrict.execute()

        # Here generate the pre score data using 'PEP' or 'Q' values
        if data.yaml_params['Parameters']['Score_Type'] == 'pep_value':
            score_setup = protein_inference.datastore.PreScorePepValue(data)
        if data.yaml_params['Parameters']['Score_Type'] == 'q_value':
            score_setup = protein_inference.datastore.PreScoreQValue(data)

        # Execute score setup...
        score_setup.execute()

        if RemoveNonUniquePeptides:
            rnup = protein_inference.datastore.RemoveNonUniquePeptides(data)
            rnup.execute()

        score_method = data.yaml_params['Parameters']['Score_Method']

        # Here select scoring
        if score_method == 'best_peptide_per_protein':
            score = protein_inference.scoring.BestPeptidePerProtein(data_class=data)
        if score_method == 'iterative_downweighted_log':
            score = protein_inference.scoring.IterativeDownweightedLog(data_class=data)
        if score_method == 'multiplicative_log':
            score = protein_inference.scoring.MultiplicativeLog(data_class=data)
        if score_method == 'downweighted_multiplicative_log':
            score = protein_inference.scoring.DownweightedMultiplicativeLog(data_class=data)
        if score_method == 'downweighted_version2':
            score = protein_inference.scoring.DownweightedVersion2(data_class=data)
        if score_method == 'top_two_combined':
            score = protein_inference.scoring.TopTwoCombined(data_class=data)
        if score_method == 'geometric_mean':
            score = protein_inference.scoring.GeometricMeanLog(data_class=data)

        # Execute scoring...
        score.execute()

        # Run protein picker on the data
        if data.yaml_params['Parameters']['Picker']:
            picker = protein_inference.picker.StandardPicker(data_class=data)
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
            group = protein_inference.grouping.SimpleSubsetting(data_class=data)
            group.execute()

        server_glpk_path = '/gne/research/apps/protchem/glpk/bin/glpsol'

        # Run GLPK setup, runner, grouper...
        if grouping_type == 'glpk':
            glpksetup = protein_inference.grouping.GlpkSetup(data_class=data,
                                                            glpkin_filename='glpkinout/glpkin_' + tag +'_'+param_tags+ '.mod')
            glpksetup.execute()
            glpkrun = protein_inference.grouping.GlpkRunner(
                path_to_glpsol=data.yaml_params['Parameters']['GLPK_Path'],
                glpkin='glpkinout/glpkin_' + tag +'_'+param_tags+ '.mod', glpkout='glpkinout/glpkout_' + tag +'_'+param_tags+ '.sol',
                file_override=False)
            glpkrun.execute()
            group = protein_inference.grouping.GlpkGrouper(data_class=data, digest_class=digest,
                                                          swissprot_override='soft',
                                                          glpksolution_filename='glpkinout/glpkout_' + tag +'_'+param_tags+ '.sol')
            group.execute()

        if grouping_type == 'multi_subsetting':
            group = protein_inference.grouping.MultiSubsetting(data_class=data)
            group.execute()

        # Next run fdrcalc on the data....
        fdr = protein_inference.fdrcalc.SetBasedFdr(data_class=data, false_discovery_rate=float(
            data.yaml_params['Parameters']['FDR']))
        fdr.execute()
        q = protein_inference.fdrcalc.QValueCalculation(data_class=data)
        q.execute()
        # Finally we have our output restricted data...
        restricted = data.fdr_restricted_grouped_scored_proteins
        # Print the len of restricted data... which is how many protein groups pass FDR threshold
        print 'Number of Proteins passing an FDR of' + str(data.yaml_params['Parameters']['FDR']) + ' = ' + str(
            len(restricted))

        export_type = data.yaml_params['Parameters']['Export']

        # Write the output to a csv...
        if 'q_value' in export_type:
            export = protein_inference.export.CsvOutCommaSepQValues(data_class=data,
                                                         filename_out=dir_name + param_tags + "_" + tag + '_' + 'CommaSep' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export.execute()

            export2 = protein_inference.export.CsvOutLeadsQValues(data_class=data,
                                                                   filename_out=dir_name + param_tags + "_" + tag + '_' + 'Leads' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export2.execute()

            export3 = protein_inference.export.CsvOutAllQValues(data_class=data,
                                                                 filename_out=dir_name + param_tags + "_" + tag + '_' + 'All' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            export3.execute()

            # export4 = protein_inference.export.CsvOutAll(data_class=data,
            #                                                    filename_out=dir_name + param_tags + "_" + tag + '_' + 'All' + '_' + data.short_score_method + '_' + data.score_type + '.csv')
            # export4.execute()


        print 'Protein Inference Finished'



        entrap = protein_inference.entrapment.GeneratePlot(data_class = data, true_db=list_of_databases[i], search_id=tag,pdf=pdf,
                                                          picked=data.yaml_params['Parameters']['Picker'],qvr=data.yaml_params['Parameters']['Restrict_Q'],pvr=data.yaml_params['Parameters']['Restrict_Pep'])
        entrap.execute()


        #youden_data = data.max_youdens_data


print 'score type = '+str(data.score_type)
print 'score method = '+str(data.score_method)
end_time = time.time()-start_time
print str(end_time)+' seconds'
print str(end_time/60)+' minutes'
print str(end_time/(60*60))+' hours'
