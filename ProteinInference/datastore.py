#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:15:15 2017

@author: hinklet
"""

import collections 
from physical import Protein

class DataStore(object):

    #Feed data_store instance the reader class...
    def __init__(self,reader_class):
        #If the reader class is from a percolator.psms then define main_data_form as reader_class.psms
        #main_data_form is the starting point for all other analyses
        if reader_class.psms:
            self.main_data_form = reader_class.psms
            self.restricted_peptides = [x.identifier.split('.')[1] for x in self.main_data_form]
        self.protein_info_dict = None
        self.potential_proteins = None
        self.main_data_restricted = None
        self.qvalues = None
        self.pepvalues = None
        self.scored_proteins = None
        self.grouped_scored_proteins = None
        self.fdr_restricted_grouped_scored_proteins = None
        self.scoring_input = None
        self.picked_proteins_scored = None
        self.protein_peptide_dictionary = None
        self.high_low_better = None
        self.glpk_protein_number_dictionary = None
        self.glpk_number_protein_dictionary = None
        self.glpk_lead_proteins = None
        self.in_silico_digest = None
        self.max_youdens_data = None
        self.score_type = None
        self.score_method = None
        self.picked_proteins_removed = None
        self.protein_group_objects = None
        
class ProteinIdentifiers(DataStore):
    
    def __init__(self,data_class):

        #Select whether to use main_data_restricted or main_data_form... restricted is used if restrict_main_data is ever called
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
            
    def execute(self):
            
        proteins = [x.possible_proteins for x in self.data_to_use]
        self.data_class.potential_proteins = proteins
        
class QValues(DataStore):
    
    def __init__(self,data_class):
        
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
            
    def execute(self):
            
        qvalues = [x.qvalue for x in self.data_to_use]
        self.data_class.qvalues = qvalues
        
class PepValues(DataStore):
    
    def __init__(self,data_class):
        
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
            
    def execute(self):
            
        pepvalues = [x.pepvalue for x in self.data_to_use]
        self.data_class.pepvalues = pepvalues
        
class ProteinInformationDictionary(DataStore):
    
    def __init__(self,data_class):
        
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
    
    def execute(self):
            
        protein_psm_score_dictionary = collections.defaultdict(list)

        #Loop through all Psms
        for psms in self.data_to_use:
            #Loop through all proteins
            for prots in psms.possible_proteins:
                protein_psm_score_dictionary[prots].append({'peptide':psms.identifier.split('.')[1],'Qvalue':psms.qvalue,'PosteriorErrorProbability':psms.pepvalue,'Percscore':psms.percscore})
        
        self.data_class.protein_info_dict = protein_psm_score_dictionary
        
    
class RestrictMainData(DataStore):

    #This is a main and central class that allows us to restrict our PSM's on serveral criteria
    #Peptide length, posterior error probability, and q value
    #This is open and customizable... Defaults are set to peptide length of 7 and pep/qvalue are left at none
    ###If this is ran, we automatically discared everything that has a PEP of 1
    #This means it has a 100% chance of being an incorrect match...
    def __init__(self,data_class,peptide_length=7,posterior_error_prob_threshold=None,q_value_threshold=None):
        
        self.main_data = data_class.main_data_form
        self.posterior_error_prob_threshold = posterior_error_prob_threshold
        self.peptide_length = peptide_length
        self.q_value_threshold = q_value_threshold
        self.data_class = data_class
        
    def execute(self):
        #If restrict_main_data is called, we automatically discard everything that has a PEP of 1
        self.main_data = [x for x in self.main_data if x.pepvalue!=1]


        #Restrict peptide length and posterior error probability
        if self.peptide_length and self.posterior_error_prob_threshold and not self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if len(psms.identifier.split('.')[1])>=self.peptide_length and psms.pepvalue<float(self.posterior_error_prob_threshold):
                    restricted_data.append(psms)

        #Restrict peptide length only
        if self.peptide_length and not self.posterior_error_prob_threshold and not self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if len(psms.identifier.split('.')[1])>=self.peptide_length:
                    restricted_data.append(psms)

        #Restrict peptide length, posterior error probability, and qvalue
        if self.peptide_length and self.posterior_error_prob_threshold and self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if len(psms.identifier.split('.')[1])>=self.peptide_length and psms.pepvalue<float(self.posterior_error_prob_threshold) and psms.qvalue<float(self.q_value_threshold):
                    restricted_data.append(psms)

        #Restrict peptide length and qvalue
        if self.peptide_length and not self.posterior_error_prob_threshold and self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if len(psms.identifier.split('.')[1])>=self.peptide_length and psms.qvalue<float(self.q_value_threshold):
                    restricted_data.append(psms)

        #Restrict posterior error probability and q value
        if not self.peptide_length and self.posterior_error_prob_threshold and self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if psms.pepvalue<float(self.posterior_error_prob_threshold) and psms.qvalue<float(self.q_value_threshold):
                    restricted_data.append(psms)     

        #Restrict qvalue only
        if not self.peptide_length and not self.posterior_error_prob_threshold and self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if psms.qvalue<float(self.q_value_threshold):
                    restricted_data.append(psms)       

        #Restrict posterior error probability only
        if not self.peptide_length and self.posterior_error_prob_threshold and not self.q_value_threshold:
            restricted_data = []
            for psms in self.main_data:
                if psms.pepvalue<float(self.posterior_error_prob_threshold):
                    restricted_data.append(psms)

        #Restrict nothing... (only PEP gets restricted - takes everything less than 1)
        if not self.peptide_length and not self.posterior_error_prob_threshold and not self.q_value_threshold:
           restricted_data = self.main_data                   
        
        self.data_class.main_data_restricted = restricted_data

        self.data_class.restricted_peptides = [x.identifier.split('.')[1] for x in restricted_data]
        
        
class PreScoreQValue(DataStore):

    #This class generates a data structure used in protein scoring
    #If this class is called, Qvalues are used for scoring...
    def __init__(self,data_class):

        #Select which data to use... restricted or main
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
    
    def execute(self):


        protein_psm_score_dictionary = collections.defaultdict(list)
        raw_peptide_dictionary = collections.defaultdict(list)

        #Loop through all Psms
        for psms in self.data_to_use:
            #Loop through all proteins
            for prots in psms.possible_proteins:
                #Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                protein_psm_score_dictionary[prots].append({'peptide':psms.identifier.split('.')[1],'score':psms.qvalue})
                raw_peptide_dictionary[prots].append({'complete_peptide':psms.identifier})
        protein_list = []
        for pkeys in protein_psm_score_dictionary.keys():
            p = Protein(identifier = pkeys)
            p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
            p.raw_peptides = [x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]
            protein_list.append(p)

        self.data_class.score_type = 'q_value'
        self.data_class.scoring_input = protein_list
        
        
class PreScorePepValue(DataStore):

    #This class generates a data structure used in protein scoring
    #If this class is called, Pepvalues are used for scoring...
    def __init__(self,data_class):
        
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
    
    def execute(self):
            
        protein_psm_score_dictionary = collections.defaultdict(list)
        raw_peptide_dictionary = collections.defaultdict(list)


        #Loop through all Psms
        for psms in self.data_to_use:
            #Loop through all proteins
            for prots in psms.possible_proteins:
                #Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                protein_psm_score_dictionary[prots].append({'peptide':psms.identifier.split('.')[1],'score':psms.pepvalue})
                raw_peptide_dictionary[prots].append({'complete_peptide': psms.identifier})
        protein_list = []
        for pkeys in protein_psm_score_dictionary.keys():
            p = Protein(identifier = pkeys)
            p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
            p.raw_peptides = [x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]
            protein_list.append(p)

        self.data_class.score_type = 'pep_value'
        self.data_class.scoring_input = protein_list
        
class ProteinToPeptideDictionary(DataStore):
    
    def __init__(self,data_class):
        #Here we create the protein to peptide mappings with peptides
        #We create a dictionary where the key is the protein and the entries are peptides
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
            
    def execute(self):
            
        dd_prots = collections.defaultdict(set)
        for peptide_objects in self.data_to_use:
            for prots in peptide_objects.possible_proteins:
                dd_prots[prots].add(peptide_objects.identifier.split('.')[1])
        
        self.data_class.protein_peptide_dictionary = dd_prots
        
class PeptideToProteinDictionary(DataStore):
    
    def __init__(self,data_class):
        #Here we create the peptide to protein mappings with unique peptides
        #This class is useful when using the glpk based grouping
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
        
    def execute(self):
            
        dd_peps = collections.defaultdict(set)
        for peptide_objects in self.data_to_use:
            for prots in peptide_objects.possible_proteins:
                dd_peps[peptide_objects.identifier.split('.')[1]].add(prots)
                
        self.data_class.peptide_protein_dictionary = dd_peps
        
        

class HigherOrLower(DataStore):
    #Important class that is useful in determining how to sort the data based on the score...
    
    def __init__(self,data_class):

        #Get the best and worst scored proteins...
        #data_class.scored_proteins is a sorted list of all scored proteins
        #These are sorted based on the method and it is pre known what the sorting should be...
        #This class is for determining downstream sorting...
        self.worst_score = data_class.scored_proteins[-1].score
        self.best_score = data_class.scored_proteins[0].score
        self.data_class = data_class
        
    def execute(self):
        
        if float(self.best_score)>float(self.worst_score):
            higher_or_lower = 'higher'
            
        if float(self.best_score)<float(self.worst_score):
            higher_or_lower='lower'

        print 'best score = ' + str(self.best_score)
        print 'worst score = ' + str(self.worst_score)

        if self.best_score==self.worst_score:
            raise ValueError('Best and Worst scores were identical, equal to '+str(self.best_score)+'. Score type '+str(self.data_class.score_type)+' produced the error, please change score type.')
        
        self.data_class.high_low_better = higher_or_lower


        
class GetProteinIdentifiers(DataStore):

    #This class gets protein identifiers from a variety of threshold spots in the algorithm
    #One can get all identifiers from (main, restricted, picked, picked_removed, and fdr_restricted)
    #This is useful as we can see what proteins are removed at what steps of the algorithm...
    def __init__(self,data_class,data_form):
        self.data_class = data_class
        self.data_form = data_form
        self.proteins = None
               
    def execute(self):
        if self.data_form == 'main':
            #All the data (unrestricted)
            data_to_select = self.data_class.main_data_form
            prots = [[x.possible_proteins] for x in data_to_select]
            self.proteins = prots
                
        if self.data_form == 'restricted':
            #Proteins that pass certain restriction criteria (peptide length, pep, qvalue)
            data_to_select = self.data_class.main_data_restricted
            prots = [[x.possible_proteins] for x in data_to_select]
            self.proteins = prots
            
        if self.data_form == 'picked':
            #Here we look at proteins that are 'picked' (aka the proteins that beat out their matching target/decoy)
            data_to_select = self.data_class.picked_proteins_scored
            prots = [x.identifier for x in data_to_select]
            self.proteins = prots
            
        if self.data_form == 'picked_removed':
            #Here we look at the proteins that were removed due to picking (aka the proteins that have a worse score than their target/decoy counterpart)
            data_to_select = self.data_class.picked_proteins_removed
            prots = [x.identifier for x in data_to_select]
            self.proteins = prots
            
        if self.data_form == 'fdr_restricted':
            #Proteins that pass fdr restriction...
            data_to_select = self.data_class.fdr_restricted_grouped_scored_proteins
            prots = [x.identifier for x in data_to_select]
            self.proteins = prots
            

class GetProteinInformation(DataStore):

    #Once proteins are scored... Here we can look at attributes of the proteins...
    #Score, groups, reviewed/unreviewed, peptides, peptide scores, picked/removed, number of peptides...
    def __init__(self,data_class):
        self.data_class = data_class
        self.scored_data = data_class.scored_proteins

    def execute(self,protein_string):

        all_scored_protein_data = self.scored_data
        identifiers = [x.identifier for x in all_scored_protein_data]
        protein_scores = [x.score for x in all_scored_protein_data]
        groups = [x.group_identification for x in all_scored_protein_data]
        reviewed = [x.reviewed for x in all_scored_protein_data]
        peptides = [x.peptides for x in all_scored_protein_data]
        #Peptide scores currently broken...
        peptide_scores = [x.peptide_scores for x in all_scored_protein_data]
        picked = [x.picked for x in all_scored_protein_data]
        num_peptides = [x.num_peptides for x in all_scored_protein_data]

        main_index = identifiers.index(protein_string)

        list_structure = [['identifier','protein_score','groups','reviewed','peptides','peptide_scores','picked','num_peptides']]
        list_structure.append([protein_string])
        list_structure[-1].append(protein_scores[main_index])
        list_structure[-1].append(groups[main_index])
        list_structure[-1].append(reviewed[main_index])
        list_structure[-1].append(peptides[main_index])
        list_structure[-1].append(peptide_scores[main_index])
        list_structure[-1].append(picked[main_index])
        list_structure[-1].append(num_peptides[main_index])

        return list_structure

            
            
            
            
        