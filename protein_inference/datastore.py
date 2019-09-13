#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:15:15 2017

@author: hinklet
"""

import collections 
from protein_inference.physical import Protein
import yaml
import re
from collections import OrderedDict


class DataStore(object):
    """
    The following Class stores data at every step of the PI analysis.
    The class serves as a central point that is accessed at virtually every PI processing step

    Example: protein_inference.datastore.DataStore(reader_class = reader)

    Where reader is a Reader class object
    """

    #Feed data_store instance the reader class...
    def __init__(self,reader_class):
        #If the reader class is from a percolator.psms then define main_data_form as reader_class.psms
        #main_data_form is the starting point for all other analyses
        if reader_class.psms:
            self.main_data_form = reader_class.psms
            self.restricted_peptides = [x.identifier.split('.')[1] for x in self.main_data_form]
            if not reader_class.search_id:
                if "Bioplex" in reader_class.target_file:
                    self.search_id = reader_class.target_file.split('_')[0]
                if "Bioplex" not in reader_class.target_file:
                    try:
                        self.search_id = reader_class.target_file.split('_')[0]
                    except AttributeError:
                        self.search_id = 'Custom'
            else:
                self.search_id = reader_class.search_id



            #This is bad because default GLPK_Path is glpsol... on rescomp this will not work...

        self.yaml_params = reader_class.yaml_params
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
        self.qvality_output = None
        
        self.decoy_symbol = "##"


    def get_sorted_identifiers(self, digest_class, scored=True):

        if scored:
            if self.picked_proteins_scored:
                proteins = set([x.identifier for x in self.picked_proteins_scored])
            else:
                proteins = set([x.identifier for x in self.scored_proteins])
        else:
            proteins = [x.identifier for x in self.scoring_input]

        all_sp_proteins = set(digest_class.swiss_prot_protein_dictionary['swiss-prot'])

        our_target_sp_proteins = sorted([x for x in proteins if x in all_sp_proteins and self.decoy_symbol not in x])
        our_decoy_sp_proteins = sorted([x for x in proteins if x in all_sp_proteins and self.decoy_symbol in x])

        our_target_tr_proteins = sorted([x for x in proteins if x not in all_sp_proteins and self.decoy_symbol not in x])
        our_decoy_tr_proteins = sorted([x for x in proteins if x not in all_sp_proteins and self.decoy_symbol in x])

        our_proteins_sorted = our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins

        return(our_proteins_sorted)
class ProteinIdentifiers(DataStore):
    """
    Class that outputs all Protein Identifiers and stores it as "potential_proteins" in the DataStore object

    Example: protein_inference.datastore.ProteinIdenfitiers(data_class = data)

    Where data is a DataStore Object.

    Running execute will store the protein list and it can be accessed as:

    data.potential_proteins - where data again is a DataStore Object
    """
    
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
    """
    Class that outputs all Qvalues and stores it as "qvalues" in the DataStore object

    Example: protein_inference.datastore.QValues(data_class = data)

    Where data is a DataStore Object.

    Running execute will store the qvalue list and it can be accessed as:

    data.qvalues - where data again is a DataStore Object
    """
    
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
    """
    Class that outputs all PepValues and stores it as "pepvalues" in the DataStore object

    Example: protein_inference.datastore.PepValues(data_class = data)

    Where data is a DataStore Object.

    Running execute will store the pepvalue list and it can be accessed as:

    data.pepvalues - where data again is a DataStore Object
    """
    
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
    """
    Class that creates a protein information dictionary and stores it as "protein_info_dict" in the DataStore object

    Example: protein_inference.datastore.ProteinInformationDictionary(data_class = data)

    Where data is a DataStore Object.

    Running execute will store the dictionary and it can be accessed as:

    data.protein_info_dict - where data again is a DataStore Object

    Stored in the dictionary is the peptide, Qvlaue, PepValue, and PercScore
    """

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
                protein_psm_score_dictionary[prots].append({'peptide':psms.identifier,'Qvalue':psms.qvalue,'PosteriorErrorProbability':psms.pepvalue,'Percscore':psms.percscore})
        
        self.data_class.protein_info_dict = protein_psm_score_dictionary
        
    
class RestrictMainData(DataStore):

    """
    This is a main and central class that allows us to restrict our PSM's on serveral criteria

    Peptide length, posterior error probability, and q value

    This is open and customizable... Defaults are set to peptide length of 7 and pep/qvalue are left at none

    If this is ran, we automatically discared everything that has a PEP of 1
    This means it has a 100% chance of being an incorrect match...

    Example: protein_inference.datastore.RestrictMainData(data_class = data,peptide_length=7,posterior_error_prob_threshold=None,q_value_threshold=None)

    Again, data is a DataStore Object

    It is highly recommended to run this class with default parameters before pre score is ran.
    """

    def __init__(self,data_class,peptide_length=7,posterior_error_prob_threshold=None,q_value_threshold=None):
        
        self.main_data = data_class.main_data_form
        self.posterior_error_prob_threshold = posterior_error_prob_threshold
        self.peptide_length = peptide_length
        self.q_value_threshold = q_value_threshold
        self.data_class = data_class
        
    def execute(self):
        print('Length of main data: '+str(len(self.main_data)))
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

        print('Length of restricted data: '+str(len(restricted_data)))


        self.data_class.restricted_peptides = [x.identifier.split('.')[1] for x in restricted_data]
        
        
class PreScoreQValue(DataStore):
    """
    This class generates the generic data structure used in protein scoring
    If this class is called, Qvalues are used for scoring...

    Example: protein_inference.datastore.PreScoreQValue(data_class = data)

    Where data is a DataStore Object
    """
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
        psmid_peptide_dictionary = collections.defaultdict(list)


        #Loop through all Psms
        for psms in self.data_to_use:
            #Loop through all proteins
            for prots in psms.possible_proteins:
                #Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                protein_psm_score_dictionary[prots].append({'peptide':psms.identifier.split('.')[1],'score':psms.qvalue})
                raw_peptide_dictionary[prots].append({'complete_peptide':psms.identifier})
                psmid_peptide_dictionary[prots].append({'peptide':psms.identifier.split('.')[1],'psm_id':psms.psm_id})

        # print 'Using Generators for Pre Scoring...'
        # def gen():
        #     for pkeys in protein_psm_score_dictionary.keys():
        #         p = Protein(identifier=pkeys)
        #         p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
        #         p.psmid_peptide_dictionary = psmid_peptide_dictionary[pkeys]
        #         p.raw_peptides = [x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]
        #         yield p

        protein_list = []
        for pkeys in sorted(protein_psm_score_dictionary.keys()):
            p = Protein(identifier = pkeys)
            p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
            p.psmid_peptide_dictionary = psmid_peptide_dictionary[pkeys]
            p.raw_peptides = set(sorted([x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]))
            protein_list.append(p)

        self.data_class.score_type = 'q_value'
        self.data_class.scoring_input = protein_list
        
        
class PreScorePepValue(DataStore):
    """
    This class generates the generic data structure used in protein scoring
    If this class is called, Pepvalues are used for scoring...

    Example: protein_inference.datastore.PreScorePepValue(data_class = data)

    Where data is a DataStore Object
    """
    def __init__(self,data_class):
        
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
    
    def execute(self):
            
        protein_psm_score_dictionary = collections.defaultdict(list)
        raw_peptide_dictionary = collections.defaultdict(list)
        psmid_peptide_dictionary = collections.defaultdict(list)


        #Loop through all Psms
        for psms in self.data_to_use:
            #Loop through all proteins
            for prots in psms.possible_proteins:
                #Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                protein_psm_score_dictionary[prots].append({'peptide':psms.identifier.split('.')[1],'score':psms.pepvalue})
                raw_peptide_dictionary[prots].append({'complete_peptide': psms.identifier})
                psmid_peptide_dictionary[prots].append({'peptide': psms.identifier.split('.')[1], 'psm_id': psms.psm_id})
        protein_list = []
        # TODO Sort the protein_psm_score_dictionary by SP, Then ##SP then TR then ##TR
        for pkeys in sorted(protein_psm_score_dictionary.keys()):
            p = Protein(identifier = pkeys)
            p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
            p.psmid_peptide_dictionary = psmid_peptide_dictionary[pkeys]
            p.raw_peptides = set(sorted([x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]))
            protein_list.append(p)

        self.data_class.score_type = 'pep_value'
        self.data_class.scoring_input = protein_list
        
class ProteinToPeptideDictionary(DataStore):
    """
    This class creates the Protein to Peptide dictionary for peptides from the search results

    Example: protein_inference.datastore.ProteinToPeptideDictionary(data_class = data)

    Where data is a DataStore Object

    This class provides a quick and easy map of Proteins to peptides from the search and is used in protein_inference.grouping
    """
    
    def __init__(self,data_class):
        #Here we create the protein to peptide mappings with peptides
        #We create a dictionary where the key is the protein and the entries are peptides
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
            
    def execute(self):

        res_pep_set = set(self.data_class.restricted_peptides)
        dd_prots = collections.defaultdict(set)
        for peptide_objects in self.data_to_use:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = peptide_objects.identifier.split('.')[1]
                if cur_peptide in res_pep_set:
                    dd_prots[prots].add(cur_peptide)
        
        self.data_class.protein_peptide_dictionary = dd_prots
        
class PeptideToProteinDictionary(DataStore):
    """
    This class creates the Peptide to Protein dictionary for peptides from the search results

    Example: protein_inference.datastore.PeptideToProteinDictionary(data_class = data)

    Where data is a DataStore Object

    This class provides a quick and easy map of Peptides to Proteins from the search and is used in protein_inference.grouping.GlpkSetup
    """
    
    def __init__(self,data_class):
        #Here we create the peptide to protein mappings with unique peptides
        #This class is useful when using the glpk based grouping
        if data_class.main_data_restricted:
            self.data_to_use = data_class.main_data_restricted
        else:
            self.data_to_use = data_class.main_data_form
        self.data_class = data_class
        
    def execute(self):


        res_pep_set = set(self.data_class.restricted_peptides)
        dd_peps = collections.defaultdict(set)
        for peptide_objects in self.data_to_use:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = peptide_objects.identifier.split('.')[1]
                if cur_peptide in res_pep_set:
                    dd_peps[cur_peptide].add(prots)
                else:
                    pass
                
        self.data_class.peptide_protein_dictionary = dd_peps
        
        

class HigherOrLower(DataStore):
    """
    Important class that is useful in determining how to sort the data based on the score...

    This class is simple in that it takes the best score from a scoring output and the worst score for a scoring output and compares them.

    It then assigns an internal variable known as high_low_better to showcase whether a higher or lower score is better based on the scoring scheme chosen.

    This is important as whenever results are sorted by score, we need to know whether a higher or a lower score is better.
    Given that a higher score or a lower score can be seen as "better" for different scoring types, we implement this class to determine if higher or lower is better

    Example: protein_inference.datastore.HigherOrLower(data_class = data)

    where data is a DataStore Object
    """

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

        print('best score = ' + str(self.best_score))
        print('worst score = ' + str(self.worst_score))

        if self.best_score==self.worst_score:
            raise ValueError('Best and Worst scores were identical, equal to '+str(self.best_score)+'. Score type '+str(self.data_class.score_type)+' produced the error, please change score type.')
        
        self.data_class.high_low_better = higher_or_lower


        
class GetProteinIdentifiers(DataStore):

    """
    This class gets protein identifiers from a variety of threshold spots in the algorithm
    One can get all identifiers from ['main', 'restricted', 'picked', 'picked_removed', and 'fdr_restricted']
    This is useful as we can see what proteins are removed from the data at what steps of the algorithm...

    Example: protein_inference.datastore.GetProteinIdentifiers(data_class = data, data_form = 'restricted')

    data_class is a DataStore object
    data_form can be one of the following = ['main', 'restricted', 'picked', 'picked_removed', and 'fdr_restricted']
    """
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
    """
    This class displays protein information once proteins are scored
    we can look at the following attributes of the proteins:
    Score, groups, reviewed/unreviewed, peptides, peptide scores, picked/removed, number of peptides...

    Example: gpi = protein_inference.datastore.GetProteinInformation(data_class = data)

             gpi.execute(protein_string = PRKDC_HUMAN|P78527)

    This code currently does not support outputting the peptide scores (unsure why)

    where data is a DataStore Object


    """
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

class RemoveMods(DataStore):
    """
    This class simply removes modifications from peptide strings
            
    """
    def __init__(self,peptide_string):
        self.peptide_string = peptide_string.upper()
        self.regex = re.compile('[^a-zA-Z]')

    def execute(self):

        # First parameter is the replacement, second parameter is your input string
        stripped_peptide = self.regex.sub('', self.peptide_string)
        self.stripped_peptide = stripped_peptide


class Exclusion(DataStore):
    """
    This class removes non unique peptides from adding to a proteins score...
    However, if the peptides from two proteins are identical, we keep all peptides
    """

    def __init__(self,data_class, digest_class, protein_subset_type = "hard"):
        self.data_class = data_class
        self.digest_class = digest_class
        self.protein_subset_type = protein_subset_type
        self.decoy_symbol = "##"

    def execute(self):
        print('Applying Exclusion Model')

        # Get all peptides and protein identifiers from scoring input
        # Get the proteins and sort them if we can...
        # Sort for SP first... then sort either number of peptides or alphabetically...
        proteins = [x.identifier for x in self.data_class.scoring_input]

        all_sp_proteins = set(self.digest_class.swiss_prot_protein_dictionary['swiss-prot'])

        our_target_sp_proteins = sorted([x for x in proteins if x in all_sp_proteins and self.decoy_symbol  not in x])
        our_decoy_sp_proteins = sorted([x for x in proteins if x in all_sp_proteins and self.decoy_symbol  in x])

        our_target_tr_proteins = sorted([x for x in proteins if x not in all_sp_proteins and self.decoy_symbol  not in x])
        our_decoy_tr_proteins = sorted([x for x in proteins if x not in all_sp_proteins and self.decoy_symbol   in x])

        our_proteins_sorted = our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins


        if self.protein_subset_type=="hard":
            # Hard protein subsetting defines protein subsets on the digest level (Entire protein is used)
            # This is how Percolator PI does subsetting
            peptides = [self.digest_class.protein_to_peptide_dictionary[x] for x in our_proteins_sorted]
        elif self.protein_subset_type=="soft":
            # Soft protein subsetting defines protein subsets on the Peptides identified from the search
            peptides = [set(x.raw_peptides) for x in self.data_class.scoring_input]
        else:
            # If neither is dfined we do "hard" exclusion
            peptides = [self.digest_class.protein_to_peptide_dictionary[x] for x in our_proteins_sorted]

        # Get frozen set of peptides....
        # We will also have a corresponding list of proteins...
        # They will have the same index...
        sets = [frozenset(e) for e in peptides]
        # Find a way to sort this list of sets...
        # We can sort the sets if we sort proteins from above...
        print(str(len(sets))+' number of peptide sets')
        us = set()
        i = 0
        # Get all peptide sets that are not a subset...
        while sets:
            i = i+1
            e = sets.pop()
            if any(e.issubset(s) for s in sets) or any(e.issubset(s) for s in us):
                continue
            else:
                us.add(e)
            if i % 10000 == 0:
                print("Parsed {} Peptide Sets".format(i))

        print("Parsed {} Peptide Sets".format(i))

        # Get their index from peptides which is the initial list of sets...
        list_of_indeces = []
        for u in us:
            ind = peptides.index(u)
            list_of_indeces.append(ind)

        non_subset_proteins = set([our_proteins_sorted[x] for x in list_of_indeces])

        print("Removing direct subset Proteins from the data")
        # Remove all proteins from scoring input that are a subset of another protein...
        self.data_class.scoring_input = [x for x in self.data_class.scoring_input if x.identifier in non_subset_proteins]

        print(str(len(self.data_class.scoring_input))+' proteins in scoring input after removing subset proteins')

        # For all the proteins that are not a complete subset of another protein...
        # Get the raw peptides...
        raw_peps = [x.raw_peptides for x in self.data_class.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        flat_peptides = [item.split(".")[1] for sublist in raw_peps for item in sublist]

        # Count the number of peptides in this list...
        # This is the number of proteins this peptide maps to....
        counted_peptides = collections.Counter(flat_peptides)

        # If the count is greater than 1... exclude the protein entirely from scoring input... :)
        raw_peps_good = set([x for x in counted_peptides.keys() if counted_peptides[x]<=1])


        current_score_input = list(self.data_class.scoring_input)
        for j in range(len(current_score_input)):
            k = j + 1
            new_psm_score_dictionary = []
            new_psmid_peptide_dictionary = []
            new_raw_peptides = []
            current_psm_score_dictionary = current_score_input[j].psm_score_dictionary
            current_psmid_peptide_dictionary = current_score_input[j].psmid_peptide_dictionary
            current_raw_peptides = current_score_input[j].raw_peptides

            for psm_scores in current_psm_score_dictionary:
                if psm_scores['peptide'] in raw_peps_good:
                    new_psm_score_dictionary.append(psm_scores)

            for psm_id in current_psmid_peptide_dictionary:
                if psm_id['peptide'] in raw_peps_good:
                    new_psmid_peptide_dictionary.append(psm_id)

            for rp in current_raw_peptides:
                if rp.split(".")[1] in raw_peps_good:
                    new_raw_peptides.append(rp)

            current_score_input[j].psm_score_dictionary = new_psm_score_dictionary
            current_score_input[j].psmid_peptide_dictionary = new_psmid_peptide_dictionary
            current_score_input[j].raw_peptides = new_raw_peptides

            if k % 10000 == 0:
                print("Redefined {} Peptide Sets".format(k))

        print("Redefined {} Peptide Sets".format(j))

        filtered_score_input = [x for x in current_score_input if x.psm_score_dictionary]

        self.data_class.scoring_input = filtered_score_input

        # Recompute the flat peptides
        raw_peps = [x.raw_peptides for x in self.data_class.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        new_flat_peptides = set([item.split(".")[1] for sublist in raw_peps for item in sublist])

        self.data_class.scoring_input = [x for x in self.data_class.scoring_input if x.psm_score_dictionary]



        self.data_class.restricted_peptides = [x for x in self.data_class.restricted_peptides if x in new_flat_peptides]


