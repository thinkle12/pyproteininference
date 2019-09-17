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
from Bio import SeqIO



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

        self.parameter_file_object = reader_class.parameter_file_object
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

    @classmethod
    def sort_proteins(cls):
        # TODO make a method that sorts entire lists of proteins or leads...
        pass

    @classmethod
    def sort_protein_sub_groups(cls):
        # TODO make a method that sorts one protein group
        pass

    @classmethod
    def sort_protein_groups(cls):
        # TODO make a method that sorts physical protein groups by their leads...
        pass

    def higher_or_lower(self):
        # TODO make this an instance method of the datastore class...
        pass

    def get_psm_data(self):
        if self.main_data_restricted:
            psm_data = self.main_data_restricted
        else:
            psm_data = self.main_data_form

        return(psm_data)

    def get_protein_data(self):
        if self.picked_proteins_scored:
            scored_proteins = self.picked_proteins_scored
        else:
            scored_proteins = self.scored_proteins

        return(scored_proteins)

    def get_protein_identifiers_from_psm_data(self):
        psm_data = self.get_psm_data()

        proteins = [x.possible_proteins for x in psm_data]

        return(proteins)

    def get_q_values(self):
        psm_data = self.get_psm_data()

        q_values = [x.qvalue for x in psm_data]

        return(q_values)

    def get_pep_values(self):
        psm_data = self.get_psm_data()

        pep_values = [x.qvalue for x in psm_data]

        return(pep_values)

    def get_protein_information_dictionary(self):
        psm_data = self.get_psm_data()

        protein_psm_score_dictionary = collections.defaultdict(list)

        # Loop through all Psms
        for psms in psm_data:
            # Loop through all proteins
            for prots in psms.possible_proteins:
                protein_psm_score_dictionary[prots].append(
                    {'peptide': psms.identifier, 'Qvalue': psms.qvalue, 'PosteriorErrorProbability': psms.pepvalue,
                     'Percscore': psms.percscore})

        self.protein_info_dict = protein_psm_score_dictionary

        return(protein_psm_score_dictionary)

    def restrict_psm_data(self,parameter_file_object):
        print('restricting data')

        peptide_length = parameter_file_object.restrict_peptide_length
        posterior_error_prob_threshold = parameter_file_object.restrict_pep
        q_value_threshold = parameter_file_object.restrict_q

        main_psm_data = self.main_data_form
        print('Length of main data: ' + str(len(self.main_data_form)))
        # If restrict_main_data is called, we automatically discard everything that has a PEP of 1
        main_psm_data = [x for x in main_psm_data if x.pepvalue != 1]

        # Restrict peptide length and posterior error probability
        if peptide_length and posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.identifier.split('.')[1]) >= peptide_length and psms.pepvalue < float(
                        posterior_error_prob_threshold):
                    restricted_data.append(psms)

        # Restrict peptide length only
        if peptide_length and not posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.identifier.split('.')[1]) >= peptide_length:
                    restricted_data.append(psms)

        # Restrict peptide length, posterior error probability, and qvalue
        if peptide_length and posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.identifier.split('.')[1]) >= peptide_length and psms.pepvalue < float(
                        posterior_error_prob_threshold) and psms.qvalue < float(q_value_threshold):
                    restricted_data.append(psms)

        # Restrict peptide length and qvalue
        if peptide_length and not posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.identifier.split('.')[1]) >= peptide_length and psms.qvalue < float(
                        q_value_threshold):
                    restricted_data.append(psms)

        # Restrict posterior error probability and q value
        if not peptide_length and posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if psms.pepvalue < float(posterior_error_prob_threshold) and psms.qvalue < float(
                        q_value_threshold):
                    restricted_data.append(psms)

                    # Restrict qvalue only
        if not peptide_length and not posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if psms.qvalue < float(q_value_threshold):
                    restricted_data.append(psms)

                    # Restrict posterior error probability only
        if not peptide_length and posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if psms.pepvalue < float(posterior_error_prob_threshold):
                    restricted_data.append(psms)

        # Restrict nothing... (only PEP gets restricted - takes everything less than 1)
        if not peptide_length and not posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = main_psm_data

        self.main_data_restricted = restricted_data

        print('Length of restricted data: ' + str(len(restricted_data)))

        self.restricted_peptides = [x.identifier.split('.')[1] for x in restricted_data]
    
    def create_scoring_input(self, score_input="pep_value"):
        psm_data = self.get_psm_data()

        protein_psm_score_dictionary = collections.defaultdict(list)
        raw_peptide_dictionary = collections.defaultdict(list)
        psmid_peptide_dictionary = collections.defaultdict(list)

        # Loop through all Psms
        for psms in psm_data:
            # Loop through all proteins
            for prots in psms.possible_proteins:
                # Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                protein_psm_score_dictionary[prots].append(
                    {'peptide': psms.identifier.split('.')[1], 'score': getattr(psms, self.score_mapper[score_input])})
                raw_peptide_dictionary[prots].append({'complete_peptide': psms.identifier})
                psmid_peptide_dictionary[prots].append(
                    {'peptide': psms.identifier.split('.')[1], 'psm_id': psms.psm_id})
        protein_list = []
        # TODO Sort the protein_psm_score_dictionary by SP, Then ##SP then TR then ##TR
        for pkeys in sorted(protein_psm_score_dictionary.keys()):
            p = Protein(identifier=pkeys)
            p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
            p.psmid_peptide_dictionary = psmid_peptide_dictionary[pkeys]
            p.raw_peptides = set(sorted([x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]))
            protein_list.append(p)

        self.score_type = score_input
        self.scoring_input = protein_list

    def protein_to_peptide_dictionary(self):
        psm_data = self.get_psm_data()

        res_pep_set = set(self.restricted_peptides)
        dd_prots = collections.defaultdict(set)
        for peptide_objects in psm_data:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = peptide_objects.identifier.split('.')[1]
                if cur_peptide in res_pep_set:
                    dd_prots[prots].add(cur_peptide)

        self.protein_peptide_dictionary = dd_prots

        return(dd_prots)

    def peptide_to_protein_dictionary(self):
        psm_data = self.get_psm_data()

        res_pep_set = set(self.restricted_peptides)
        dd_peps = collections.defaultdict(set)
        for peptide_objects in psm_data:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = peptide_objects.identifier.split('.')[1]
                if cur_peptide in res_pep_set:
                    dd_peps[cur_peptide].add(prots)
                else:
                    pass

        self.peptide_protein_dictionary = dd_peps

        return(dd_peps)
        
    def higher_or_lower(self):
        worst_score = self.scored_proteins[-1].score
        best_score = self.scored_proteins[0].score

        if float(best_score) > float(worst_score):
            higher_or_lower = 'higher'

        if float(best_score) < float(worst_score):
            higher_or_lower = 'lower'

        print('best score = ' + str(best_score))
        print('worst score = ' + str(worst_score))

        if best_score == worst_score:
            raise ValueError(
                'Best and Worst scores were identical, equal to ' + str(best_score) + '. Score type ' + str(
                    self.score_type) + ' produced the error, please change score type.')

        self.high_low_better = higher_or_lower

        return(higher_or_lower)
        

    def get_protein_identifiers(self, data_form):
        if data_form == 'main':
            # All the data (unrestricted)
            data_to_select = self.main_data_form
            prots = [[x.possible_proteins] for x in data_to_select]
            proteins = prots

        if data_form == 'restricted':
            # Proteins that pass certain restriction criteria (peptide length, pep, qvalue)
            data_to_select = self.main_data_restricted
            prots = [[x.possible_proteins] for x in data_to_select]
            proteins = prots

        if data_form == 'picked':
            # Here we look at proteins that are 'picked' (aka the proteins that beat out their matching target/decoy)
            data_to_select = self.picked_proteins_scored
            prots = [x.identifier for x in data_to_select]
            proteins = prots

        if data_form == 'picked_removed':
            # Here we look at the proteins that were removed due to picking (aka the proteins that have a worse score than their target/decoy counterpart)
            data_to_select = self.picked_proteins_removed
            prots = [x.identifier for x in data_to_select]
            proteins = prots

        if data_form == 'fdr_restricted':
            # Proteins that pass fdr restriction...
            data_to_select = self.fdr_restricted_grouped_scored_proteins
            prots = [x.identifier for x in data_to_select]
            proteins = prots

        return(proteins)
        
    def get_protein_information(self, protein_string):
        all_scored_protein_data = self.scored_proteins
        identifiers = [x.identifier for x in all_scored_protein_data]
        protein_scores = [x.score for x in all_scored_protein_data]
        groups = [x.group_identification for x in all_scored_protein_data]
        reviewed = [x.reviewed for x in all_scored_protein_data]
        peptides = [x.peptides for x in all_scored_protein_data]
        # Peptide scores currently broken...
        peptide_scores = [x.peptide_scores for x in all_scored_protein_data]
        picked = [x.picked for x in all_scored_protein_data]
        num_peptides = [x.num_peptides for x in all_scored_protein_data]

        main_index = identifiers.index(protein_string)

        list_structure = [['identifier', 'protein_score', 'groups', 'reviewed', 'peptides', 'peptide_scores', 'picked',
                           'num_peptides']]
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


