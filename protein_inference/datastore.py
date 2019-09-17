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
        self.decoy_symbol = self.parameter_file_object.decoy_symbol
        self.score_mapper = {"q_value":"qvalue",
                             "pep_value": "pepvalue",
                             "perc_score":"percscore"}


    def get_sorted_identifiers(self, digest_class, scored=True):

        if scored:
            if self.picked_proteins_scored:
                proteins = set([x.identifier for x in self.picked_proteins_scored])
            else:
                proteins = set([x.identifier for x in self.scored_proteins])
        else:
            proteins = [x.identifier for x in self.scoring_input]

        all_sp_proteins = set(digest_class.swiss_prot_protein_set)

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

    def exclude_non_distinguishing_peptides(self, digest_class, protein_subset_type = "hard"):
        # TODO Move this into param file
        decoy_symbol = "##"

        print('Applying Exclusion Model')

        our_proteins_sorted = self.get_sorted_identifiers(digest_class=digest_class, scored=False)

        if protein_subset_type == "hard":
            # Hard protein subsetting defines protein subsets on the digest level (Entire protein is used)
            # This is how Percolator PI does subsetting
            peptides = [digest_class.protein_to_peptide_dictionary[x] for x in our_proteins_sorted]
        elif protein_subset_type == "soft":
            # Soft protein subsetting defines protein subsets on the Peptides identified from the search
            peptides = [set(x.raw_peptides) for x in self.scoring_input]
        else:
            # If neither is dfined we do "hard" exclusion
            peptides = [digest_class.protein_to_peptide_dictionary[x] for x in our_proteins_sorted]

        # Get frozen set of peptides....
        # We will also have a corresponding list of proteins...
        # They will have the same index...
        sets = [frozenset(e) for e in peptides]
        # Find a way to sort this list of sets...
        # We can sort the sets if we sort proteins from above...
        print(str(len(sets)) + ' number of peptide sets')
        us = set()
        i = 0
        # Get all peptide sets that are not a subset...
        while sets:
            i = i + 1
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
        self.scoring_input = [x for x in self.scoring_input if
                                         x.identifier in non_subset_proteins]

        print(str(len(self.scoring_input)) + ' proteins in scoring input after removing subset proteins')

        # For all the proteins that are not a complete subset of another protein...
        # Get the raw peptides...
        raw_peps = [x.raw_peptides for x in self.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        flat_peptides = [item.split(".")[1] for sublist in raw_peps for item in sublist]

        # Count the number of peptides in this list...
        # This is the number of proteins this peptide maps to....
        counted_peptides = collections.Counter(flat_peptides)

        # If the count is greater than 1... exclude the protein entirely from scoring input... :)
        raw_peps_good = set([x for x in counted_peptides.keys() if counted_peptides[x] <= 1])

        current_score_input = list(self.scoring_input)
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

        self.scoring_input = filtered_score_input

        # Recompute the flat peptides
        raw_peps = [x.raw_peptides for x in self.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        new_flat_peptides = set([item.split(".")[1] for sublist in raw_peps for item in sublist])

        self.scoring_input = [x for x in self.scoring_input if x.psm_score_dictionary]

        self.restricted_peptides = [x for x in self.restricted_peptides if x in new_flat_peptides]


    def protein_picker(self):
        # Use higher or lower class to determine if a higher protein score or lower protein score is better based on the scoring method used
        higher_or_lower = self.higher_or_lower()
        # Here we determine if a lower or higher score is better
        # Since all input is ordered from best to worst we can do the following

        index_to_remove = []
        # data_class.scored_proteins is simply a list of Protein objects...
        # Create list of all decoy proteins
        decoy_proteins = [x.identifier for x in self.scored_proteins if self.decoy_symbol in x.identifier]
        # Create a list of all potential matching targets (some of these may not exist in the search)
        # TODO change this to being a .replace('##','')
        # TODO we should also really put the decoy symbol in the param file...
        matching_targets = [x.replace(self.decoy_symbol,"") for x in decoy_proteins]

        # Create a list of all the proteins from the scored data
        all_proteins = [x.identifier for x in self.scored_proteins]
        print(str(len(all_proteins)) + ' proteins scored')

        total_targets = []
        total_decoys = []
        decoys_removed = []
        targets_removed = []
        # Loop over all decoys identified in the search
        print('Picking Proteins...')
        for i in range(len(decoy_proteins)):
            cur_decoy_index = all_proteins.index(decoy_proteins[i])
            cur_decoy_protein_object = self.scored_proteins[cur_decoy_index]
            total_decoys.append(cur_decoy_protein_object.identifier)

            # Try, Except here because the matching target to the decoy may not be a result from the search
            try:
                cur_target_index = all_proteins.index(matching_targets[i])
                cur_target_protein_object = self.scored_proteins[cur_target_index]
                total_targets.append(cur_target_protein_object.identifier)

                if higher_or_lower == 'higher':
                    if cur_target_protein_object.score > cur_decoy_protein_object.score:
                        index_to_remove.append(cur_decoy_index)
                        decoys_removed.append(cur_decoy_index)
                        cur_target_protein_object.picked = True
                        cur_decoy_protein_object.picked = False
                    else:
                        index_to_remove.append(cur_target_index)
                        targets_removed.append(cur_target_index)
                        cur_decoy_protein_object.picked = True
                        cur_target_protein_object.picked = False

                if higher_or_lower == 'lower':
                    if cur_target_protein_object.score < cur_decoy_protein_object.score:
                        index_to_remove.append(cur_decoy_index)
                        decoys_removed.append(cur_decoy_index)
                        cur_target_protein_object.picked = True
                        cur_decoy_protein_object.picked = False
                    else:
                        index_to_remove.append(cur_target_index)
                        targets_removed.append(cur_target_index)
                        cur_decoy_protein_object.picked = True
                        cur_target_protein_object.picked = False
            except ValueError:
                pass

        print(str(len(total_decoys)) + ' total decoy proteins')
        print(str(len(total_targets)) + ' matching target proteins also found in search')
        print(str(len(decoys_removed)) + ' decoy proteins to be removed')
        print(str(len(targets_removed)) + ' target proteins to be removed')

        print('Removing Lower Scoring Proteins...')
        picked_list = []
        removed_proteins = []
        for protein_objects in self.scored_proteins:
            if protein_objects.picked == True:
                picked_list.append(protein_objects)
            else:
                removed_proteins.append(protein_objects)
        self.picked_proteins_scored = picked_list
        self.picked_proteins_removed = removed_proteins
        print('Finished Removing Proteins')

    def set_based_fdr(self, false_discovery_rate = 0.01, regular=True):
        """
        function calculates set based FDR on the lead protein in the group
        Input is a DataStore object as well as an integer false discovery rate

        Example: protein_inference.fdrcalc.SetBasedFDR(data_class = data,false_discovery_rate=XX))

        FDR is calculated As (2*decoys)/total
        """
        # pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x[0] for x in self.grouped_scored_proteins]
        # Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        # Reverse the list (best to worst) -> (worst to best)
        lead_proteins.reverse()

        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_decoy_target_list = [1 if elem.startswith(self.decoy_symbol) else 0 for elem in
                                        lead_proteins]
            total = len(lead_proteins)
            decoys = sum(binary_decoy_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            if regular:
                fdr = (2 * decoys) / (float(total))
            else:
                fdr = (decoys) / (float(total))
            fdr_list.append(fdr)
            # print(fdr)
            if fdr < false_discovery_rate:
                break
            else:
                # Here we delete the worst score every time... thus making our list smaller and smaller
                del lead_proteins[0]

        lead_proteins.reverse()

        self.fdr_list = fdr_list

        fdr_restricted_set = [self.grouped_scored_proteins[x] for x in range(len(lead_proteins))]

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        print('Protein Group leads that pass with more than 1 PSM with a ' + str(
            false_discovery_rate) + ' FDR = ' + str(len(fdr_restricted_set) - len(onehitwonders)))
        print('Protein Group lead One hit Wonders that pass ' + str(false_discovery_rate) + ' FDR = ' + str(
            len(onehitwonders)))

        print('Number of Protein groups that pass a ' + str(false_discovery_rate * 100) + ' percent FDR: ' + str(
            len(fdr_restricted_set)))
        self.fdr_restricted_grouped_scored_proteins = fdr_restricted_set


