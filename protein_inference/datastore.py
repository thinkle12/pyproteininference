#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:15:15 2017

@author: hinklet
"""

import collections
from Bio import SeqIO
from logging import getLogger
from protein_inference.physical import Protein, Psm




class DataStore(object):
    """
    The following Class stores data at every step of the PI analysis.
    The class serves as a central point that is accessed at virtually every PI processing step

    Example: protein_inference.datastore.DataStore(reader_class = reader)

    Where reader is a Reader class object
    """

    #Feed data_store instance the reader class...
    def __init__(self,reader_class, digest_class, validate=True):
        #If the reader class is from a percolator.psms then define main_data_form as reader_class.psms
        #main_data_form is the starting point for all other analyses
        if reader_class.psms:
            self.main_data_form = reader_class.psms
            self.restricted_peptides = [Psm.split_peptide(peptide_string=x.identifier) for x in self.main_data_form]
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
        self.score = None
        self.score_method = None
        self.picked_proteins_removed = None
        self.protein_group_objects = None
        self.qvality_output = None
        self.decoy_symbol = self.parameter_file_object.decoy_symbol
        self.digest_class = digest_class
        self.SCORE_MAPPER = {"q_value":"qvalue",
                             "pep_value": "pepvalue",
                             "perc_score":"percscore",
                             "q-value": "qvalue",
                             "posterior_error_prob":"pepvalue",
                             "posterior_error_probability":"pepvalue"}
        self.CUSTOM_SCORE_KEY = "custom_score"

        self.logger = getLogger('protein_inference.datastore.DataStore')

        # Run Checks and Validations
        if validate:
            self.validate_psm_data()
            self.validate_digest(digest_class)
            self.check_data_consistency(digest_class)


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

    def restrict_psm_data(self,parameter_file_object, remove1pep=True):

        logger = getLogger('protein_inference.datastore.DataStore.restrict_psm_data')

        logger.info('Restricting PSM data')

        peptide_length = parameter_file_object.restrict_peptide_length
        posterior_error_prob_threshold = parameter_file_object.restrict_pep
        q_value_threshold = parameter_file_object.restrict_q

        main_psm_data = self.main_data_form
        logger.info('Length of main data: ' + str(len(self.main_data_form)))
        # If restrict_main_data is called, we automatically discard everything that has a PEP of 1
        if remove1pep and posterior_error_prob_threshold:
            main_psm_data = [x for x in main_psm_data if x.pepvalue != 1]

        # Restrict peptide length and posterior error probability
        if peptide_length and posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(Psm.split_peptide(peptide_string=psms.identifier)) >= peptide_length and psms.pepvalue < float(
                        posterior_error_prob_threshold):
                    restricted_data.append(psms)

        # Restrict peptide length only
        if peptide_length and not posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(Psm.split_peptide(peptide_string=psms.identifier)) >= peptide_length:
                    restricted_data.append(psms)

        # Restrict peptide length, posterior error probability, and qvalue
        if peptide_length and posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(Psm.split_peptide(peptide_string=psms.identifier)) >= peptide_length and psms.pepvalue < float(
                        posterior_error_prob_threshold) and psms.qvalue < float(q_value_threshold):
                    restricted_data.append(psms)

        # Restrict peptide length and qvalue
        if peptide_length and not posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(Psm.split_peptide(peptide_string=psms.identifier)) >= peptide_length and psms.qvalue < float(
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

        logger.info('Length of restricted data: ' + str(len(restricted_data)))

        self.restricted_peptides = [Psm.split_peptide(peptide_string=x.identifier) for x in restricted_data]
    
    def create_scoring_input(self, score_input="posterior_error_prob"):

        logger = getLogger('protein_inference.datastore.DataStore.create_scoring_input')

        logger.info("Creating Scoring Input")
        psm_data = self.get_psm_data()

        protein_psm_score_dictionary = collections.defaultdict(list)
        raw_peptide_dictionary = collections.defaultdict(list)
        psmid_peptide_dictionary = collections.defaultdict(list)

        try:
            score_key = self.SCORE_MAPPER[score_input]
        except KeyError:
            score_key = self.CUSTOM_SCORE_KEY

        if self.parameter_file_object.inference_type!="peptide_centric":
            # Loop through all Psms
            for psms in psm_data:
                # Loop through all proteins
                for prots in psms.possible_proteins:
                    # Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                    protein_psm_score_dictionary[prots].append(
                        {'peptide': Psm.split_peptide(peptide_string=psms.identifier), 'score': getattr(psms, score_key)})
                    raw_peptide_dictionary[prots].append({'complete_peptide': psms.identifier})
                    psmid_peptide_dictionary[prots].append(
                        {'peptide': Psm.split_peptide(peptide_string=psms.identifier), 'psm_id': psms.psm_id})


        else:
            self.peptide_to_protein_dictionary()

            sp_proteins = self.digest_class.swiss_prot_protein_set
            for psms in psm_data:
                protein_set = self.peptide_protein_dictionary[Psm.split_peptide(peptide_string=psms.identifier)]
                # Sort protein_set by sp-alpha, decoy-sp-alpha, tr-alpha, decoy-tr-alpha
                sorted_protein_list = self.sort_protein_strings(protein_string_list=protein_set, sp_proteins=sp_proteins)
                # Restrict the number of identifiers by the value in param file max_identifiers_peptide_centric
                sorted_protein_list = sorted_protein_list[:self.parameter_file_object.max_identifiers_peptide_centric]
                protein_name = ";".join(sorted_protein_list)

                # Generate a protein psm score dictionary for each protein... here peptides are listed as well as the selected score
                protein_psm_score_dictionary[protein_name].append(
                    {'peptide': Psm.split_peptide(peptide_string=psms.identifier), 'score': getattr(psms, score_key)})
                raw_peptide_dictionary[protein_name].append({'complete_peptide': psms.identifier})
                psmid_peptide_dictionary[protein_name].append(
                    {'peptide': Psm.split_peptide(peptide_string=psms.identifier), 'psm_id': psms.psm_id})

        protein_list = []
        # TODO Sort the protein_psm_score_dictionary by SP, Then ##SP then TR then ##TR
        for pkeys in sorted(protein_psm_score_dictionary.keys()):
            p = Protein(identifier=pkeys)
            p.psm_score_dictionary = protein_psm_score_dictionary[pkeys]
            p.psmid_peptide_dictionary = psmid_peptide_dictionary[pkeys]
            p.raw_peptides = set(sorted([x['complete_peptide'] for x in raw_peptide_dictionary[pkeys]]))
            protein_list.append(p)

        self.score = score_input
        self.scoring_input = protein_list

    def protein_to_peptide_dictionary(self):
        psm_data = self.get_psm_data()

        res_pep_set = set(self.restricted_peptides)
        dd_prots = collections.defaultdict(set)
        for peptide_objects in psm_data:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = Psm.split_peptide(peptide_string=peptide_objects.identifier)
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
                cur_peptide = Psm.split_peptide(peptide_string=peptide_objects.identifier)
                if cur_peptide in res_pep_set:
                    dd_peps[cur_peptide].add(prots)
                else:
                    pass

        self.peptide_protein_dictionary = dd_peps

        return(dd_peps)

    def unique_to_leads_peptides(self):
        if self.grouped_scored_proteins:
            lead_peptides = [list(x[0].peptides) for x in self.grouped_scored_proteins]
            flat_peptides = [item for sublist in lead_peptides for item in sublist]
            counted_peps = collections.Counter(flat_peptides)
            unique_to_leads_peptides = set([x for x in counted_peps if counted_peps[x] == 1])
        else:
            unique_to_leads_peptides=set()

        return(unique_to_leads_peptides)
        
    def higher_or_lower(self):

        logger = getLogger('protein_inference.datastore.DataStore.higher_or_lower')

        logger.info("Determining If a higher or lower score is better based on scored proteins")
        worst_score = self.scored_proteins[-1].score
        best_score = self.scored_proteins[0].score

        if float(best_score) > float(worst_score):
            higher_or_lower = 'higher'

        if float(best_score) < float(worst_score):
            higher_or_lower = 'lower'

        logger.info('best score = ' + str(best_score))
        logger.info('worst score = ' + str(worst_score))

        if best_score == worst_score:
            raise ValueError(
                'Best and Worst scores were identical, equal to ' + str(best_score) + '. Score type ' + str(
                    self.score) + ' produced the error, please change score type.')

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

        logger = getLogger('protein_inference.datastore.DataStore.exclude_non_distinguishing_peptides')

        logger.info('Applying Exclusion Model')

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
        logger.info(str(len(sets)) + ' number of peptide sets')
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
                logger.info("Parsed {} Peptide Sets".format(i))

        logger.info("Parsed {} Peptide Sets".format(i))

        # Get their index from peptides which is the initial list of sets...
        list_of_indeces = []
        for u in us:
            ind = peptides.index(u)
            list_of_indeces.append(ind)

        non_subset_proteins = set([our_proteins_sorted[x] for x in list_of_indeces])

        self.non_subset_proteins = non_subset_proteins

        logger.info("Removing direct subset Proteins from the data")
        # Remove all proteins from scoring input that are a subset of another protein...
        self.scoring_input = [x for x in self.scoring_input if
                                         x.identifier in non_subset_proteins]

        logger.info(str(len(self.scoring_input)) + ' proteins in scoring input after removing subset proteins')

        # For all the proteins that are not a complete subset of another protein...
        # Get the raw peptides...
        raw_peps = [x.raw_peptides for x in self.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        flat_peptides = [Psm.split_peptide(peptide_string=item) for sublist in raw_peps for item in sublist]

        # Count the number of peptides in this list...
        # This is the number of proteins this peptide maps to....
        counted_peptides = collections.Counter(flat_peptides)

        self.counted_peptides = counted_peptides

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
                if Psm.split_peptide(peptide_string=rp) in raw_peps_good:
                    new_raw_peptides.append(rp)

            current_score_input[j].psm_score_dictionary = new_psm_score_dictionary
            current_score_input[j].psmid_peptide_dictionary = new_psmid_peptide_dictionary
            current_score_input[j].raw_peptides = new_raw_peptides

            if k % 10000 == 0:
                logger.info("Redefined {} Peptide Sets".format(k))

        logger.info("Redefined {} Peptide Sets".format(j))

        filtered_score_input = [x for x in current_score_input if x.psm_score_dictionary]

        self.scoring_input = filtered_score_input

        # Recompute the flat peptides
        raw_peps = [x.raw_peptides for x in self.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        new_flat_peptides = set([Psm.split_peptide(peptide_string=item) for sublist in raw_peps for item in sublist])

        self.scoring_input = [x for x in self.scoring_input if x.psm_score_dictionary]

        self.restricted_peptides = [x for x in self.restricted_peptides if x in new_flat_peptides]


    def protein_picker(self):

        logger = getLogger('protein_inference.datastore.DataStore.protein_picker')

        logger.info("Running Protein Picker")

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
        logger.info(str(len(all_proteins)) + ' proteins scored')

        total_targets = []
        total_decoys = []
        decoys_removed = []
        targets_removed = []
        # Loop over all decoys identified in the search
        logger.info('Picking Proteins...')
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

        logger.info(str(len(total_decoys)) + ' total decoy proteins')
        logger.info(str(len(total_targets)) + ' matching target proteins also found in search')
        logger.info(str(len(decoys_removed)) + ' decoy proteins to be removed')
        logger.info(str(len(targets_removed)) + ' target proteins to be removed')

        logger.info('Removing Lower Scoring Proteins...')
        picked_list = []
        removed_proteins = []
        for protein_objects in self.scored_proteins:
            if protein_objects.picked == True:
                picked_list.append(protein_objects)
            else:
                removed_proteins.append(protein_objects)
        self.picked_proteins_scored = picked_list
        self.picked_proteins_removed = removed_proteins
        logger.info('Finished Removing Proteins')

    def set_based_fdr(self, false_discovery_rate = 0.01, regular=True):
        """
        function calculates set based FDR on the lead protein in the group
        Input is a DataStore object as well as an integer false discovery rate

        Example: protein_inference.fdrcalc.SetBasedFDR(data_class = data,false_discovery_rate=XX))

        FDR is calculated As (2*decoys)/total
        """

        logger = getLogger('protein_inference.datastore.DataStore.set_based_fdr')

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

        logger.info('Protein Group leads that pass with more than 1 PSM with a ' + str(
            false_discovery_rate) + ' FDR = ' + str(len(fdr_restricted_set) - len(onehitwonders)))
        logger.info('Protein Group lead One hit Wonders that pass ' + str(false_discovery_rate) + ' FDR = ' + str(
            len(onehitwonders)))

        logger.info('Number of Protein groups that pass a ' + str(false_discovery_rate * 100) + ' percent FDR: ' + str(
            len(fdr_restricted_set)))
        self.fdr_restricted_grouped_scored_proteins = fdr_restricted_set


    def calculate_q_values(self):
        """
        Class calculates Q values on the lead protein in the groups
        Input is a DataStore object

        Example: protein_inference.fdrcalc.SetBasedFDR(data_class = data)

        Q values are calculated As (2*decoys)/total
        """

        logger = getLogger('protein_inference.datastore.DataStore.calculate_q_values')

        logger.info("Calculating Q values from the protein group objects")

        #pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x.proteins[0] for x in self.protein_group_objects]
        #Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        lead_proteins.reverse()

        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_decoy_target_list = [1 if elem.startswith(self.decoy_symbol) else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            decoys = sum(binary_decoy_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            fdr = (decoys * 2) / (float(total))
            fdr_list.append(fdr)
            del lead_proteins[0]

        qvalue_list = []
        new_fdr_list = []
        for fdrs in fdr_list:
            new_fdr_list.append(fdrs)
            # qvalue = min(new_fdr_list)
            qvalue = fdrs
            qvalue_list.append(qvalue)

        qvalue_list.reverse()

        for k in range(len(self.protein_group_objects)):
            self.protein_group_objects[k].q_value = qvalue_list[k]

        logger.info("Finished Q value Calculation")


    def entrapment_fdr(self, true_database, false_discovery_rate=.05):
        """
        Class calculates Entrapment FDR on the lead protein in the groups.
        Input is a DataStore object, an entrapment database, and a false discovery rate

        Example: protein_inference.fdrcalc.EntrapFdr(data_class = data, entrapment_database = "example_entrap.fasta", other_database = None, false_discovery_rate=.05)

        FDR values are calculated As (entrapped proteins)/total

        This class is useful for calculating an entrapment FDR IF using a benchmark dataset with KNOWN protein content.
        Entrapment DB would be target proteins known to NOT be in the sample. However, these entrapped proteins should be in the main database that
        the search was searched against via comet/mascot
        """

        logger = getLogger('protein_inference.datastore.DataStore.entrapment_fdr')

        true_handle = SeqIO.parse(true_database, 'fasta')
        true_proteins = []
        for records in true_handle:
            if not records.id.startsiwith(self.decoy_symbol):  # TODO this should not be an in... should be a startswith and use params
                true_proteins.append(records.id)

        protein_data = [x[0].identifier for x in self.grouped_scored_proteins]
        false_true_positives = []
        decoys = []
        for i in range(len(protein_data)):
            # Get a list of false true positives from the protein data... basically if its not a decoy and if its not in true_proteins
            if not protein_data[i].startswith(self.decoy_symbol) and protein_data[i] not in true_proteins:
                false_true_positives.append(protein_data[i])
            if protein_data[i].startswith(self.decoy_symbol):
                decoys.append(protein_data[i])

        entrapment_proteins_set = list(false_true_positives)
        entrapment_proteins_set = set(entrapment_proteins_set)

        # pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x[0] for x in self.grouped_scored_proteins]
        # Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        # Reverse the list (best to worst) -> (worst to best)
        lead_proteins.reverse()

        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_entrap_target_list = [1 if elem in entrapment_proteins_set else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            entraped = sum(binary_entrap_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            fdr = (entraped) / (float(total))
            fdr_list.append(fdr)
            # print(fdr)
            if fdr < false_discovery_rate:
                break
            else:
                # Here we delete the worst score every time... thus making our list smaller and smaller
                del lead_proteins[0]

        lead_proteins.reverse()

        fdr_restricted_set = [self.grouped_scored_proteins[x] for x in range(len(lead_proteins))]

        self.fdr_list = fdr_list

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        logger.info('Protein Group leads that pass with more than 1 PSM with a ' + str(
            false_discovery_rate) + ' Entrapment FDR = ' + str(len(fdr_restricted_set) - len(onehitwonders)))
        logger.info('Protein Group lead One hit Wonders that pass ' + str(
            false_discovery_rate) + ' Entrapment FDR = ' + str(
            len(onehitwonders)))

        logger.info(
            'Number of Protein groups that pass a ' + str(false_discovery_rate * 100) + ' Entrapment FDR: ' + str(
                len(fdr_restricted_set)))
        self.fdr_restricted_grouped_scored_proteins = fdr_restricted_set

        self.restricted_proteins = fdr_restricted_set

        self.entrapment_proteins = false_true_positives

    def validate_psm_data(self):
        self._validate_decoys_from_data()
        self._validate_isoform_from_data()


    def validate_digest(self, digest_class):
        self._validate_reviewed_v_unreviewed(digest_class)
        self._check_target_decoy_split(digest_class)


    def check_data_consistency(self, digest_class):
        self._check_data_digest_overlap_psms(digest_class)
        self._check_data_digest_overlap_proteins(digest_class)


    def _check_data_digest_overlap_psms(self, digest_class):
        logger = getLogger('protein_inference.datastore.DataStore._check_data_digest_overlap_psms')

        ## TODO write a function here that looks at the peptides we have and checks how many of these peptides we do not find in our Digest...
        peptides = [Psm.split_peptide(peptide_string=x.identifier) for x in self.main_data_form]
        peptides_in_digest = set(digest_class.peptide_to_protein_dictionary.keys())
        peptides_from_search_in_digest = [x for x in peptides if x in peptides_in_digest]
        percentage = float(len(peptides))/float(len(peptides_from_search_in_digest))
        logger.info("{} PSMs identified from input files".format(len(peptides)))
        logger.info("{} PSMs identified from input files that are also present in database digestion".format(len(peptides_from_search_in_digest)))
        logger.info("{} percent of PSMs identified from input files that are also present in database digestion".format(percentage))

    def _check_data_digest_overlap_proteins(self, digest_class):
        logger = getLogger('protein_inference.datastore.DataStore._check_data_digest_overlap_proteins')

        proteins = [x.possible_proteins for x in self.main_data_form]
        flat_proteins = set([item for sublist in proteins for item in sublist])
        proteins_in_digest = set(digest_class.protein_to_peptide_dictionary.keys())
        proteins_from_search_in_digest = [x for x in flat_proteins if x in proteins_in_digest]
        percentage = float(len(flat_proteins))/float(len(proteins_from_search_in_digest))
        logger.info("{} proteins identified from input files".format(len(flat_proteins)))
        logger.info("{} proteins identified from input files that are also present in database digestion".format(len(proteins_from_search_in_digest)))
        logger.info("{} percent of proteins identified from input files that are also present in database digestion".format(percentage))

    def _check_target_decoy_split(self, digest_class):
        logger = getLogger('protein_inference.datastore.DataStore._check_target_decoy_split')

        # Check the number of targets vs the number of decoys from the digest
        targets = [x for x in digest_class.protein_to_peptide_dictionary.keys() if self.parameter_file_object.decoy_symbol not in x]
        decoys = [x for x in digest_class.protein_to_peptide_dictionary.keys() if self.parameter_file_object.decoy_symbol in x]
        ratio = float(len(targets)) / float(len(decoys))
        logger.info("Number of Target Proteins in Digest: {}".format(len(targets)))
        logger.info("Number of Decoy Proteins in Digest: {}".format(len(decoys)))
        logger.info("Ratio of Targets Proteins to Decoy Proteins: {}".format(ratio))

    def _validate_score(self):
        # Make sure the score is actually in our data file header... Not sure if we can do this...?
        pass

    def _validate_decoys_from_data(self):
        logger = getLogger('protein_inference.datastore.DataStore._validate_decoys_from_data')

        # Check to see if we find decoys from our input files
        proteins = [x.possible_proteins for x in self.main_data_form]
        flat_proteins = set([item for sublist in proteins for item in sublist])
        targets = [x for x in flat_proteins if self.parameter_file_object.decoy_symbol not in x]
        decoys = [x for x in flat_proteins if self.parameter_file_object.decoy_symbol in x]
        logger.info("Number of Target Proteins in Data Files: {}".format(len(targets)))
        logger.info("Number of Decoy Proteins in Data Files: {}".format(len(decoys)))


    def _validate_isoform_from_data(self):
        logger = getLogger('protein_inference.datastore.DataStore._validate_isoform_from_data')

        # Check to see if we find any proteins with isoform info in name in our input files
        proteins = [x.possible_proteins for x in self.main_data_form]
        flat_proteins = set([item for sublist in proteins for item in sublist])
        non_iso = [x for x in flat_proteins if self.parameter_file_object.isoform_symbol not in x]
        iso = [x for x in flat_proteins if self.parameter_file_object.isoform_symbol in x]
        logger.info("Number of Non Isoform Labeled Proteins in Data Files: {}".format(len(non_iso)))
        logger.info("Number of Isoform Labeled Proteins in Data Files: {}".format(len(iso)))


    def _validate_reviewed_v_unreviewed(self, digest_class):
        logger = getLogger('protein_inference.datastore.DataStore._validate_reviewed_v_unreviewed')

        # Check to see if we get reviewed prots in digest...
        reviewed_proteins = len(digest_class.swiss_prot_protein_set)
        proteins_in_digest = len(set(digest_class.protein_to_peptide_dictionary.keys()))
        unreviewed_proteins = proteins_in_digest - reviewed_proteins
        logger.info("Number of Total Proteins in from Digest: {}".format(proteins_in_digest))
        logger.info("Number of Reviewed Proteins in from Digest: {}".format(reviewed_proteins))
        logger.info("Number of Unreviewed Proteins in from Digest: {}".format(unreviewed_proteins))

    @classmethod
    def sort_protein_strings(cls,protein_string_list, sp_proteins, decoy_symbol):

        our_target_sp_proteins = sorted(
            [x for x in protein_string_list if
             x in sp_proteins and decoy_symbol not in x])
        our_decoy_sp_proteins = sorted(
            [x for x in protein_string_list if
             x in sp_proteins and decoy_symbol in x])

        our_target_tr_proteins = sorted(
            [x for x in protein_string_list if
             x not in sp_proteins and decoy_symbol not in x])
        our_decoy_tr_proteins = sorted(
            [x for x in protein_string_list if
             x not in sp_proteins and decoy_symbol in x])

        identifiers_sorted = our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins

        return(identifiers_sorted)

