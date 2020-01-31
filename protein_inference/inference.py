#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:51:50 2017

@author: hinklet
"""

import os
import subprocess
import collections
###Should I just import ProteinInference instead?
from protein_inference import datastore
from protein_inference.physical import ProteinGroup
from protein_inference.physical import Psm
import re
from collections import OrderedDict
from logging import getLogger
import pulp


class Inference(object):
    """
    Parent Inference class for all grouper subset classes
    """

    INFERENCE_TYPES = ["parsimony", "inclusion", "exclusion", "none", "peptide_centric", None]
    GROUPING_TYPES = ["subset_peptides", "shared_peptides", "none", None]
    LP_SOLVERS = ["pulp", "glpk", "none", "None"]
    
    def __init__(self):
        None

    def run_inference(self, data_class, digest_class):

        logger = getLogger('protein_inference.inference.Inference.run_inference')

        inference_type = data_class.parameter_file_object.inference_type

        logger.info("Running Inference with Inference Type: {}".format(inference_type))

        # For parsimony... Run GLPK setup, runner, grouper...
        if inference_type == 'parsimony':
            group = Parsimony(data_class=data_class, digest_class=digest_class)
            group.infer_proteins()

        if inference_type == "inclusion":
            group = Inclusion(data_class=data_class, digest_class=digest_class)
            group.infer_proteins()

        if inference_type == "exclusion":
            group = Exclusion(data_class=data_class, digest_class=digest_class)
            group.infer_proteins()

        if inference_type == "none":
            group = FirstProtein(data_class=data_class, digest_class=digest_class)
            group.infer_proteins()

    def _group_by_peptides(self, scored_data, data_class, digest_class, inference_type="parsimony", lead_protein_objects=None, grouping_type="shared_peptides"):

        logger = getLogger('protein_inference.inference.Inference._group_by_peptides')

        logger.info("Grouping Peptides with Grouping Type: {}".format(grouping_type))
        logger.info("Grouping Peptides with Inference Type: {}".format(inference_type))

        scored_data = sorted(scored_data, key=lambda k: (len(k.raw_peptides),k.identifier), reverse=True)

        if inference_type== "parsimony":
            scored_proteins = lead_protein_objects
            scored_proteins = sorted(scored_proteins, key=lambda k: (len(k.raw_peptides),k.identifier), reverse=True)
        else:
            scored_proteins = list(scored_data)
            scored_proteins = sorted(scored_proteins, key=lambda k: (len(k.raw_peptides),k.identifier), reverse=True)
            # Add some code here to remove complete subset proteins...
        protein_finder = [x.identifier for x in scored_data]

        prot_pep_dict = data_class.protein_to_peptide_dictionary()

        protein_tracker = set()
        restricted_peptides_set = set(data_class.restricted_peptides)
        try:
            picked_removed = set([x.identifier for x in data_class.picked_proteins_removed])
        except TypeError:
            picked_removed = set([])
        list_of_prots_not_in_db = []
        list_of_peps_not_in_db = []
        missing_proteins = set()
        in_silico_peptides_to_proteins = digest_class.peptide_to_protein_dictionary
        grouped_proteins = []
        for protein_objects in scored_proteins:
            if protein_objects not in protein_tracker:
                protein_tracker.add(protein_objects)
                cur_protein_identifier = protein_objects.identifier

                # Set peptide variable if the peptide is in the restricted peptide set
                # Sort the peptides alphabetically
                protein_objects.peptides = set(sorted([x for x in prot_pep_dict[cur_protein_identifier] if
                                                       x in restricted_peptides_set]))
                protein_list_group = [protein_objects]
                current_peptides = prot_pep_dict[cur_protein_identifier]

                current_grouped_proteins = set()
                for peptide in current_peptides:  # Probably put an if here... if peptide is in the list of peptide after being restricted by datastore.RestrictMainData
                    if peptide in restricted_peptides_set:
                        # Get the proteins that map to the current peptide using in_silico_peptides_to_proteins
                        # First make sure our peptide is formatted properly...
                        if not peptide.isupper() or not peptide.isalpha():
                            # If the peptide is not all upper case or if its not all alphabetical...
                            peptide = Psm.remove_peptide_mods(peptide)
                        potential_protein_list = in_silico_peptides_to_proteins[peptide]
                        if not potential_protein_list:
                            list_of_prots_not_in_db.append(peptide)
                            list_of_peps_not_in_db.append(protein_objects.identifier)
                            logger.warning('Protein ' + str(protein_objects.identifier) + ' and Peptide ' + str(
                                peptide) + ' is not in database...')
                        # Assign proteins to groups based on shared peptide... unless the protein is equivalent to the current identifier
                        for protein in potential_protein_list:
                            # If statement below to avoid grouping the same protein twice and to not group the lead
                            if protein not in current_grouped_proteins and protein != cur_protein_identifier and protein not in picked_removed and protein not in missing_proteins and inference_type!="none":
                                try:
                                    # Try to find its object using protein_finder (list of identifiers) and scored_proteins (list of Protein Objects)
                                    cur_index = protein_finder.index(protein)
                                    current_protein_object = scored_data[cur_index]
                                    if not current_protein_object.peptides:
                                        current_protein_object.peptides = set(
                                            sorted([x for x in prot_pep_dict[current_protein_object.identifier] if
                                                    x in restricted_peptides_set]))
                                    if grouping_type=="shared_peptides":
                                        current_grouped_proteins.add(current_protein_object)
                                    elif grouping_type=="subset_peptides":
                                        if current_protein_object.peptides.issubset(protein_objects.peptides):
                                            current_grouped_proteins.add(current_protein_object)
                                            # if inference_type!="inclusion":
                                            protein_tracker.add(current_protein_object)
                                        else:
                                            pass
                                    else:
                                        pass
                                except ValueError:
                                    logger.warning("Protein from DB {} not found with protein finder for peptide {}".format(protein,
                                                                                                                   peptide))
                                    missing_proteins.add(protein)

                            else:
                                pass
                # Add the proteins to the lead if they share peptides...
                protein_list_group = protein_list_group + list(current_grouped_proteins)
                # protein_list_group at first is just the lead protein object...
                # We then try apply grouping by looking at all peptide from the lead...
                # For all of these peptide look at all other non lead proteins and try to assign them to the group...
                # We assign the entire protein object as well... in the above try/except
                # Then append this sub group to the main list
                # The variable grouped_proteins is now a list of lists which each element being a Protein object and each list of protein objects corresponding to a group
                grouped_proteins.append(protein_list_group)

        return_dict = {"grouped_proteins": grouped_proteins, "missing_proteins": list_of_prots_not_in_db,
                       "missing_peptides": list_of_peps_not_in_db}

        return (return_dict)

    def _apply_group_ids_no_groups(self, grouped_protein_objects, digest_class, data_class):

        logger = getLogger('protein_inference.inference.Inference._apply_group_ids_no_groups')

        sp_protein_set = set(digest_class.swiss_prot_protein_set)

        prot_pep_dict = data_class.protein_to_peptide_dictionary()

        score_dd = collections.defaultdict(list)

        # Here we create group ID's
        group_id = 0
        list_of_group_objects = []
        for protein_group in grouped_protein_objects:
            protein_list = []
            group_id = group_id + 1
            pg = ProteinGroup(group_id)
            logger.info("Created Protein Group with ID: {}".format(str(group_id)))
            for prot in protein_group:
                cur_protein = prot
                # The following loop assigns group_id's, reviewed/unreviewed status, and number of unique peptides...
                if group_id not in cur_protein.group_identification:
                    cur_protein.group_identification.add(group_id)
                if cur_protein.identifier in sp_protein_set:
                    cur_protein.reviewed = True
                else:
                    cur_protein.unreviewed = True
                cur_identifier = cur_protein.identifier
                cur_protein.num_peptides = len(prot_pep_dict[cur_identifier])
                # Here append the number of unique peptides... so we can use this as secondary sorting...
                protein_list.append(cur_protein)
                # Sorted protein_groups then becomes a list of lists... of protein objects

            pg.proteins = protein_list
            list_of_group_objects.append(pg)


        return_dict = {"scores_grouped": grouped_protein_objects, "group_objects":list_of_group_objects}

        return(return_dict)


    def _apply_group_ids_peptide_centric(self, data_class):

        logger = getLogger('protein_inference.inference.Inference._apply_group_ids_peptide_centric')

        grouped_protein_objects = data_class.get_protein_data()

        # Here we create group ID's
        group_id = 0
        list_of_proteins_grouped = []
        list_of_group_objects = []
        for protein_group in grouped_protein_objects:
            # TODO, This split needs to get moved to an external function or method and should be backed by param file
            protein_group.peptides = set([x.split(".")[1] for x in list(protein_group.raw_peptides)])
            protein_list = []
            group_id = group_id + 1
            pg = ProteinGroup(group_id)
            logger.info("Created Protein Group with ID: {}".format(str(group_id)))
            # The following loop assigns group_id's, reviewed/unreviewed status, and number of unique peptides...
            if group_id not in protein_group.group_identification:
                protein_group.group_identification.add(group_id)
            protein_group.num_peptides = len(protein_group.peptides)
            # Here append the number of unique peptides... so we can use this as secondary sorting...
            protein_list.append(protein_group)
            # Sorted protein_groups then becomes a list of lists... of protein objects

            pg.proteins = protein_list
            list_of_group_objects.append(pg)
            list_of_proteins_grouped.append([protein_group])


        return_dict = {"scores_grouped": list_of_proteins_grouped, "group_objects":list_of_group_objects}

        return(return_dict)


    def _apply_protein_group_ids(self, grouped_protein_objects, data_class, digest_class):
        logger = getLogger('protein_inference.inference.Inference._apply_protein_group_ids')

        sp_protein_set = set(digest_class.swiss_prot_protein_set)

        prot_pep_dict = data_class.protein_to_peptide_dictionary()

        # Here we create group ID's
        group_id = 0
        list_of_group_objects = []
        for protein_group in grouped_protein_objects:
            protein_list = []
            group_id = group_id + 1
            pg = ProteinGroup(group_id)
            logger.info("Created Protein Group with ID: {}".format(str(group_id)))
            for protein in protein_group:
                cur_protein = protein
                # The following loop assigns group_id's, reviewed/unreviewed status, and number of unique peptides...
                if group_id not in cur_protein.group_identification:
                    cur_protein.group_identification.add(group_id)
                if protein.identifier in sp_protein_set:
                    cur_protein.reviewed = True
                else:
                    cur_protein.unreviewed = True
                cur_identifier = protein.identifier
                cur_protein.num_peptides = len(prot_pep_dict[cur_identifier])
                # Here append the number of unique peptides... so we can use this as secondary sorting...
                protein_list.append(cur_protein)
                # Sorted protein_groups then becomes a list of lists... of protein objects

            pg.proteins = protein_list
            list_of_group_objects.append(pg)

        return_dict = {"scores_grouped": grouped_protein_objects, "group_objects":list_of_group_objects}


        return(return_dict)

    def _swissprot_and_isoform_override(self, scored_data, grouped_proteins, data_class, digest_class, override_type="soft", isoform_override=True):
        logger = getLogger('protein_inference.inference.Inference._swissprot_and_isoform_override')

        sp_protein_set = set(digest_class.swiss_prot_protein_set)
        scored_proteins = list(scored_data)
        protein_finder = [x.identifier for x in scored_proteins]

        prot_pep_dict = data_class.protein_to_peptide_dictionary()


        # Get the higher or lower variable
        if not data_class.high_low_better:
            higher_or_lower = data_class.higher_or_lower()
        else:
            higher_or_lower = data_class.high_low_better


        logger.info('Applying Group IDs... and Executing {} Swissprot Override...'.format(override_type))
        # Here we create group ID's for all groups and do some sorting
        scores_grouped = []
        group_id = 0
        leads = set()
        lead_replaced_prot_pairs = []
        list_of_group_objects = []
        for protein_group in grouped_proteins:
            protein_list = []
            group_id = group_id + 1
            # Make a protein group
            pg = ProteinGroup(group_id)
            logger.info("Created Protein Group with ID: {}".format(str(group_id)))
            for prots in protein_group:
                # Loop over all proteins in the original group
                try:
                    # The following loop assigns group_id's, reviewed/unreviewed status, and number of unique peptides...
                    pindex = protein_finder.index(prots.identifier)
                    # Attempt to find the protein object by identifier
                    cur_protein = scored_proteins[pindex]
                    if group_id not in cur_protein.group_identification:
                        cur_protein.group_identification.add(group_id)
                    if prots.identifier in sp_protein_set:
                        cur_protein.reviewed = True
                    else:
                        cur_protein.unreviewed = True
                    cur_identifier = prots.identifier
                    cur_protein.num_peptides = len(prot_pep_dict[cur_identifier])
                    # Here append the number of unique peptides... so we can use this as secondary sorting...
                    protein_list.append(cur_protein)
                    # Sorted groups then becomes a list of lists... of protein objects


                except ValueError:
                    # Here we pass if the protein does not have a score...
                    # Potentially it got 'picked' (removed) by protein picker...
                    pass

            # Sort the groups based on higher or lower indication, secondarily sort the groups based on number of unique peptides
            # We use the index [1:] as we do not wish to sort the lead protein... from GLPK
            if higher_or_lower == 'lower':
                protein_list[1:] = sorted(protein_list[1:],
                                        key=lambda k: (float(k.score), -float(k.num_peptides)), reverse=False)
            if higher_or_lower == 'higher':
                protein_list[1:] = sorted(protein_list[1:],
                                        key=lambda k: (float(k.score), float(k.num_peptides)), reverse=True)

            # scores_grouped is the MAIN list of lists with grouped protein objects
            scores_grouped.append(protein_list)
            # If the lead is reviewed append it to leads and do nothing else...
            # If the lead is unreviewed then try to replace it with the best reviewed hit
            if not protein_list[0].reviewed:
                # If the lead is unreviewed attempt to replace it...
                # Start to loop through protein_list which is the current group...
                for hits in protein_list[1:]:
                    # Find the first reviewed it... if its not a lead protein already then do score swap and break...
                    if hits.reviewed:
                        best_swiss_prot_prot = hits

                        if override_type=="soft":
                            # If the lead proteins peptides are a subset of the best swissprot.... then swap the proteins... (meaning equal peptides or the swissprot completely covers the tremble reference)
                            if best_swiss_prot_prot.identifier not in leads and set(protein_list[0].peptides).issubset(
                                    set(best_swiss_prot_prot.peptides)):
                                # We use -1 as the idex of scores_grouped because the current 'protein_list' is the last entry appended to scores grouped
                                # Essentially scores_grouped[-1]==protein_list
                                # We need this syntax so we can switch the location of the unreviewed lead identifier with the best reviewed identifier in scores_grouped
                                swiss_prot_override_index = scores_grouped[-1].index(best_swiss_prot_prot)
                                cur_tr_lead = scores_grouped[-1][0]
                                logger.info(cur_tr_lead.identifier)
                                scores_grouped[-1][0], scores_grouped[-1][swiss_prot_override_index] = scores_grouped[-1][
                                                                                                           swiss_prot_override_index], \
                                                                                                       scores_grouped[-1][0]
                                scores_grouped[-1][swiss_prot_override_index], scores_grouped[-1][0]
                                new_sp_lead = scores_grouped[-1][0]
                                logger.info(new_sp_lead.identifier)
                                lead_replaced_prot_pairs.append([cur_tr_lead, new_sp_lead])
                                # Append new_sp_lead protein to leads, to make sure we dont repeat leads
                                leads.add(new_sp_lead.identifier)
                                break
                            else:
                                # If no reviewed and none not in leads then pass...
                                pass

                        if override_type=="hard":
                            if best_swiss_prot_prot.identifier not in leads:
                                # We use -1 as the index of scores_grouped because the current 'protein_list' is the last entry appended to scores_grouped
                                # Essentially scores_grouped[-1]==protein_list
                                # We need this syntax so we can switch the location of the unreviewed lead identifier with the best reviewed identifier in scores_grouped
                                swiss_prot_override_index = scores_grouped[-1].index(best_swiss_prot_prot)
                                cur_tr_lead = scores_grouped[-1][0]
                                logger.info(cur_tr_lead.identifier)
                                scores_grouped[-1][0], scores_grouped[-1][swiss_prot_override_index] = \
                                scores_grouped[-1][swiss_prot_override_index], scores_grouped[-1][0]
                                new_sp_lead = scores_grouped[-1][0]
                                logger.info(new_sp_lead.identifier)
                                lead_replaced_prot_pairs.append([cur_tr_lead, new_sp_lead])
                                # Append new_sp_lead protein to leads, to make sure we dont repeat leads
                                leads.add(new_sp_lead.identifier)
                                break
                            else:
                                # If no reviewed and none not in leads then pass...
                                pass

                    else:
                        pass
            ###NEW ISOFORM OVERRIDE
            if isoform_override:
                if protein_list[0].reviewed:
                    if data_class.parameter_file_object.isoform_symbol in protein_list[0].identifier:
                        pure_id = protein_list[0].identifier.split(data_class.parameter_file_object.isoform_symbol)[0]
                        # Start to loop through protein_list which is the current group...
                        for potential_replacement in protein_list[1:]:
                            isoform_override = potential_replacement
                            if isoform_override.identifier == pure_id and isoform_override.identifier not in leads and set(
                                    protein_list[0].peptides).issubset(set(isoform_override.peptides)):
                                isoform_override_index = scores_grouped[-1].index(isoform_override)
                                cur_iso_lead = scores_grouped[-1][0]
                                logger.info(cur_iso_lead.identifier)
                                scores_grouped[-1][0], scores_grouped[-1][isoform_override_index] = scores_grouped[-1][
                                                                                                        isoform_override_index], \
                                                                                                    scores_grouped[-1][0]
                                scores_grouped[-1][isoform_override_index], scores_grouped[-1][0]
                                new_iso_lead = scores_grouped[-1][0]
                                logger.info(new_iso_lead.identifier)
                                lead_replaced_prot_pairs.append([cur_iso_lead, new_iso_lead])
                                leads.add(protein_list[0].identifier)

            pg.proteins = protein_list
            list_of_group_objects.append(pg)

        return_dict = {"lead_replaced": lead_replaced_prot_pairs, "scores_grouped": scores_grouped, "group_objects":list_of_group_objects}

        return(return_dict)


    def _reassign_leads(self, list_of_group_objects, data_class):

        logger = getLogger('protein_inference.inference.Inference._reassign_leads')

        # Get the higher or lower variable
        if not data_class.high_low_better:
            higher_or_lower = data_class.higher_or_lower()
        else:
            higher_or_lower = data_class.high_low_better

        # TODO why is this only occuring on the list_of group_objects...
        # TODO This should be happening to the list of protein objects tooo.... aka scores_grouped
        # Sometimes we have cases where:
        # protein a maps to peptides 1,2,3
        # protein b maps to peptides 1,2
        # protein c maps to a bunch of peptides and peptide 3
        # Therefore, in the model proteins a and b are equivalent in that they map to 2 peptides together - 1 and 2. peptide 3 maps to a but also to c...
        # Sometimes the model (glpk) will spit out protein b as the lead... we wish to swap protein b as the lead with protein a because it will likely have a better score...
        logger.info('Potentially Reassigning leads...')
        lead_protein_set = set([x.proteins[0].identifier for x in list_of_group_objects])
        for i in range(len(list_of_group_objects)):
            for j in range(1, len(list_of_group_objects[i].proteins)):  # Loop over all sub proteins in the group...
                # if the lead proteins peptides are a subset of one of its proteins in the group, and the secondary protein is not a lead protein and its score is better than the leads... and it has more peptides...
                new_lead = list_of_group_objects[i].proteins[j]
                old_lead = list_of_group_objects[i].proteins[0]
                if higher_or_lower == 'higher':
                    if set(old_lead.peptides).issubset(set(new_lead.peptides)) and new_lead.identifier not in lead_protein_set and old_lead.score <= new_lead.score and len(old_lead.peptides) < len(new_lead.peptides):
                        logger.info('protein {} will replace protein {} as lead, with index {}, New Num Peptides: {}, Old Num Peptides: {}'.format(str(new_lead.identifier),str(old_lead.identifier),str(j),str(len(new_lead.peptides)),str(len(old_lead.peptides))))
                        lead_protein_set.add(new_lead.identifier)
                        lead_protein_set.remove(old_lead.identifier)
                        # Swap their positions in the list
                        list_of_group_objects[i].proteins[0], list_of_group_objects[i].proteins[j] = new_lead, old_lead
                        logger.info(j)
                        break


                if higher_or_lower == 'lower':
                    if set(old_lead.peptides).issubset(set(new_lead.peptides)) and new_lead.identifier not in lead_protein_set and old_lead.score >= new_lead.score and len(old_lead.peptides) < len(new_lead.peptides):
                        logger.info('protein {} will replace protein {} as lead, with index {}, New Num Peptides: {}, Old Num Peptides: {}'.format(str(new_lead.identifier),str(old_lead.identifier),str(j),str(len(new_lead.peptides)),str(len(old_lead.peptides))))
                        lead_protein_set.add(new_lead.identifier)
                        lead_protein_set.remove(old_lead.identifier)
                        # Swap their positions in the list
                        list_of_group_objects[i].proteins[0], list_of_group_objects[i].proteins[j] = new_lead, old_lead
                        break

        return(list_of_group_objects)


class Inclusion(Inference):

    def __init__(self,data_class,digest_class):
        self.data_class = data_class
        self.digest_class = digest_class
        self.scored_data = self.data_class.get_protein_data()
        self.data_class = data_class
        self.lead_protein_set = None

    def infer_proteins(self):

        logger = getLogger('protein_inference.inference.Inclusion.infer_proteins')

        group_dict = self._group_by_peptides(scored_data=self.scored_data, data_class=self.data_class,
                                             digest_class=self.digest_class, inference_type="inclusion", grouping_type=self.data_class.parameter_file_object.grouping_type)

        self.list_of_prots_not_in_db = group_dict["missing_proteins"]
        self.list_of_peps_not_in_db = group_dict["missing_peptides"]
        grouped_proteins = group_dict["grouped_proteins"]

        # Get the higher or lower variable
        if not self.data_class.high_low_better:
            hl = self.data_class.higher_or_lower()
            higher_or_lower = self.data_class.high_low_better
        else:
            higher_or_lower = self.data_class.high_low_better

        logger.info("Applying Group ID's for the Inclusion Method")
        # regrouped_proteins = self._apply_protein_group_ids(grouped_protein_objects = grouped_proteins,
        #                                                    data_class = self.data_class, digest_class = self.digest_class)

        regrouped_proteins = self._apply_group_ids_no_groups(grouped_protein_objects = grouped_proteins,
                                                           data_class = self.data_class, digest_class = self.digest_class)



        scores_grouped = regrouped_proteins["scores_grouped"]
        list_of_group_objects = regrouped_proteins["group_objects"]

        logger.info('Sorting Results based on lead Protein Score')
        if higher_or_lower == 'lower':
            scores_grouped = sorted(scores_grouped, key=lambda k: (float(k[0].score), -float(k[0].num_peptides)), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: (float(k.proteins[0].score), -float(k.proteins[0].num_peptides)),
                                           reverse=False)
        if higher_or_lower == 'higher':
            scores_grouped = sorted(scores_grouped, key=lambda k: (float(k[0].score),float(k[0].num_peptides)), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: (float(k.proteins[0].score),float(k.proteins[0].num_peptides)),
                                           reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped
        self.data_class.protein_group_objects = list_of_group_objects


class Exclusion(Inference):

    def __init__(self, data_class, digest_class):
        self.data_class = data_class
        self.digest_class = digest_class
        self.scored_data = self.data_class.get_protein_data()
        self.data_class = data_class
        self.lead_protein_set = None

    def infer_proteins(self):
        logger = getLogger('protein_inference.inference.Exclusion.infer_proteins')

        group_dict = self._group_by_peptides(scored_data=self.scored_data, data_class=self.data_class,
                                             digest_class=self.digest_class, inference_type="exclusion", grouping_type=self.data_class.parameter_file_object.grouping_type)

        self.list_of_prots_not_in_db = group_dict["missing_proteins"]
        self.list_of_peps_not_in_db = group_dict["missing_peptides"]
        grouped_proteins = group_dict["grouped_proteins"]

        # Get the higher or lower variable
        if not self.data_class.high_low_better:
            self.data_class.higher_or_lower()
        else:
            higher_or_lower = self.data_class.high_low_better

        logger.info("Applying Group ID's for the Exclusion Method")
        regrouped_proteins = self._apply_protein_group_ids(grouped_protein_objects=grouped_proteins,
                                                           data_class=self.data_class, digest_class=self.digest_class)

        scores_grouped = regrouped_proteins["scores_grouped"]
        list_of_group_objects = regrouped_proteins["group_objects"]

        logger.info('Sorting Results based on lead Protein Score')
        if higher_or_lower == 'lower':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=False)
        if higher_or_lower == 'higher':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped
        self.data_class.protein_group_objects = list_of_group_objects


class Parsimony(Inference):

    def __init__(self,data_class,digest_class):
        self.data_class = data_class
        self.digest_class = digest_class
        self.scored_data = self.data_class.get_protein_data()
        self.lead_protein_set = None
        self.parameter_file_object = data_class.parameter_file_object

    def _setup_glpk(self, glpkin_filename='glpkin.mod'):
        """
        This class is used to setup the glpk file for analysis.

        Example: protein_inference.grouping.GlpkSetup(data_class = data ,glpkin_filename='glpkin_example.mod'))

        Class will use attributes from DataStore object data.
        Class will also write the glpkin filename to be used in protein_inference.grouping.GlpkRunner()

        The Bulk of the glpk input file looks as follows:
        s.t. c1: y[5658] >=1;
        s.t. c2: y[14145]+y[4857]+y[4858]+y[10143]+y[2966] >=1;
        s.t. c3: y[320]+y[4893]+y[4209]+y[911]+y[2767]+y[2296]+y[10678]+y[3545] >=1
        """

        # Here we get the peptide to protein dictionary
        pep_prot_dict = self.data_class.peptide_to_protein_dictionary()

        # Here we get the protein to peptide dictionary...
        prot_pep_dict = self.data_class.protein_to_peptide_dictionary()

        identifiers_sorted = self.data_class.get_sorted_identifiers(digest_class=self.digest_class, scored=True)

        # Get all the proteins that we scored and the ones picked if picker was ran...
        data_proteins = [x for x in self.data_class.protein_peptide_dictionary.keys() if x in identifiers_sorted]
        # Get the set of peptides for each protein...
        data_peptides = [set(self.data_class.protein_peptide_dictionary[x]) for x in data_proteins]
        flat_peptides_in_data = set([item for sublist in data_peptides for item in sublist])
        peptide_sets = []
        # Loop over the list of peptides...
        for k in range(len(data_peptides)):
            raw_peptides = data_peptides[k]
            t_set = set()
            # Loop over each individual peptide per protein...
            for peps in raw_peptides:
                peptide = peps

                # Remove mods...
                new_peptide = Psm.remove_peptide_mods(peptide)
                # Add it to a temporary set...
                t_set.add(new_peptide)
            # Append this set to a new list...
            peptide_sets.append(t_set)
            # Set that proteins peptides to be the unmodified ones...
            data_peptides[k] = t_set

        # Get them all...
        all_peptides = [x for x in data_peptides]
        # Remove redundant sets...
        restrictedlst = [set(i) for i in OrderedDict.fromkeys(frozenset(item) for item in peptide_sets)]

        # Loop over  the restricted list of peptides...
        ind_list = []
        for p in restrictedlst:
            # Get its index in terms of the overall list...
            ind_list.append(all_peptides.index(p))

        # Get the protein based on the index
        restricted_proteins = [data_proteins[x] for x in range(len(data_peptides)) if x in ind_list]

        # Here we get the list of all proteins
        plist = []
        for peps in pep_prot_dict.keys():
            for prots in list(pep_prot_dict[peps]):
                if prots in restricted_proteins and peps in flat_peptides_in_data:
                    plist.append(prots)

        # Here we get the unique proteins
        unique_prots = list(set(plist).union())
        unique_protein_set = set(unique_prots)

        unique_prots_sorted = [x for x in identifiers_sorted if x in unique_prots]

        if len(unique_prots)!=len(unique_prots_sorted):
            raise ValueError("Sorted proteins length is not equal to unsorted length...")

        # Here we get a subset of the unique proteins
        # Where the peptides from the protein cannot be a subset of the peptides from any other protein that has already been added to the list...
        # unique_prots_restricted = []
        # unique_prots_peptides_restricted = []
        # For p in range(len(unique_prots))
        #   # At Index 0, simply append the protein and peptides...
        #   if p==0:
        #       unique_prots_restricted.append(unique_prots[p])
        #       unique_prots_peptides_restricted.append(pep_prot_dict[unique_prots[p]])
        #   else:
        #       current_peptides = pep_prot_dict[unique_prots[p]]
        #       for s in unique_prots_peptides_restricted:
        #           if current_peptides.issubset(s):
        #               issubset = True
        #               break
        #       if issubset==True:
        #           pass
        #       if issubset==False:
        #           unique_prots_restricted.append(unique_prots[p])
        #           unique_prots_peptides_restricted.append(pep_prot_dict[unique_prots[p]])
        #
        # unique_prots = unique_prots_restricted

        # Setup default dictionaries
        dd_num = collections.defaultdict(list)
        dd_prot_nums = collections.defaultdict(list)

        # For all the unique proteins from the search create a number to protein dictionary and a protein to number dictionary
        # Here we essentially assign a number to each protein
        # This is important as the glpk analysis just treats proteins as numbers...
        for p in range(len(unique_prots_sorted)):
            dd_num[unique_prots_sorted[p]].append(p)
            dd_prot_nums[p].append(unique_prots_sorted[p])

        # Store this data as glpk_protein_number_dictionary and glpk_number_protein_dictionary
        # The numbers are important as they are used in the GLPK input and we need to know what number in the GLPK output corresponds with which protein from the search
        self.data_class.glpk_protein_number_dictionary = dd_num
        self.data_class.glpk_number_protein_dictionary = dd_prot_nums
        # Create the GLPK input file
        fileout = open(glpkin_filename, 'w')

        # Not sure if this header string is correct or if it needs to be here...
        fileout.write('/* sets */' + '\n' + 'set PROTEINS;' + '\n' + '\n' + '\n')
        fileout.write('/* decision variables: yi, i in {1,..,5}. yi = 1 -> protein i is selected */' + '\n')
        fileout.write('var y {i in PROTEINS} binary >=0;' + '\n')
        fileout.write('/* objective function */' + '\n')
        fileout.write('minimize z: sum{i in PROTEINS} y[i];' + '\n' + '\n')
        fileout.write('/* Constraints */' + '\n')

        # Here we create the bulk of the input file which needs to look as follows:
        # s.t. c1: y[5658] >=1;
        # s.t. c2: y[14145]+y[4857]+y[4858]+y[10143]+y[2966] >=1;
        # s.t. c3: y[320]+y[4893]+y[4209]+y[911]+y[2767]+y[2296]+y[10678]+y[3545] >=1;
        # Each of the lines (constants, c1,c2,c3) is a peptide and each of the y[x] is a protein
        tot_peps = sorted(list(set(flat_peptides_in_data))) # Sort peptides alphabetically first...
        for j in range(len(tot_peps)):
            combine = ['y[' + str(dd_num[x][0]) + ']' for x in sorted(pep_prot_dict[tot_peps[j]]) if
                       x in unique_protein_set]
            fileout.write('s.t. c' + str(j + 1) + ': ' + '+'.join(combine) + ' >=1;' + '\n')

        # Finish writing the rest of the file and close it
        fileout.write('\n')
        fileout.write('data;' + '\n')
        numlist = [str(dd_num[x][0]) for x in sorted(unique_prots)]
        strlist = ' '.join(numlist)
        # End the file with listing the entire set of proteins... (as its number identifier)
        fileout.write('set PROTEINS := ' + strlist + ' ;' + '\n' + '\n')

        fileout.write('end;')
        fileout.close()

    def _glpk_runner(self, path_to_glpsol = 'glpsol',glpkin='glpkin.mod',glpkout='glpkout.sol', skip_running = False):
        """
        The GlpkRunner class takes a path to glpsol, the glpk input file from protein_inference.grouping.GlpkSetup(), a glpkout filename as well as a file_override option

        Example: protein_inference.grouping.GlpkRunner(path_to_glpsol = '/glpsol',glpkin='glpkin_example.mod',glpkout='glpkout_example.sol',file_override=False)

        path to glpsol on rescomp3 is: '/gne/research/apps/protchem/glpk/bin/glpsol'

        Typically set file_override to false unless you know what you are doing (IE you have a specific glpk solution file you want to use)

        Important output of this class is the glpk output solution file to be used in protein_inference.grouping.GlpkGrouper
        """
        logger = getLogger('protein_inference.inference.Parsimony._glpk_runner')

        #If there is no file_override (mainly for offline testing)
        #Run GLPK with the following command
        if not skip_running:
            p = subprocess.Popen(str(path_to_glpsol)+' -m '+str(glpkin)+' -o '+str(glpkout), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output = p.communicate()

            logger.info('Start Command line Stdout')
            logger.info(output[0])
            logger.info('End Command line Stdout')
            logger.info('Start Command line Stderr')
            logger.info(output[1])
            logger.info('End Command line Stderr')

            if output[0]=='':
                raise ValueError('Glpk did not produce any output... See potential error output above')
        else:
            logger.info("Not running GLPK, File {} will be used downstream in grouping".format(glpkout))




        ###Define indicies better and make this code more readable...
        ###Define indicies in init and show commented examples of how the data looks...

    def _glpk_grouper(self, data_class, digest_class, swissprot_override="soft", glpksolution_filename='glpkout.sol'):
        """
        This function takes a digest class object, a glpk solution file, as well as an option for swissprot override (protein naming convention).

        Example: protein_inference.parsimony._glpk_grouper(data_class = data,digest_class = digest,swissprot_override="soft",glpksolution_filename='glpkout_example.sol')

        Where data is a DataStore object, digest is a Digest.insilicodigest.InSilicoDigest() class object and glpksolution_filename is a glpk solution file.

        Finally, swissprot_override is a lead protein override for naming convention which can have 3 options: False, "soft", or "hard".
        selecting False skips the override.

        selecting "soft" will cycle through all protein leads and groups. If a protein lead is an unreviewed hit (trembl) then all proteins in the lead proteins group are inspected.
        If any of these intra group proteins are a reviewed hit (swissprot) and the lead trembl proteins peptides are a complete subset of the swissprot proteins peptides, then we select
        to swap the trembl and swissprot hits so that the swissprot is now the new lead protein.

        We opt to use this soft override scheme as default because given the redundancy of our databases. Out of GLPK we see a lot of unreviewed hits being called lead proteins.
        If proteins share the same peptides or share a very similar set of peptides we are unsure how GLPK selects which protein to supply as the lead.
        As such, we get many reviewed and unreviewed proteins as lead proteins.
        Performing this "soft" override switches many of these lead trembl hits to reviewed swissprot hits.
        This is important as group members are used to seeing swissprot identifiers as opposed to trembl identifiers

        selecting "hard" will perform the same swapping as the "soft" override. However,
        the "hard" override will swap a trembl lead with the swissprot protein with the highest number of peptides in its group (so long as its already not a lead protein itself)
        even if the unreviewed has a unique peptide relative to the peptides of all proteins in its group.
        This setting is not recommended given that you can potentially lose out on Unique peptides

        """
        logger = getLogger('protein_inference.inference.Parsimony._glpk_grouper')

        scored_data = data_class.get_protein_data()


        glpk_out = open(glpksolution_filename, 'r')

        # Get the number protein dictionary from glpk_setup
        dd_prot_nums = self.data_class.glpk_number_protein_dictionary

        glpk_out = glpk_out.read()
        glpk_out = glpk_out.split('\n')

        # Cant find a better way to do this... there are modules out there that work with glpk...
        start = glpk_out.index('   No. Column name       Activity     Lower bound   Upper bound')

        newlist = []
        # Fix this -13 and +2... not really sure how
        # Based on the output file we should start two lines after the start index and stop 13 before the end of the file...
        for lines in range(start + 2, len(glpk_out) - 13):
            new = [x.strip() for x in glpk_out[lines].split(' ')]
            res = []
            for stuff in new:
                if stuff != '':
                    res.append(stuff)
            newlist.append(res)

        # Use 1 here as 1 is the location in the line of the protein number
        # 3 is the location of the binary (which indicates whether or not the protein number has a unique peptide making it a lead protein)
        numbers = [x[1].split(']')[0].split('[')[-1] for x in newlist]
        binary = [x[3] for x in newlist]

        self.numbers = numbers

        # Here we extract out the lead proteins...
        lead_proteins = []
        for k in range(len(numbers)):
            if binary[k] == '1':
                try:
                    passing_protein_number = int(numbers[k])
                    lead_proteins.append(dd_prot_nums[passing_protein_number][0])
                except IndexError:
                    logger.warning("No Protein for Protein Number" + str(passing_protein_number))

        lead_protein_set = set(lead_proteins)
        self.lead_protein_set = lead_protein_set

        logger.info('Number of lead proteins = ' + str(len(lead_proteins)))

        scored_proteins = list(scored_data)
        protein_finder = [x.identifier for x in scored_proteins]

        lead_protein_objects = []
        lead_protein_identifiers = []
        for proteins in lead_proteins:
            if proteins in protein_finder:
                p_ind = protein_finder.index(proteins)
                protein_object = scored_proteins[p_ind]
                lead_protein_objects.append(protein_object)
                lead_protein_identifiers.append(protein_object.identifier)
            else:
                # Why are some proteins not being found when we run exclusion???
                logger.warning("Protein {} not found with protein finder...".format(proteins))

        self.lead_protein_objects = lead_protein_objects

        # Now we have the lead Proteins so we need to get the peptides for each lead protein
        # Then we find all proteins that share at least 1 peptide with each lead protein
        # If they share at least 1 peptide then assign that protein to the group...
        self.data_class.glpk_lead_proteins = lead_protein_objects

        group_dict = self._group_by_peptides(scored_data=scored_data, data_class=self.data_class,
                                             digest_class=self.digest_class, inference_type="parsimony",
                                             lead_protein_objects=self.lead_protein_objects, grouping_type=self.data_class.parameter_file_object.grouping_type)

        self.list_of_prots_not_in_db = group_dict["missing_proteins"]
        self.list_of_peps_not_in_db = group_dict["missing_peptides"]
        grouped_proteins = group_dict["grouped_proteins"]

        regrouped_proteins = self._swissprot_and_isoform_override(scored_data=scored_data,
                                                                  grouped_proteins=grouped_proteins,
                                                                  data_class=self.data_class,
                                                                  digest_class=self.digest_class,
                                                                  override_type="soft", isoform_override=True)

        scores_grouped = regrouped_proteins["scores_grouped"]
        list_of_group_objects = regrouped_proteins["group_objects"]
        lead_replaced_prot_pairs = regrouped_proteins["lead_replaced"]

        self.data_class.lead_replaced_proteins = lead_replaced_prot_pairs

        # Get the higher or lower variable
        if not self.data_class.high_low_better:
            higher_or_lower = self.data_class.higher_or_lower()
        else:
            higher_or_lower = self.data_class.high_low_better

        logger.info('Sorting Results based on lead Protein Score')
        if higher_or_lower == 'lower':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=False)
        if higher_or_lower == 'higher':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=True)

        list_of_group_objects = self._reassign_leads(list_of_group_objects, data_class=self.data_class)

        logger.info('Re Sorting Results based on lead Protein Score')
        if higher_or_lower == 'lower':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=False)
        if higher_or_lower == 'higher':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped
        self.data_class.protein_group_objects = list_of_group_objects

    def infer_proteins(self, glpkinout_directory = "glpkinout/", skip_running_glpk = False):

        logger = getLogger('protein_inference.inference.Parsimony.infer_proteins')

        try:
            os.mkdir(glpkinout_directory)
        except OSError:
            logger.warning("Directory {} cannot be created or already exists".format(glpkinout_directory))

        self._setup_glpk(glpkin_filename=os.path.join(glpkinout_directory,'glpkin_' + self.parameter_file_object.tag + '.mod'))

        # For reference server_glpk_path = '/gne/research/apps/protchem/glpk/bin/glpsol'
        self._glpk_runner(path_to_glpsol=self.parameter_file_object.glpk_path,
                          glpkin=os.path.join(glpkinout_directory,'glpkin_' + self.parameter_file_object.tag + '.mod'),
                          glpkout=os.path.join(glpkinout_directory,'glpkout_' + self.parameter_file_object.tag + '.sol'),
                          skip_running=skip_running_glpk)

        self._glpk_grouper(data_class=self.data_class, digest_class=self.digest_class,
                           swissprot_override='soft',
                           glpksolution_filename=os.path.join(glpkinout_directory,'glpkout_' + self.parameter_file_object.tag + '.sol'))

class FirstProtein(Inference):

    def __init__(self, data_class, digest_class):
        self.data_class = data_class
        self.digest_class = digest_class
        self.scored_data = self.data_class.get_protein_data()
        self.data_class = data_class
        self.lead_protein_set = None

    def infer_proteins(self):

        logger = getLogger('protein_inference.inference.FirstProtein.infer_proteins')

        group_dict = self._group_by_peptides(scored_data=self.scored_data, data_class=self.data_class,
                                             digest_class=self.digest_class, inference_type="none", grouping_type=self.data_class.parameter_file_object.grouping_type)

        self.list_of_prots_not_in_db = group_dict["missing_proteins"]
        self.list_of_peps_not_in_db = group_dict["missing_peptides"]
        grouped_proteins = group_dict["grouped_proteins"]

        # Get the higher or lower variable
        if not self.data_class.high_low_better:
            self.data_class.higher_or_lower()
        else:
            higher_or_lower = self.data_class.high_low_better

        logger.info("Applying Group ID's for the First Protein Method")
        regrouped_proteins = self._apply_protein_group_ids(grouped_protein_objects=grouped_proteins,
                                                           data_class=self.data_class, digest_class=self.digest_class)

        scores_grouped = regrouped_proteins["scores_grouped"]
        list_of_group_objects = regrouped_proteins["group_objects"]

        logger.info('Sorting Results based on lead Protein Score')
        if higher_or_lower == 'lower':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=False)
        if higher_or_lower == 'higher':
            scores_grouped = sorted(scores_grouped, key=lambda k: float(k[0].score), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score),
                                           reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped
        self.data_class.protein_group_objects = list_of_group_objects


class PeptideCentric(Inference):

    def __init__(self, data_class, digest_class):

        self.data_class = data_class
        self.digest_class = digest_class
        self.scored_data = self.data_class.get_protein_data()
        self.lead_protein_set = None
        self.parameter_file_object = data_class.parameter_file_object


    def infer_proteins(self):

        logger = getLogger('protein_inference.inference.PeptideCentric.infer_proteins')

        # Get the higher or lower variable
        if not self.data_class.high_low_better:
            hl = self.data_class.higher_or_lower()
            higher_or_lower = self.data_class.high_low_better
        else:
            higher_or_lower = self.data_class.high_low_better

        logger.info("Applying Group ID's for the Peptide Centric Method")

        regrouped_proteins = self._apply_group_ids_peptide_centric(data_class=self.data_class)

        scores_grouped = regrouped_proteins["scores_grouped"]
        list_of_group_objects = regrouped_proteins["group_objects"]

        logger.info('Sorting Results based on lead Protein Score')
        if higher_or_lower == 'lower':
            scores_grouped = sorted(scores_grouped, key=lambda k: (float(k[0].score), -float(k[0].num_peptides)),
                                    reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: (
            float(k.proteins[0].score), -float(k.proteins[0].num_peptides)),
                                           reverse=False)
        if higher_or_lower == 'higher':
            scores_grouped = sorted(scores_grouped, key=lambda k: (float(k[0].score), float(k[0].num_peptides)),
                                    reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: (
            float(k.proteins[0].score), float(k.proteins[0].num_peptides)),
                                           reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped
        self.data_class.protein_group_objects = list_of_group_objects