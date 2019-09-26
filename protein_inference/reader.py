#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 12:52:14 2017

@author: hinklet
"""

from protein_inference.physical import Psm
import yaml
import csv
import re
import itertools
##Here we create the reader class which is capable of reading in files as input for ProteinInference
##The reader class creates a list of psm objects with various attributes...

class Reader(object):
    # TODO make this an abstract base class... IF I can
    """
    Main Reader Class which is parent to all reader subclasses
    """
    
    def __init__(self):
        None

    def remap(self, fieldnames):
        price_count = itertools.count(1)
        return ['alternative_protein_{}'.format(next(price_count)) if f=="" else f
                for f in fieldnames]
        
class PercolatorReader(Reader):
    """
    The following class takes a percolator target file and a percolator decoy file and creates standard protein_inference.physical.Psm() objects.
    This reader class is used as input for the DataStore class which is used in all further protein_inference Classes

    Example: protein_inference.reader.PercolatorReader(target_file = "example_target.txt", decoy_file = "example_decoy.txt")

    Percolator Output is formatted as follows:
    with each entry being tabbed delimited
    PSMId	 score	q-value	posterior_error_prob	peptide	proteinIds
    116108.15139.15139.6.dta	 3.44016	 0.000479928	7.60258e-10	K.MVVSMTLGLHPWIANIDDTQYLAAK.R	CNDP1_HUMAN|Q96KN2	B4E180_HUMAN|B4E180	A8K1K1_HUMAN|A8K1K1	J3KRP0_HUMAN|J3KRP0

    """
    def __init__(self,target_file,decoy_file,digest_class,parameter_file_object):
        self.target_file = target_file
        self.decoy_file = decoy_file
        #Define Indicies based on input
        # TODO Change this from index based to reading in as a dict... similar to percolator
        self.psmid_index = 0
        self.perc_score_index = 1
        self.q_value_index = 2
        self.posterior_error_prob_index = 3
        self.peptide_index = 4
        self.proteinIDs_index = 5
        self.psms = None
        self.search_id = None
        self.digest_class = digest_class

        self.parameter_file_object = parameter_file_object

        
    def read_psms(self):
        #Read in and split by line
        # If target_file is a list... read them all in and concatenate...
        if isinstance(self.target_file, (list,)):
            all_target = []
            for t_files in self.target_file:
                print(t_files)
                ptarg = []
                with open(t_files, 'r') as perc_target_file:
                    spamreader = csv.reader(perc_target_file, delimiter='\t')
                    for row in spamreader:
                        ptarg.append(row)
                del ptarg[0]
                all_target = all_target + ptarg
        else:
            # If not just read the file...
            ptarg = []
            with open(self.target_file, 'r') as perc_target_file:
                spamreader = csv.reader(perc_target_file, delimiter='\t')
                for row in spamreader:
                    ptarg.append(row)
            del ptarg[0]
            all_target = ptarg

        
        #Repeat for decoy file
        if isinstance(self.decoy_file, (list,)):
            all_decoy = []
            for d_files in self.decoy_file:
                print(d_files)
                pdec = []
                with open(d_files, 'r') as perc_decoy_file:
                    spamreader = csv.reader(perc_decoy_file, delimiter='\t')
                    for row in spamreader:
                        pdec.append(row)
                del pdec[0]
                all_decoy = all_decoy + pdec
        else:
            pdec = []
            with open(self.decoy_file, 'r') as perc_decoy_file:
                spamreader = csv.reader(perc_decoy_file, delimiter='\t')
                for row in spamreader:
                    pdec.append(row)
            del pdec[0]
            all_decoy = pdec

        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary
        #Combine the lists
        perc_all = all_target+all_decoy

        perc_all_filtered = []
        for psms in perc_all:
            try:
                float(psms[self.posterior_error_prob_index])
                perc_all_filtered.append(psms)
            except ValueError as e:
                pass

        # Filter by pep
        perc_all = sorted(perc_all_filtered, key=lambda x: float(x[self.posterior_error_prob_index]), reverse = False)

        # TODO
        # TRY TO GET PERC_ALL AS A GENERATOR
        # Can do this... just give the option to feed a combined file... and loop over all the files...
        # Hmm... problem is it still needs to be sorted by perc score....

        list_of_psm_objects = []
        peptide_tracker = set()
        all_sp_proteins = set(self.digest_class.swiss_prot_protein_set)
        #We only want to get unique peptides... using all messes up scoring...
        #Create Psm objects with the identifier, percscore, qvalue, pepvalue, and possible proteins...

        # TODO
        # make this for loop a generator...

        print(len(perc_all))
        for psm_info in perc_all:
            current_peptide = psm_info[self.peptide_index]
            #Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier = current_peptide)
                #Add all the attributes
                p.percscore = float(psm_info[self.perc_score_index])
                p.qvalue = float(psm_info[self.q_value_index])
                p.pepvalue = float(psm_info[self.posterior_error_prob_index])
                poss_proteins = list(set(psm_info[self.proteinIDs_index:self.proteinIDs_index+50]))
                p.possible_proteins = poss_proteins # Restrict to 50 total possible proteins...
                p.psm_id = psm_info[self.psmid_index]

                if "." in current_peptide:
                    # If we have full peptides split it and take the middle...
                    current_peptide = current_peptide.split(".")[1]
                if not current_peptide.isupper() or not current_peptide.isalpha():
                    # If we have mods remove them...
                    peptide_string = current_peptide.upper()
                    stripped_peptide = Psm.remove_peptide_mods(peptide_string)
                    current_peptide = stripped_peptide
                # Add the other possible_proteins from insilicodigest here...
                current_alt_proteins = list(peptide_to_protein_dictionary[current_peptide]) # TODO This peptide needs to be scrubbed of Mods...
                # Sort Alt Proteins by Swissprot then Trembl...
                our_target_sp_proteins = sorted(
                    [x for x in current_alt_proteins if x in all_sp_proteins and self.parameter_file_object.decoy_symbol not in x])
                our_decoy_sp_proteins = sorted(
                    [x for x in current_alt_proteins if x in all_sp_proteins and self.parameter_file_object.decoy_symbol in x])

                our_target_tr_proteins = sorted(
                    [x for x in current_alt_proteins if x not in all_sp_proteins and self.parameter_file_object.decoy_symbol not in x])
                our_decoy_tr_proteins = sorted(
                    [x for x in current_alt_proteins if x not in all_sp_proteins and self.parameter_file_object.decoy_symbol in x])

                identifiers_sorted = our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins

                # Restrict to 50 possible proteins
                for alt_proteins in identifiers_sorted[:50]:
                    if alt_proteins not in p.possible_proteins:
                        p.possible_proteins.append(alt_proteins)

                p.possible_proteins = [x for x in p.possible_proteins if x]
                if not current_alt_proteins:
                    print("Peptide {} was not found in the supplied DB".format(current_peptide))




                list_of_psm_objects.append(p)
                peptide_tracker.add(current_peptide)


        self.psms = list_of_psm_objects

        print(len(self.psms))

        
        #return perc
        
class ProteologicPostSearchReader(Reader):
    """
    Potential Future class to read in from another source.
    Potentially from a database, another Psm scoring source, or potentially from Percolator XML source.
    This Method will be used to read from post processing proteologic logical object... which can either be LDA or Percolator Results...
    Essentially it will just be a searchID and an LDA/Percolator ID
    """
    
    def __init__(self, proteologic_object, search_id, postsearch_id, digest_class, parameter_file_object):
        self.proteologic_object=proteologic_object
        self.search_id = search_id
        self.postsearch_id = postsearch_id

        self.psms = None
        self.digest_class = digest_class

        self.parameter_file_object = parameter_file_object




    def read_psms(self):
        print('Reading in data...')
        if isinstance(self.proteologic_object, (list,)):
            list_of_psms = []
            for p_objs in self.proteologic_object:
                for psms in p_objs.physical_object.psm_sets:
                    list_of_psms.append(psms)
        else:
            list_of_psms = self.proteologic_object.physical_object.psm_sets

        # Sort this by posterior error prob...
        list_of_psms = sorted(list_of_psms, key=lambda x: float(x.psm_filter.pepvalue))




        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary

        list_of_psm_objects = []
        peptide_tracker = set()
        all_sp_proteins = set(self.digest_class.swiss_prot_protein_set)
        #Peptide tracker is used because we only want UNIQUE peptides...
        #The data is sorted by percolator score... or at least it should be...
        #Or sorted by posterior error probability

        # TODO, We may need to try the quant thing on ALL PSMs, not just unique Peptides...
        for peps in list_of_psms:
            self.peps = peps
            current_peptide = peps.peptide.sequence
            # Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier=current_peptide)
                # Add all the attributes
                p.percscore = float(0) # Will be stored in table in future I think...
                p.qvalue = float(peps.psm_filter.q_value)
                p.pepvalue = float(peps.psm_filter.pepvalue)
                if peps.peptide.protein not in peps.alternative_proteins:
                    p.possible_proteins = [peps.peptide.protein] + peps.alternative_proteins
                else:
                    p.possible_proteins = peps.alternative_proteins

                p.possible_proteins = list(filter(None, p.possible_proteins))
                p.psm_id = peps.spectrum.spectrum_identifier

                if "." in current_peptide:
                    # If we have full peptides split it and take the middle...
                    current_peptide = current_peptide.split(".")[1]
                if not current_peptide.isupper() or not current_peptide.isalpha():
                    # If we have mods remove them...
                    peptide_string = current_peptide.upper()
                    stripped_peptide = Psm.remove_peptide_mods(peptide_string)
                    current_peptide = stripped_peptide
                # Add the other possible_proteins from insilicodigest here...
                current_alt_proteins = list(peptide_to_protein_dictionary[current_peptide])
                # Sort Alt Proteins by Swissprot then Trembl...
                our_target_sp_proteins = sorted(
                    [x for x in current_alt_proteins if x in all_sp_proteins and self.parameter_file_object.decoy_symbol not in x])
                our_decoy_sp_proteins = sorted(
                    [x for x in current_alt_proteins if x in all_sp_proteins and self.parameter_file_object.decoy_symbol in x])

                our_target_tr_proteins = sorted(
                    [x for x in current_alt_proteins if x not in all_sp_proteins and self.parameter_file_object.decoy_symbol not in x])
                our_decoy_tr_proteins = sorted(
                    [x for x in current_alt_proteins if x not in all_sp_proteins and self.parameter_file_object.decoy_symbol in x])

                identifiers_sorted = our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins

                # Restrict to 50 possible proteins...
                for alt_proteins in identifiers_sorted[:50]:
                    if alt_proteins not in p.possible_proteins:
                        p.possible_proteins.append(alt_proteins)
                if not current_alt_proteins:
                    print("Peptide {} was not found in the supplied DB".format(current_peptide))




                list_of_psm_objects.append(p)
                peptide_tracker.add(current_peptide)

            # TODO, Here we will Keep track of a large set of all PSMs..
            # We will keep track of this new PSMid -> search_id.peptide_id


        self.psms = list_of_psm_objects
        print('Finished reading in data...')


