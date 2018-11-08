#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 12:52:14 2017

@author: hinklet
"""

from physical import Psm
import yaml
##Here we create the reader class which is capable of reading in files as input for ProteinInference
##The reader class creates a list of psm objects with various attributes...

class Reader(object):
    """
    Main Reader Class which is parent to all reader subclasses
    """
    
    def __init__(self):
        None 
        
class PercolatorRead(Reader):
    """
    The following class takes a percolator target file and a percolator decoy file and creates standard protein_inference.physical.Psm() objects.
    This reader class is used as input for the DataStore class which is used in all further protein_inference Classes

    Example: protein_inference.reader.PercolatorRead(target_file = "example_target.txt", decoy_file = "example_decoy.txt")

    Percolator Output is formatted as follows:
    with each entry being tabbed delimited
    PSMId	 score	q-value	posterior_error_prob	peptide	proteinIds
    116108.15139.15139.6.dta	 3.44016	 0.000479928	7.60258e-10	K.MVVSMTLGLHPWIANIDDTQYLAAK.R	CNDP1_HUMAN|Q96KN2	B4E180_HUMAN|B4E180	A8K1K1_HUMAN|A8K1K1	J3KRP0_HUMAN|J3KRP0

    """
    def __init__(self,target_file,decoy_file,digest_class,yaml_param_file):
        self.target_file = target_file
        self.decoy_file = decoy_file
        #Define Indicies based on input
        self.psmid_index = 0
        self.perc_score_index = 1
        self.q_value_index = 2
        self.posterior_error_prob_index = 3
        self.peptide_index = 4
        self.proteinIDs_index = 5
        self.psms = None
        self.search_id = None
        self.digest_class = digest_class

        self.yaml_param_file = yaml_param_file
        if self.yaml_param_file:
            with open(self.yaml_param_file, 'r') as stream:
                yaml_params = yaml.load(stream)
            #Read param file..... params read in as dict...
            #We will then have a bunch of... if param in param_dict.... do... if not pass....
        else:
            #Default params to be used if no Yaml file is provided...
            yaml_params = {'Parameters': {'Picker': True,
                                      'Restrict_Q': False,
                                      'Restrict_Pep': False,
                                      'Restrict_Peptide_Length': 7,
                                          'Group': "glpk",
                                          'Export': "q_value_leads",
                                          'FDR': .01,
                                          'Score_Type': "q_value",
                                          'Score_Method': "downweighted_multiplicative_log",
                                          'Missed_Cleavages': 2,
                                          'Digest_Type': "trypsin",
                                          'GLPK_Path':'glpsol'}}

        self.yaml_params = yaml_params

        
    def execute(self):
        #Read in and split by line
        # If target_file is a list... read them all in and concatenate...
        if isinstance(self.target_file, (list,)):
            all_target = []
            for t_files in self.target_file:
                print(t_files)
                t = open(t_files)
                t = t.read()
                lines = t.split('\n')

                #Split each line by tab
                ptarg = [x.split('\t') for x in lines]
                #Remove empty ending line in file
                del ptarg[-1]
                #Remove header
                del ptarg[0]
                all_target = all_target + ptarg
        else:
            # If not just read the file...
            t = open(self.target_file,)
            t = t.read()
            lines = t.split('\n')

            # Split each line by tab
            ptarg = [x.split('\t') for x in lines]
            # Remove empty ending line in file
            del ptarg[-1]
            # Remove header
            del ptarg[0]
            all_target = ptarg

        
        #Repeat for decoy file
        if isinstance(self.decoy_file, (list,)):
            all_decoy = []
            for d_files in self.decoy_file:
                print(d_files)
                d = open(d_files)
                d = d.read()
                dlines = d.split('\n')

                pdec = [x.split('\t') for x in dlines]
                #Remove empty ending in file
                del pdec[-1]
                #Remove header
                del pdec[0]
                all_decoy = all_decoy + pdec
        else:
            d = open(self.decoy_file)
            d = d.read()
            dlines = d.split('\n')

            pdec = [x.split('\t') for x in dlines]
            # Remove empty ending in file
            del pdec[-1]
            # Remove header
            del pdec[0]
            all_decoy = pdec

        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary
        #Combine the lists
        perc_all = all_target+all_decoy
        perc_all = sorted(perc_all, key=lambda x: float(x[1]), reverse = True)

        # TODO
        # TRY TO GET PERC_ALL AS A GENERATOR
        # Can do this... just give the option to feed a combined file... and loop over all the files...
        # Hmm... problem is it still needs to be sorted by perc score....

        list_of_psm_objects = []
        peptide_tracker = set()
        #We only want to get unique peptides... using all messes up scoring...
        #Create Psm objects with the identifier, percscore, qvalue, pepvalue, and possible proteins...

        # TODO
        # make this for loop a generator...

        print(len(perc_all))
        i = 0
        for psm_info in perc_all:
            i = i+1
            print(i)
            current_peptide = psm_info[self.peptide_index]
            #Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier = current_peptide)
                #Add all the attributes
                p.percscore = float(psm_info[self.perc_score_index])
                p.qvalue = float(psm_info[self.q_value_index])
                p.pepvalue = float(psm_info[self.posterior_error_prob_index])
                p.possible_proteins = psm_info[self.proteinIDs_index:]
                p.psm_id = psm_info[self.psmid_index]

                # Add the other possible_proteins from insilicodigest here...
                current_alt_proteins = list(peptide_to_protein_dictionary[current_peptide])
                for alt_proteins in current_alt_proteins:
                    if alt_proteins not in p.possible_proteins:
                        p.possible_proteins.append(alt_proteins)

                list_of_psm_objects.append(p)
                peptide_tracker.add(current_peptide)

        self.psms = list_of_psm_objects

        
        #return perc
        
class ProteologicPostSearch(Reader):
    """
    Potential Future class to read in from another source.
    Potentially from a database, another Psm scoring source, or potentially from Percolator XML source.
    This Method will be used to read from post processing proteologic logical object... which can either be LDA or Percolator Results...
    Essentially it will just be a searchID and an LDA/Percolator ID
    """
    
    def __init__(self,proteologic_object,search_id,percolator_id,digest_class,yaml_param_file):
        self.proteologic_object=proteologic_object
        self.search_id = search_id
        self.percolator_id = percolator_id

        self.psms = None
        self.digest_class = digest_class

        self.yaml_param_file = yaml_param_file
        if self.yaml_param_file:
            with open(self.yaml_param_file, 'r') as stream:
                yaml_params = yaml.load(stream)
            #Read param file..... params read in as dict...
            #We will then have a bunch of... if param in param_dict.... do... if not pass....
        else:
            #Default params to be used if no Yaml file is provided...
            yaml_params = {'Parameters': {'Picker': True,
                                      'Restrict_Q': False,
                                      'Restrict_Pep': False,
                                      'Restrict_Peptide_Length': 7,
                                          'Group': "glpk",
                                          'Export': "q_value_leads",
                                          'FDR': .01,
                                          'Score_Type': "q_value",
                                          'Score_Method': "downweighted_multiplicative_log",
                                          'Missed_Cleavages': 2,
                                          'Digest_Type': "trypsin",
                                          'GLPK_Path':'glpsol'}}

        self.yaml_params = yaml_params


    def execute(self):
        print('Reading in data...')
        #I want this to be logical object but it doesnt work... weird...
        list_of_psms = self.proteologic_object.physical_object.psm_sets
        #Sort this by posterior error prob...

        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary

        list_of_psm_objects = []
        peptide_tracker = []
        #Peptide tracker is used because we only want UNIQUE peptides...
        #The data is sorted by percolator score... or at least it should be...
        #Or sorted by posterior error probability
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

                # Add the other possible_proteins from insilicodigest here...
                current_alt_proteins = list(peptide_to_protein_dictionary[current_peptide])
                for alt_proteins in current_alt_proteins:
                    if alt_proteins not in p.possible_proteins:
                        p.possible_proteins.append(alt_proteins)


                list_of_psm_objects.append(p)
                peptide_tracker.append(current_peptide)

        self.psms = list_of_psm_objects
        print('Finished reading in data...')


