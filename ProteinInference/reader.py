#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 12:52:14 2017

@author: hinklet
"""

from physical import Psm
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
    The following class takes a percolator target file and a percolator decoy file and creates standard ProteinInference.physical.Psm() objects.
    This reader class is used as input for the DataStore class which is used in all further ProteinInference Classes

    Example: ProteinInference.reader.PercolatorRead(target_file = "example_target.txt", decoy_file = "example_decoy.txt")

    Percolator Output is formatted as follows:
    with each entry being tabbed delimited
    PSMId	 score	q-value	posterior_error_prob	peptide	proteinIds
    116108.15139.15139.6.dta	 3.44016	 0.000479928	7.60258e-10	K.MVVSMTLGLHPWIANIDDTQYLAAK.R	CNDP1_HUMAN|Q96KN2	B4E180_HUMAN|B4E180	A8K1K1_HUMAN|A8K1K1	J3KRP0_HUMAN|J3KRP0

    """
    def __init__(self,target_file,decoy_file):
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
        
    def execute(self):
        #Read in and split by line
        t = open(self.target_file)
        t = t.read()
        lines = t.split('\n')
        
        #Split each line by tab
        ptarg = [x.split('\t') for x in lines]
        #Remove empty ending line in file
        del ptarg[-1]
        #Remove header
        del ptarg[0]
        
        #Repeat for decoy file
        d = open(self.decoy_file)
        d = d.read()
        dlines = d.split('\n')
        
        pdec = [x.split('\t') for x in dlines]
        #Remove empty ending in file
        del pdec[-1]
        #Remove header
        del pdec[0]
        
        #Combine the lists
        perc_all = ptarg+pdec
        perc_all = sorted(perc_all, key=lambda x: float(x[1]), reverse = True)
        list_of_psm_objects = []
        peptide_tracker = []
        #We only want to get unique peptides... using all messes up scoring...
        #Create Psm objects with the identifier, percscore, qvalue, pepvalue, and possible proteins...
        for psm_info in perc_all:
            current_peptide = psm_info[self.peptide_index]
            #Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier = current_peptide)
                #Add all the attributes
                p.percscore = float(psm_info[self.perc_score_index])
                p.qvalue = float(psm_info[self.q_value_index])
                p.pepvalue = float(psm_info[self.posterior_error_prob_index])
                p.possible_proteins = psm_info[self.proteinIDs_index:]
                list_of_psm_objects.append(p)
                peptide_tracker.append(current_peptide)

        self.psms = list_of_psm_objects

        
        #return perc
        
class OtherMethod(Reader):
    """
    Potential Future class to read in from another source.
    Potentially from a database, another Psm scoring source, or potentially from Percolator XML source.
    This Method will be used to read from post processing proteologic logical object... which can either be LDA or Percolator Results...
    Essentially it will just be a searchID and an LDA/Percolator ID
    """
    
    def __init__(self,something_else):
        self.something_else=something_else
        
    def print_stuff(self):
        print self.something_else