#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:43:35 2017

@author: hinklet
"""



class Protein(object):
    """
    The following class is a physical Protein class that stores characteristics of a protein for the entire analysis.
    We use __slots__ to predefine the attributes the Protein Object can have
    This is done to speed up runtime of the PI algorithm

    Example: protein_inference.physical.Protein(identifier = "PRKDC_HUMAN|P78527")
    """
    __slots__ = ['identifier','score', 'group_identification', 'reviewed', 'unreviewed', 'peptides', 'peptide_scores', 'picked', 'psm_score_dictionary',
                 'num_peptides','unique_peptides','num_unique_peptides','raw_peptides','psmid_peptide_dictionary']
    
    def __init__(self,identifier):
        self.identifier = identifier
        self.score = None
        self.group_identification = set()
        self.reviewed = False
        self.unreviewed = False
        self.peptides = None
        self.peptide_scores = None
        self.picked = True
        self.psm_score_dictionary = None
        self.num_peptides = None
        self.unique_peptides = None
        self.num_unique_peptides = None
        self.raw_peptides = set()
        self.psmid_peptide_dictionary = None
        
        
class Psm(object):
    """
    The following class is a physical Psm class that stores characteristics of a psm for the entire analysis.
    We use __slots__ to predefine the attributes the Psm Object can have
    This is done to speed up runtime of the PI algorithm

    Example: Psm(identifier = "K.DLIDEGHAATQLVNQLHDVVVENNLSDK.Q")
    """
    __slots__ = ['identifier','percscore','qvalue','pepvalue','possible_proteins','psm_id']
    
    def __init__(self,identifier):
        self.identifier = identifier
        self.percscore = None
        self.qvalue = None
        self.pepvalue = None
        self.possible_proteins = None
        self.psm_id = None
        

class ProteinGroup(object):
    """
    The following class is a physical Protein Group class that stores characteristics of a Protein Group for the entire analysis.
    We use __slots__ to predefine the attributes the Psm Object can have
    This is done to speed up runtime of the PI algorithm

    Example: ProteinGroup(number_id = 1)
    """
    __slots__ = ['proteins','number_id','q_value']

    def __init__(self,number_id):
        self.proteins = []
        self.number_id = number_id
        self.q_value = None
