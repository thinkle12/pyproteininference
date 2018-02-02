#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:43:35 2017

@author: hinklet
"""



class Protein(object):
    __slots__ = ['identifier','score', 'group_identification', 'reviewed', 'unreviewed', 'peptides', 'peptide_scores', 'picked', 'psm_score_dictionary',
                 'num_peptides','unique_peptides','num_unique_peptides','raw_peptides']
    
    def __init__(self,identifier):
        self.identifier = identifier
        self.score = None
        self.group_identification = []
        self.reviewed = False
        self.unreviewed = False
        self.peptides = None
        self.peptide_scores = None
        self.picked = True
        self.psm_score_dictionary = None
        self.num_peptides = None
        self.unique_peptides = None
        self.num_unique_peptides = None
        self.raw_peptides = None
        
        
class Psm(object):
    __slots__ = ['identifier','percscore','qvalue','pepvalue','possible_proteins']
    
    def __init__(self,identifier):
        self.identifier = identifier
        self.percscore = None
        self.qvalue = None
        self.pepvalue = None
        self.possible_proteins = None
        

class ProteinGroup(object):
    __slots__ = ['proteins','number_id','q_value']

    def __init__(self,number_id):
        self.proteins = []
        self.number_id = number_id
        self.q_value = None
