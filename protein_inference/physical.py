#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:43:35 2017

@author: hinklet
"""
import re


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
    # TODO need to be able to handle identifiers with and without trailing/leading AA with the .'s
    """
    __slots__ = ['identifier','percscore','qvalue','pepvalue','possible_proteins','psm_id','custom_score']

    # The regex removes anything between parantheses including parenthases - \([^()]*\)
    # The regex removes anything between brackets including parenthases - \[.*?\]
    # And the regex removes anything that is not an A-Z character [^A-Z]
    MOD_REGEX = re.compile('\([^()]*\)|\[.*?\]|[^A-Z]')

    FRONT_FLANKING_REGEX = re.compile('^[A-Z|-][.]')
    BACK_FLANKING_REGEX = re.compile('[.][A-Z|-]$')
    
    def __init__(self,identifier):
        self.identifier = identifier
        self.percscore = None
        self.qvalue = None
        self.pepvalue = None
        self.possible_proteins = None
        self.psm_id = None
        self.custom_score = None

    @classmethod
    def remove_peptide_mods(cls, peptide_string):
        stripped_peptide = cls.MOD_REGEX.sub('', peptide_string)
        return(stripped_peptide)

    @classmethod
    def split_peptide(cls, peptide_string, delimiter = "."):
        peptide_split = peptide_string.split(delimiter)
        if len(peptide_split)==3:
            # If we get 3 chunks it will usually be ['A', 'ADGSDFGSS', 'F']
            # So take index 1
            peptide = peptide_split[1]
        elif len(peptide_split)==1:
            # If we get 1 chunk it should just be ['ADGSDFGSS']
            # So take index 0
            peptide = peptide_split[0]
        else:
            # If we split the peptide and it is not length 1 or 3 then try to split with pro
            peptide = cls.split_peptide_pro(peptide_string=peptide_string, delimiter=delimiter)

        return(peptide)

    @classmethod
    def split_peptide_pro(cls, peptide_string, delimiter = "."):

        if delimiter!=".":
            front_regex = '^[A-Z|-][{}]'.format(delimiter)
            cls.FRONT_FLANKING_REGEX = re.compile(front_regex)
            back_regex = '[{}][A-Z|-]$'.format(delimiter)
            cls.BACK_FLANKING_REGEX = re.compile(back_regex)

        # Replace the front flanking with nothing
        peptide_string = cls.FRONT_FLANKING_REGEX.sub('', peptide_string)

        # Replace the back flanking with nothing
        peptide_string = cls.BACK_FLANKING_REGEX.sub('', peptide_string)

        return(peptide_string)





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
