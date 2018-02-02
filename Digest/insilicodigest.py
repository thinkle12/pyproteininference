#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:25:00 2017

@author: hinklet
"""

import collections
from Bio import SeqIO
import os
import cPickle as pickle

class InSilicoDigest(object):
    
    def __init__(self,database_path,num_miss_cleavs=2,digest_type='trypsin'):
        self.peptide_to_protein_dictionary = None
        self.protein_to_peptide_dictionary = None
        self.protein_peptide_complete_set = None
        self.swiss_prot_protein_dictionary = None
        self.database_path = database_path
        self.num_miss_cleavs = num_miss_cleavs
        self.list_of_digest_types = ['trypsin']
        if digest_type in self.list_of_digest_types:
            self.digest_type = digest_type
        else:
            raise ValueError('digest_type must be equal to one of the following'+str(self.list_of_digest_types)+' or... (List more digest types here in the future...)')
        
        
    def digest(self,proseq,miss_cleavage):
        peptides=[]
        cut_sites=[0]
        if self.digest_type == 'trypsin':
            for i in range(0,len(proseq)-1):
                if proseq[i]=='K' and proseq[i+1]!='P':
                    cut_sites.append(i+1)
                elif proseq[i]=='R' and proseq[i+1]!='P':
                    cut_sites.append(i+1)

        if self.digest_type=='other_type':
            pass
            #Here write more code for other types of digest types....
        
        if cut_sites[-1]!=len(proseq):
            cut_sites.append(len(proseq))
    
        if len(cut_sites)>2:
            if  miss_cleavage==0:
                for j in range(0,len(cut_sites)-1):
                    peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
    
            elif miss_cleavage==1:
                for j in range(0,len(cut_sites)-2):
                    peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                    peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
                
                peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    
            elif miss_cleavage==2:
                for j in range(0,len(cut_sites)-3):
                    peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                    peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
                    peptides.append(proseq[cut_sites[j]:cut_sites[j+3]])
                
                peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
                peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
                peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
        else: #there is no tryptic site in the protein sequence
            peptides.append(proseq)
        return peptides



    def execute(self):
        db_path_only = '/'.join(self.database_path.split('/')[:-1])+'/'
        pickle_filename_tag = self.database_path.split('.')[0] + '_' + str(self.digest_type) + '_' + str(self.num_miss_cleavs)
        if pickle_filename_tag.split('/')[-1]+'_pep_to_prot.pickle' not in os.listdir(db_path_only) or pickle_filename_tag.split('/')[-1]+'_prot_to_pep.pickle' not in os.listdir(db_path_only) or pickle_filename_tag.split('/')[-1]+'_sp_dict.pickle' not in os.listdir(db_path_only):
            handle=SeqIO.parse(self.database_path,'fasta')

            print 'Starting Digest...'
            pep_dict = collections.defaultdict(set)
            prot_dict = collections.defaultdict(set)
            sp_dict = collections.defaultdict(list)





            # We use [:2] because that is the first two letters of the protein identifier
            # We use [3:] below also because we need to remove "sp|" or "tr|" to match what the search results show....
            for record in handle:
                identifier_stripped = record.id[3:]

                if record.id[:2] == 'sp':
                    sp_dict['swiss-prot'].append(identifier_stripped)
                proseq = str(record.seq)
                peptide_list = InSilicoDigest(self.database_path,self.num_miss_cleavs,self.digest_type).digest(proseq, self.num_miss_cleavs)
                for peptide in peptide_list:
                    pep_dict[peptide].add(identifier_stripped)
                    prot_dict[identifier_stripped].add(peptide)

            handle.close()


            # print 'Starting to pickle dictionaries...'
            #
            # pickle_out = open(pickle_filename_tag + '_pep_to_prot.pickle', "wb")
            # pickle.dump(pep_dict, pickle_out)
            # pickle_out.close()
            #
            # print 'pep to prot dictionary has been pickled...'
            #
            # pickle_out = open(pickle_filename_tag + '_prot_to_pep.pickle', "wb")
            # pickle.dump(prot_dict, pickle_out)
            # pickle_out.close()
            #
            # print 'prot to pep dictionary has been pickled...'
            #
            # pickle_out = open(pickle_filename_tag + '_sp_dict.pickle', "wb")
            # pickle.dump(sp_dict, pickle_out)
            # pickle_out.close()
            #
            # print 'sp dictionary has been pickled...'

            print 'Digest finished, peptide and protein dictionaries created based on the provided database'


        else:
            print 'Skipping Digest, Importing digest information from'+'\n'+pickle_filename_tag

            pickle_in = open(pickle_filename_tag + '_pep_to_prot.pickle', "rb")
            pep_dict = pickle.load(pickle_in)
            pickle_in.close()

            pickle_in = open(pickle_filename_tag + '_prot_to_pep.pickle', "rb")
            prot_dict = pickle.load(pickle_in)
            pickle_in.close()

            pickle_in = open(pickle_filename_tag + '_sp_dict.pickle', "rb")
            sp_dict = pickle.load(pickle_in)
            pickle_in.close()

            print 'Pickle Files Loaded...'


        self.swiss_prot_protein_dictionary = sp_dict
        self.peptide_to_protein_dictionary = pep_dict
        self.protein_to_peptide_dictionary = prot_dict




                
