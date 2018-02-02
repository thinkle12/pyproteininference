#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 10:49:47 2017

@author: hinklet
"""

import csv

class Export(object):
    
    def __init__(self):
        None
        
class CsvOutAll(Export):
    
    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.fdr_restricted_grouped_scored_proteins
        self.filename_out = filename_out
        
    def execute(self):
        ungrouped_list = [['Protein','Score','Number_of_Peptides','Identifier_Type','GroupID','Peptides']]
        for groups in self.data_to_write:
            for prots in groups:
                ungrouped_list.append([prots.identifier])
                ungrouped_list[-1].append(prots.score)
                ungrouped_list[-1].append(prots.num_peptides)
                if prots.reviewed==True:
                    ungrouped_list[-1].append('Reviewed')
                else:
                    ungrouped_list[-1].append('Unreviewed')
                ungrouped_list[-1].append(prots.group_identification)
                ungrouped_list[-1].append(prots.peptides)
        
        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
            
class CsvOutLeads(Export):
    
    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.fdr_restricted_grouped_scored_proteins
        self.filename_out = filename_out
        
    def execute(self):
        ungrouped_list = [['Protein','Score','Number_of_Peptides','Identifier_Type','GroupID','Peptides']]
        for groups in self.data_to_write:
            ungrouped_list.append([groups[0].identifier])
            ungrouped_list[-1].append(groups[0].score)
            ungrouped_list[-1].append(groups[0].num_peptides)
            if groups[0].reviewed == True:
                ungrouped_list[-1].append('Reviewed')
            else:
                ungrouped_list[-1].append('Unreviewed')
            ungrouped_list[-1].append(groups[0].group_identification)
            ungrouped_list[-1].append(groups[0].peptides)
        
        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
            
class CsvOutCommaSep(Export):
    
    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.fdr_restricted_grouped_scored_proteins
        self.filename_out = filename_out
        
    def execute(self):
        ungrouped_list = [['Protein','Score','Number_of_Peptides','GroupID','Other_Potential_Identifiers','Identifier_Type']]
        for groups in self.data_to_write:
            for prots in groups:
                if prots==groups[0]:
                    ungrouped_list.append([prots[0]])
                    ungrouped_list[-1].append(prots[1])
                    ungrouped_list[-1].append(prots[2])                
                    ungrouped_list[-1].append(prots[3])
                    ungrouped_list[-1].append(prots[4])
                else:
                    ungrouped_list[-1].append(prots[0])
        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)


class CsvOutLeadsQValues(Export):

    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.protein_group_objects
        self.filename_out = filename_out

    def execute(self):
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_to_write:
            lead_protein = groups.proteins[0]
            ungrouped_list.append([lead_protein.identifier])
            ungrouped_list[-1].append(lead_protein.score)
            ungrouped_list[-1].append(groups.q_value)
            ungrouped_list[-1].append(lead_protein.num_peptides)
            if lead_protein.reviewed == True:
                ungrouped_list[-1].append('Reviewed')
            else:
                ungrouped_list[-1].append('Unreviewed')
            ungrouped_list[-1].append(groups.number_id)
            ungrouped_list[-1].append(lead_protein.peptides)

        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
