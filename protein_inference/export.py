#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 10:49:47 2017

@author: hinklet
"""

import csv

class Export(object):
    # TODO make this an abstract base class...
    """
    Main Class for exporting sorted protein lists to CSV
    """
    
    def __init__(self):
        None
        
class CsvOutAll(Export):
    """
    Class that outputs a subset of the passing proteins based on FDR.
    Output proteins will be proteins that pass fdrcalc.SetBasedFdr(data_class = data,false_discovery_rate=.XX)
    Only Proteins that pass FDR will be output and ALL proteins
    will be output not just leads

    Example: protein_inference.export.CsvOutAll(data_class = data, filename_out = "example.csv")

    if we assume data = datastore.Datastore(reader_class = reader)
    and if we assume reader is a class object from reader.Reader()
    """
    
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
                for peps in prots.peptides:
                    ungrouped_list[-1].append(peps)
        
        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
            
class CsvOutLeads(Export):
    """
    Class that outputs a subset of the passing proteins based on FDR.
    Output proteins will be proteins that pass fdrcalc.SetBasedFdr(data_class = data,false_discovery_rate=.XX)
    Only Proteins that pass FDR will be output and only Lead proteins will be output

    Example: protein_inference.export.CsvOutLeads(data_class = data, filename_out = "example.csv")

    if we assume data = datastore.Datastore(reader_class = reader)
    and if we assume reader is a class object from reader.Reader()
    """
    
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
            for peps in groups[0].peptides:
                ungrouped_list[-1].append(peps)
        
        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
            
class CsvOutCommaSep(Export):
    """
    Class that outputs a subset of the passing proteins based on FDR.
    Output proteins will be proteins that pass fdrcalc.SetBasedFdr(data_class = data,false_discovery_rate=.XX)
    Only Proteins that pass FDR will be output and only Lead proteins will be output.
    Proteins in the groups of lead proteins will also be output in the same row as the lead

    Example: protein_inference.export.CsvOutCommaSep(data_class = data, filename_out = "example.csv")

    if we assume data = datastore.Datastore(reader_class = reader)
    and if we assume reader is a class object from reader.Reader()
    """
    
    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.fdr_restricted_grouped_scored_proteins
        self.filename_out = filename_out
        
    def execute(self):
        ungrouped_list = [['Protein','Score','Number_of_Peptides','Identifier_Type','GroupID','Other_Potential_Identifiers']]
        for groups in self.data_to_write:
            for prots in groups:
                if prots==groups[0]:
                    ungrouped_list.append([prots.identifier])
                    ungrouped_list[-1].append(prots.score)
                    ungrouped_list[-1].append(prots.num_peptides)
                    if prots.reviewed == True:
                        ungrouped_list[-1].append('Reviewed')
                    else:
                        ungrouped_list[-1].append('Unreviewed')
                    ungrouped_list[-1].append(prots.group_identification)
                else:
                    ungrouped_list[-1].append(prots.identifier)
        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)


class CsvOutLeadsQValues(Export):
    """
    Class that outputs all lead proteins with Q values.

    Example: protein_inference.export.CsvOutLeadsQValues(data_class = data, filename_out = "example.csv")
    if we assume data = datastore.Datastore(reader_class = reader)

    and if we assume reader is a class object from reader.Reader()
    """

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
            for peps in lead_protein.peptides:
                ungrouped_list[-1].append(peps)

        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)

class CsvOutCommaSepQValues(Export):
    """
    Class that outputs all lead proteins with Q values.
    Proteins in the groups of lead proteins will also be output in the same row as the lead

    Example: protein_inference.export.CsvOutCommaSepQValues(data_class = data, filename_out = "example.csv")
    if we assume data = datastore.Datastore(reader_class = reader)

    and if we assume reader is a class object from reader.Reader()
    """

    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.protein_group_objects
        self.filename_out = filename_out

    def execute(self):
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Other_Potential_Identifiers']]
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
            for other_prots in groups.proteins[1:]:
                ungrouped_list[-1].append(other_prots.identifier)

        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)

class CsvOutAllQValues(Export):
    """
    Class that outputs all proteins with Q values.
    Non Lead proteins are also output - entire group gets output
    Proteins in the groups of lead proteins will also be output in the same row as the lead

    Example: protein_inference.export.CsvOutCommaSepQValues(data_class = data, filename_out = "example.csv")

    if we assume data = datastore.Datastore(reader_class = reader)
    and if we assume reader is a class object from reader.Reader()
    """

    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.protein_group_objects
        self.filename_out = filename_out

    def execute(self):
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_to_write:
            for proteins in groups.proteins:
                ungrouped_list.append([proteins.identifier])
                ungrouped_list[-1].append(proteins.score)
                ungrouped_list[-1].append(groups.q_value)
                ungrouped_list[-1].append(proteins.num_peptides)
                if proteins.reviewed == True:
                    ungrouped_list[-1].append('Reviewed')
                else:
                    ungrouped_list[-1].append('Unreviewed')
                ungrouped_list[-1].append(groups.number_id)
                for peps in proteins.peptides:
                    ungrouped_list[-1].append(peps)

        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)

class ProteologicAllQValues(Export):

    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.protein_group_objects
        self.filename_out = filename_out

    def execute(self):
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_to_write:
            for proteins in groups.proteins:
                ungrouped_list.append([proteins.identifier])
                ungrouped_list[-1].append(proteins.score)
                ungrouped_list[-1].append(groups.q_value)
                ungrouped_list[-1].append(proteins.num_peptides)
                if proteins.reviewed == True:
                    ungrouped_list[-1].append('Reviewed')
                else:
                    ungrouped_list[-1].append('Unreviewed')
                ungrouped_list[-1].append(groups.number_id)
                for peps in proteins.peptides:
                    ungrouped_list[-1].append(peps)

        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)


class CsvOutLeadsQValuesLong(Export):
    """
    Class that outputs all lead proteins with Q values.

    Example: protein_inference.export.CsvOutLeadsQValues(data_class = data, filename_out = "example.csv")
    if we assume data = datastore.Datastore(reader_class = reader)

    and if we assume reader is a class object from reader.Reader()
    """

    def __init__(self,data_class,filename_out):
        self.data_to_write = data_class.protein_group_objects
        self.filename_out = filename_out

    def execute(self):
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_to_write:
            lead_protein = groups.proteins[0]
            for peps in lead_protein.peptides:
                ungrouped_list.append([lead_protein.identifier])
                ungrouped_list[-1].append(lead_protein.score)
                ungrouped_list[-1].append(groups.q_value)
                ungrouped_list[-1].append(lead_protein.num_peptides)
                if lead_protein.reviewed == True:
                    ungrouped_list[-1].append('Reviewed')
                else:
                    ungrouped_list[-1].append('Unreviewed')
                ungrouped_list[-1].append(groups.number_id)
                ungrouped_list[-1].append(peps)

        with open(self.filename_out, "wb") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)