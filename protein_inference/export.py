#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 10:49:47 2017

@author: hinklet
"""

import os
import csv
from logging import getLogger


class Export(object):
    """
    Main Class for exporting sorted protein lists to CSV
    """

    EXPORT_TYPES = ['leads', 'all', 'comma_sep', 'q_value_comma_sep', 'q_value', 'q_value_all']
    
    def __init__(self,data_class):
        self.data_class = data_class

    def export_to_csv(self, directory, export_type = "q_value"):

        logger = getLogger('protein_inference.export.Export.export_to_csv')

        data = self.data_class
        tag = data.parameter_file_object.tag

        if 'leads' in export_type:
            filename = '{}_leads_{}_{}.csv'.format(tag, data.short_score_method, data.score)
            logger.info("Exporting Protein Inference Data to File: {}".format(filename))
            complete_filepath = os.path.join(directory, filename)
            self.csv_export_leads_restricted(filename_out=complete_filepath)

        if 'all' in export_type:
            filename = '{}_all_{}_{}.csv'.format(tag, data.short_score_method, data.score)
            logger.info("Exporting Protein Inference Data to File: {}".format(filename))
            complete_filepath = os.path.join(directory, filename)
            self.csv_export_all_restricted(complete_filepath)

        if 'comma_sep' in export_type:
            filename = '{}_comma_sep_{}_{}.csv'.format(tag, data.short_score_method, data.score)
            logger.info("Exporting Protein Inference Data to File: {}".format(filename))
            complete_filepath = os.path.join(directory, filename)
            self.csv_export_comma_sep_restricted(complete_filepath)

        if 'q_value_comma_sep' in export_type:
            filename = '{}_q_value_comma_sep_{}_{}.csv'.format(tag, data.short_score_method, data.score)
            logger.info("Exporting Protein Inference Data to File: {}".format(filename))
            complete_filepath = os.path.join(directory, filename)
            self.csv_export_q_value_comma_sep(complete_filepath)

        if 'q_value' in export_type:
            filename = '{}_q_value_leads_{}_{}.csv'.format(tag, data.short_score_method, data.score)
            logger.info("Exporting Protein Inference Data to File: {}".format(filename))
            complete_filepath = os.path.join(directory, filename)
            self.csv_export_q_value_leads(complete_filepath)

        if 'q_value_all' in export_type:
            filename = '{}_q_value_all_{}_{}.csv'.format(tag, data.short_score_method, data.score)
            logger.info("Exporting Protein Inference Data to File: {}".format(filename))
            complete_filepath = os.path.join(directory, filename)
            self.csv_export_q_value_all(complete_filepath)

    def csv_export_all_restricted(self, filename_out):
        """
        Class that outputs a subset of the passing proteins based on FDR.
        Output proteins will be proteins that pass fdrcalc.SetBasedFdr(data_class = data,false_discovery_rate=.XX)
        Only Proteins that pass FDR will be output and ALL proteins
        will be output not just leads

        Example: protein_inference.export.CsvOutAll(data_class = data, filename_out = "example.csv")

        if we assume data = datastore.Datastore(reader_class = reader)
        and if we assume reader is a class object from reader.Reader()
        """
        data_to_write = self.data_class.fdr_restricted_grouped_scored_proteins
        ungrouped_list = [['Protein', 'Score', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in data_to_write:
            for prots in groups:
                ungrouped_list.append([prots.identifier])
                ungrouped_list[-1].append(prots.score)
                ungrouped_list[-1].append(prots.num_peptides)
                if prots.reviewed == True:
                    ungrouped_list[-1].append('Reviewed')
                else:
                    ungrouped_list[-1].append('Unreviewed')
                ungrouped_list[-1].append(prots.group_identification)
                for peps in prots.peptides:
                    ungrouped_list[-1].append(peps)

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)

    def csv_export_leads_restricted(self, filename_out):
        """
        Class that outputs a subset of the passing proteins based on FDR.
        Output proteins will be proteins that pass fdrcalc.SetBasedFdr(data_class = data,false_discovery_rate=.XX)
        Only Proteins that pass FDR will be output and only Lead proteins will be output

        Example: protein_inference.export.CsvOutLeads(data_class = data, filename_out = "example.csv")

        if we assume data = datastore.Datastore(reader_class = reader)
        and if we assume reader is a class object from reader.Reader()
        """
        ungrouped_list = [['Protein', 'Score', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_class.fdr_restricted_grouped_scored_proteins:
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

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)


    def csv_export_comma_sep_restricted(self, filename_out):
        """
        Class that outputs a subset of the passing proteins based on FDR.
        Output proteins will be proteins that pass fdrcalc.SetBasedFdr(data_class = data,false_discovery_rate=.XX)
        Only Proteins that pass FDR will be output and only Lead proteins will be output.
        Proteins in the groups of lead proteins will also be output in the same row as the lead

        Example: protein_inference.export.CsvOutCommaSep(data_class = data, filename_out = "example.csv")

        if we assume data = datastore.Datastore(reader_class = reader)
        and if we assume reader is a class object from reader.Reader()
        """
        ungrouped_list = [['Protein','Score','Number_of_Peptides','Identifier_Type','GroupID','Other_Potential_Identifiers']]
        for groups in self.data_class.fdr_restricted_grouped_scored_proteins:
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
        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
            
    def csv_export_q_value_leads(self, filename_out):
        """
        Class that outputs all lead proteins with Q values.

        Example: protein_inference.export.CsvOutLeadsQValues(data_class = data, filename_out = "example.csv")
        if we assume data = datastore.Datastore(reader_class = reader)

        and if we assume reader is a class object from reader.Reader()
        """
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_class.protein_group_objects:
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
            peptides = lead_protein.peptides
            for peps in sorted(peptides):
                ungrouped_list[-1].append(peps)

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)

    def csv_export_q_value_comma_sep(self, filename_out):
        """
        Class that outputs all lead proteins with Q values.
        Proteins in the groups of lead proteins will also be output in the same row as the lead

        Example: protein_inference.export.CsvOutCommaSepQValues(data_class = data, filename_out = "example.csv")
        if we assume data = datastore.Datastore(reader_class = reader)

        and if we assume reader is a class object from reader.Reader()
        """
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Other_Potential_Identifiers']]
        for groups in self.data_class.protein_group_objects:
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

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)


    def csv_export_q_value_all(self, filename_out):
        """
        Class that outputs all proteins with Q values.
        Non Lead proteins are also output - entire group gets output
        Proteins in the groups of lead proteins will also be output in the same row as the lead

        Example: protein_inference.export.CsvOutCommaSepQValues(data_class = data, filename_out = "example.csv")

        if we assume data = datastore.Datastore(reader_class = reader)
        and if we assume reader is a class object from reader.Reader()
        """
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_class.protein_group_objects:
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
                for peps in sorted(proteins.peptides):
                    ungrouped_list[-1].append(peps)

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)


    def csv_export_q_value_all_proteologic(self, filename_out):
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_class.protein_group_objects:
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

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)
            
    def csv_export_q_value_all_long(self, filename_out):
        """
        Class that outputs all lead proteins with Q values.

        Example: protein_inference.export.CsvOutLeadsQValues(data_class = data, filename_out = "example.csv")
        if we assume data = datastore.Datastore(reader_class = reader)

        and if we assume reader is a class object from reader.Reader()
        """
        ungrouped_list = [['Protein', 'Score', 'Q_Value', 'Number_of_Peptides', 'Identifier_Type', 'GroupID', 'Peptides']]
        for groups in self.data_class.protein_group_objects:
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

        with open(filename_out, "w") as f:
            writer = csv.writer(f)
            writer.writerows(ungrouped_list)