#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 09:55:41 2017

@author: hinklet
"""

from Bio import SeqIO


class Fdr(object):
    
    def __init__(self):
        None
        
class SetBasedFdr(Fdr):
    """
    Class calculates set based FDR on the lead protein in the group
    Input is a DataStore object as well as an integer false discovery rate

    Example: py_protein_inference.fdrcalc.SetBasedFDR(data_class = data,false_discovery_rate=XX))

    FDR is calculated As (2*decoys)/total
    """
    def __init__(self,data_class,false_discovery_rate=.01):
        self.grouped_scored_data = data_class.grouped_scored_proteins
        self.false_discovery_rate = false_discovery_rate
        self.data_class = data_class
        
    def execute(self):
        
        #pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x[0] for x in self.grouped_scored_data]
        #Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        #Reverse the list (best to worst) -> (worst to best)
        lead_proteins.reverse()
        
        
        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_decoy_target_list = [ 1 if '#' in elem else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            decoys = sum(binary_decoy_target_list)
            #Calculate FDR at every step starting with the entire list...
            #Delete first entry (worst score) every time we go through a cycle
            fdr = (2*decoys)/(float(total))
            fdr_list.append(fdr)
            print fdr
            if fdr<self.false_discovery_rate:
                break
            else:
                #Here we delete the worst score every time... thus making our list smaller and smaller
                del lead_proteins[0]
        
        lead_proteins.reverse()

        self.fdr_list = fdr_list

        fdr_restricted_set = [self.grouped_scored_data[x] for x in range(len(lead_proteins))]

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        print 'Protein Group leads that pass with more than 1 PSM with a ' + str(self.false_discovery_rate) + ' FDR = ' + str(len(fdr_restricted_set)-len(onehitwonders))
        print 'Protein Group lead One hit Wonders that pass '+str(self.false_discovery_rate)+' FDR = '+str(len(onehitwonders))

        print 'Number of Protein groups that pass a '+str(self.false_discovery_rate*100)+' percent FDR: '+str(len(fdr_restricted_set))
        self.data_class.fdr_restricted_grouped_scored_proteins = fdr_restricted_set


class QValueCalculation(Fdr):
    """
    Class calculates Q values on the lead protein in the groups
    Input is a DataStore object

    Example: py_protein_inference.fdrcalc.SetBasedFDR(data_class = data)

    Q values are calculated As (2*decoys)/total
    """
    def __init__(self,data_class):
        self.group_objects = data_class.protein_group_objects
        self.data_class = data_class

    def execute(self):
        #pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x.proteins[0] for x in self.group_objects]
        #Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        lead_proteins.reverse()

        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_decoy_target_list = [1 if '#' in elem else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            decoys = sum(binary_decoy_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            fdr = (decoys * 2) / (float(total))
            fdr_list.append(fdr)
            del lead_proteins[0]

        qvalue_list = []
        new_fdr_list = []
        for fdrs in fdr_list:
            new_fdr_list.append(fdrs)
            qvalue = min(new_fdr_list)
            qvalue_list.append(qvalue)

        qvalue_list.reverse()

        for k in range(len(self.group_objects)):
            self.group_objects[k].q_value = qvalue_list[k]


class EntrapFdr(Fdr):
    """
    Class calculates Entrapment FDR on the lead protein in the groups.
    Input is a DataStore object, an entrapment database, and a false discovery rate

    Example: py_protein_inference.fdrcalc.EntrapFdr(data_class = data, entrapment_database = "example_entrap.fasta", other_database = None, false_discovery_rate=.05)

    FDR values are calculated As (entrapped proteins)/total

    This class is useful for calculating an entrapment FDR IF using a benchmark dataset with KNOWN protein content.
    Entrapment DB would be target proteins known to NOT be in the sample. However, these entrapped proteins should be in the main database that
    the search was searched against via comet/mascot
    """
    def __init__(self, data_class, true_database, other_database = None, false_discovery_rate=.05):
        self.grouped_scored_data = data_class.grouped_scored_proteins
        self.false_discovery_rate = false_discovery_rate
        self.data_class = data_class
        self.other_database = other_database
        self.true_database = true_database

    def execute(self):

        # entrapment_handle = SeqIO.parse(self.entrapment_database, 'fasta')
        #
        # entrapment_proteins = []
        # for records in entrapment_handle:
        #     if '#' not in records.id:
        #         entrapment_proteins.append(records.id)

        entrapment_proteins = []
        if self.other_database:
            other_handle = SeqIO.parse(self.other_database, 'fasta')
            other_proteins = []
            for records in other_handle:
                if '#' not in records.id:
                    other_proteins.append(records.id)
            entrapment_proteins = entrapment_proteins+other_proteins

        true_handle = SeqIO.parse(self.true_database, 'fasta')
        true_proteins = []
        for records in true_handle:
            if '#' not in records.id:
                true_proteins.append(records.id)

        protein_data = [x[0].identifier for x in self.grouped_scored_data]
        false_true_positives = []
        decoys = []
        decoy_fdr = []
        entrapment_fdr = []
        for i in range(len(protein_data)):
            # Get a list of false true positives from the protein data... basically if its not a decoy and if its not in true_proteins
            if '#' not in protein_data[i] and protein_data[i] not in true_proteins:
                false_true_positives.append(protein_data[i])
            if '#' in protein_data[i]:
                decoys.append(protein_data[i])


        entrapment_proteins_set = list(false_true_positives)
        entrapment_proteins_set = set(entrapment_proteins_set)

        # pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x[0] for x in self.grouped_scored_data]
        # Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        # Reverse the list (best to worst) -> (worst to best)
        lead_proteins.reverse()

        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_entrap_target_list = [1 if elem in entrapment_proteins_set else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            entraped = sum(binary_entrap_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            fdr = (entraped) / (float(total))
            fdr_list.append(fdr)
            print fdr
            if fdr < self.false_discovery_rate:
                break
            else:
                # Here we delete the worst score every time... thus making our list smaller and smaller
                del lead_proteins[0]

        lead_proteins.reverse()

        fdr_restricted_set = [self.grouped_scored_data[x] for x in range(len(lead_proteins))]

        self.fdr_list = fdr_list

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        print 'Protein Group leads that pass with more than 1 PSM with a ' + str(
            self.false_discovery_rate) + ' Entrapment FDR = ' + str(len(fdr_restricted_set) - len(onehitwonders))
        print 'Protein Group lead One hit Wonders that pass ' + str(self.false_discovery_rate) + ' Entrapment FDR = ' + str(
            len(onehitwonders))

        print 'Number of Protein groups that pass a ' + str(self.false_discovery_rate * 100) + ' Entrapment FDR: ' + str(
            len(fdr_restricted_set))
        self.data_class.fdr_restricted_grouped_scored_proteins = fdr_restricted_set

        self.restricted_proteins = fdr_restricted_set

        self.entrapment_proteins = entrapment_proteins


class PureDecoyFdr(Fdr):
    """
    Class calculates set based Decoy FDR on the lead protein in the group
    Input is a DataStore object as well as an integer false discovery rate

    Example: py_protein_inference.fdrcalc.PureDecoyFdr(data_class = data,false_discovery_rate=XX))

    FDR is calculated As (decoys)/total

    This class is used with py_protein_inference.fdrcalc.EntrapFdr in py_protein_inference.entrapment.GeneratePlot()
    The prior mentioned class generates plots of Decoy FDR vs Entrapment FDR to test accuracy of PI methods
    """
    def __init__(self, data_class, false_discovery_rate=.01):
        self.grouped_scored_data = data_class.grouped_scored_proteins
        self.false_discovery_rate = false_discovery_rate
        self.data_class = data_class

    def execute(self):

        # pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x[0] for x in self.grouped_scored_data]
        # Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        # Reverse the list (best to worst) -> (worst to best)
        lead_proteins.reverse()

        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_decoy_target_list = [1 if '#' in elem else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            decoys = sum(binary_decoy_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            fdr = (decoys) / (float(total))
            fdr_list.append(fdr)
            print fdr
            if fdr < self.false_discovery_rate:
                break
            else:
                # Here we delete the worst score every time... thus making our list smaller and smaller
                del lead_proteins[0]

        lead_proteins.reverse()

        self.fdr_list = fdr_list

        fdr_restricted_set = [self.grouped_scored_data[x] for x in range(len(lead_proteins))]

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        print 'Protein Group leads that pass with more than 1 PSM with a ' + str(
            self.false_discovery_rate) + ' FDR = ' + str(len(fdr_restricted_set) - len(onehitwonders))
        print 'Protein Group lead One hit Wonders that pass ' + str(self.false_discovery_rate) + ' FDR = ' + str(
            len(onehitwonders))

        print 'Number of Protein groups that pass a ' + str(self.false_discovery_rate * 100) + ' percent FDR: ' + str(
            len(fdr_restricted_set))
        self.data_class.fdr_restricted_grouped_scored_proteins = fdr_restricted_set
