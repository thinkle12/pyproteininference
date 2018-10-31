#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:14:28 2017

@author: hinklet
"""

import datastore

class Picker(object):
    """
    Main Parent Picker class
    """
    
    def __init__(self):
        None
        
class StandardPicker(Picker):
    """
    The following class executes Protein Picker.

    Protein picker compares scores of target and decoy protein pairs, IE:
    PRKDC_HUMAN|P78527 vs ##PRKDC_HUMAN|##P78527.

    Scores of the Target and Decoy are compared and the Protein with the worse score of the two is completely removed from further analysis.

    This sort of analysis is built into other protein inference tools such as Percolator Built in Inference.
    It seems to be a good way of dealing with the TITIN problem as for a lot of searches we will see TITIN as well as Decoy ##TITIN

    Example: py_protein_inference.picker.StandardPicker(data_class = data)

    Where data is a DataStore Object.
    """
    
    def __init__(self,data_class):
        self.scored_data = data_class.scored_proteins
        self.protein_index = 0
        self.protein_score_index = 1
        self.data_class = data_class
        
    def execute(self):

        #Use higher or lower class to determine if a higher protein score or lower protein score is better based on the scoring method used
        high_low = datastore.HigherOrLower(self.data_class)
        high_low.execute()
        higher_or_lower = self.data_class.high_low_better
        #Here we determine if a lower or higher score is better
        #Since all input is ordered from best to worst we can do the following

        
        index_to_remove = []
        #data_class.scored_proteins is simply a list of Protein objects...
        #Create list of all decoy proteins
        decoy_proteins = [x.identifier for x in self.scored_data if '##' in x.identifier]
        #Create a list of all potential matching targets (some of these may not exist in the search)
        matching_targets = [''.join(x.split('##')) for x in decoy_proteins]
        
        #Create a list of all the proteins from the scored data
        all_proteins = [x.identifier for x in self.scored_data]
        print str(len(all_proteins))+' proteins scored'
        
        total_targets = []
        total_decoys = []
        decoys_removed = []
        targets_removed = []
        #Loop over all decoys identified in the search
        print 'Picking Proteins...'
        for i in range(len(decoy_proteins)):
            cur_decoy_index = all_proteins.index(decoy_proteins[i])
            cur_decoy_protein_object = self.scored_data[cur_decoy_index]
            total_decoys.append(cur_decoy_protein_object.identifier)
            
            #Try, Except here because the matching target to the decoy may not be a result from the search
            try:
                cur_target_index = all_proteins.index(matching_targets[i])
                cur_target_protein_object = self.scored_data[cur_target_index]
                total_targets.append(cur_target_protein_object.identifier)
                
                if higher_or_lower=='higher':
                    if cur_target_protein_object.score>cur_decoy_protein_object.score:
                        index_to_remove.append(cur_decoy_index)
                        decoys_removed.append(cur_decoy_index)
                        cur_target_protein_object.picked=True
                        cur_decoy_protein_object.picked=False
                    else:
                        index_to_remove.append(cur_target_index)
                        targets_removed.append(cur_target_index)
                        cur_decoy_protein_object.picked=True
                        cur_target_protein_object.picked=False
                        
                if higher_or_lower=='lower':
                    if cur_target_protein_object.score<cur_decoy_protein_object.score:
                        index_to_remove.append(cur_decoy_index)
                        decoys_removed.append(cur_decoy_index)
                        cur_target_protein_object.picked=True
                        cur_decoy_protein_object.picked=False
                    else:
                        index_to_remove.append(cur_target_index)
                        targets_removed.append(cur_target_index)
                        cur_decoy_protein_object.picked=True
                        cur_target_protein_object.picked=False
            except ValueError:
                pass
        
        print str(len(total_decoys))+' total decoy proteins'
        print str(len(total_targets))+' matching target proteins also found in search'
        print str(len(decoys_removed))+' decoy proteins to be removed'
        print str(len(targets_removed))+' target proteins to be removed'
        
        
        print 'Removing Lower Scoring Proteins...'
        picked_list = []
        removed_proteins = []
        for protein_objects in self.scored_data:
            if protein_objects.picked==True:
                picked_list.append(protein_objects)
            else:
                removed_proteins.append(protein_objects)
        self.data_class.picked_proteins_scored = picked_list
        self.data_class.picked_proteins_removed = removed_proteins
        print 'Finished Removing Proteins'


