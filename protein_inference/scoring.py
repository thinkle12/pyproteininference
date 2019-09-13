#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:25:00 2017

@author: hinklet
"""

import math
import numpy
import sys
from functools import reduce

class Score(object):
    # TODO make this an abstract base class...

    """
    Parent Score class for all scoring subset classes
    """
    
    def __init__(self):
        None
        
class BestPeptidePerProtein(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a best peptide per protein scoring scheme

    Example: protein_inference.scoring.BestPeptidePerProtein(data_class = data)

    Where data is a DataStore Object
     """
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError('scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):
    
        
        all_scores = []

        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = []
            temp_peps = []
            for vals in protein.psm_score_dictionary:
                temp_peps.append(vals['peptide'])
                val_list.append(float(vals['score']))
            score = min([float(x) for x in val_list])

            protein.score = score

            
            all_scores.append(protein)
        #Here do ascending sorting because a lower pep or q value is better
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=False)


        self.data_class.score_method = 'best_peptide_per_protein'
        self.data_class.short_score_method = 'bppp'
        self.data_class.scored_proteins = all_scores
    
class FishersMethod(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Fishers method scoring scheme

    Example: protein_inference.scoring.FishersMethod(data_class = data)

    Where data is a DataStore Object
     """
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):


        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = (x['score'] for x in protein.psm_score_dictionary)
            score = -2*sum([math.log(x) for x in val_list])

            protein.score = score

            all_scores.append(protein)
            
        #Here reverse the sorting to descending because a higher score is better
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)
            
        self.data_class.score_method = 'fishers_method'
        self.data_class.short_score_method = 'fm'
        self.data_class.scored_proteins = all_scores
            
class MultiplicativeLog(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Multiplicative Log scoring scheme.
    All the Qvalues/PepValues from all the peptides per protein are multiplied together and we take -Log(X) of the multiplied Peptide scores

    Example: protein_inference.scoring.MultiplicativeLog(data_class = data)

    Where data is a DataStore Object
     """
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    def execute(self):            
        # Instead of making all_scores a list... make it a generator??

        all_scores = []
        print('Scoring Proteins...')
        print('Using Generators')
        for protein in self.pre_score_data:
            # We create a generator of val_list...
            val_list = (x['score'] for x in protein.psm_score_dictionary)
            
            combine = reduce(lambda x, y: x*y, val_list)
            if combine==0:
                combine=sys.float_info.min
            score = -math.log(combine)
            protein.score = score

            all_scores.append(protein)


        #Higher score is better as a smaller q or pep in a -log will give a larger value
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)
        
        self.data_class.score_method = 'multiplicative_log'
        self.data_class.short_score_method = 'ml'
        self.data_class.scored_proteins = all_scores
        

    
class DownweightedMultiplicativeLog(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Multiplicative Log scoring scheme.
    All the Qvalues/PepValues from all the peptides per protein are multiplied together and
    then this number is divided by the set QValue/PepValue mean raised to the number of peptides for that protein
    then we take -Log(X) of the following value

    Example: protein_inference.scoring.DownweightedMultiplicativeLog(data_class = data)

    Where data is a DataStore Object
     """
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):


        score_list = []
        for proteins in self.pre_score_data:
            cur_score_dict = proteins.psm_score_dictionary
            for scores in cur_score_dict:
                score_list.append(float(scores['score']))
        score_mean = numpy.mean(score_list)


        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = [x['score'] for x in protein.psm_score_dictionary]
            #Divide by the score mean raised to the length of the number of unique peptides for the protein
            #This is an attempt to normalize for number of peptides per protein
            combine = reduce(lambda x, y: x*y, val_list)
            if combine==0:
                combine=sys.float_info.min
            score = -math.log(combine/(score_mean**len(val_list)))
            protein.score = score

            all_scores.append(protein)

        #Higher score is better as a smaller q or pep in a -log will give a larger value
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)

        self.data_class.score_method = 'downweighted_multiplicative_log'
        self.data_class.short_score_method = 'dwml'
        self.data_class.scored_proteins = all_scores
    
    
class TopTwoCombined(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Top Two scoring scheme.
    The top two scores for each protein are multiplied together and we take -Log(X) of the  multiplied value.
    If a protein only has 1 score/peptide, then we only do -Log(X) of the 1 peptide score

    Example: protein_inference.scoring.TopTwoCombined(data_class = data)

    Where data is a DataStore Object
     """
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):
            

        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            try:
                #Try to combine the top two scores
                #Divide by 2 to attempt to normalize the value
                score = -math.log((val_list[0]*val_list[1])/2)
            except IndexError:
                #If there is only 1 score/1 peptide then just use the 1 peptide provided
                score = -math.log(val_list[0])


            protein.score = score
            all_scores.append(protein)
            
        #Higher score is better as a smaller q or pep in a -log will give a larger value    
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)
        
        self.data_class.score_method = 'top_two_combined'
        self.data_class.short_score_method = 'ttc'
        self.data_class.scored_proteins = all_scores
    
class DownweightedVersion2(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Downweighted Multiplicative Log scoring scheme.
    Each peptide is iteratively downweighted by raising the peptide QValue or PepValue to the following power (1/(1+index_number)).
    Where index_number is the peptide number per protein...
    Each score for a protein provides less and less weight iteratively

    We also take -Log(X) of the final score here

    Example: protein_inference.scoring.DownweightedVersion2(data_class = data)

    Where data is a DataStore Object
     """
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):


        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            #Here take each score and raise it to the power of (1/(1+index_number)).
            #This downweights each successive score by reducing its weight in a decreasing fashion
            #Basically, each score for a protein will provide less and less weight iteratively
            val_list = [val_list[x]**(1/float(1+x)) for x in range(len(val_list))]
            #val_list = [val_list[x]**(1/float(1+(float(x)/10))) for x in range(len(val_list))]
            score = -math.log(reduce(lambda x, y: x*y, val_list))
                
                
            protein.score = score
            all_scores.append(protein)

        #Higher score is better as a smaller q or pep in a -log will give a larger value                
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)
        
        self.data_class.score_method = 'downweighted_version2'
        self.data_class.short_score_method = 'dw2'
        self.data_class.scored_proteins = all_scores


class IterativeDownweightedLog(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Downweighted Multiplicative Log scoring scheme.
    Each peptide is iteratively downweighted by multiplying the peptide QValue or PepValue to the following  (1+index_number).
    Where index_number is the peptide number per protein...
    Each score for a protein provides less and less weight iteratively

    We also take -Log(X) of the final score here

    Example: protein_inference.scoring.IterativeDownweightedLog(data_class = data)

    Where data is a DataStore Object
     """

    def __init__(self, data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class

    def execute(self):

        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            mean = numpy.mean(val_list)
            # Here take each score and multiply it by its index number).
            # This downweights each successive score by reducing its weight in a decreasing fashion
            # Basically, each score for a protein will provide less and less weight iteratively
            val_list = [val_list[x] * (float(1 + x)) for x in range(len(val_list))]
            # val_list = [val_list[x]**(1/float(1+(float(x)/10))) for x in range(len(val_list))]
            combine = reduce(lambda x, y: x*y, val_list)
            if combine==0:
                combine=sys.float_info.min
            score = -math.log(combine)
            protein.score = score

            protein.score = score
            all_scores.append(protein)

        # Higher score is better as a smaller q or pep in a -log will give a larger value
        all_scores = sorted(all_scores, key=lambda k: k.score, reverse=True)

        self.data_class.score_method = 'iterative_downweighted_log'
        self.data_class.short_score_method = 'idwl'
        self.data_class.scored_proteins = all_scores


class GeometricMeanLog(Score):
    """
    Class scores all proteins from data_class.scoring_input.
    This is either Pep or Q values from an unrestricted or restricted subset of the input scores from Percolator

    This class uses a Geometric Mean scoring scheme.

    We also take -Log(X) of the final score here

    Example: protein_inference.scoring.GeometricMeanLog(data_class = data)

    Where data is a DataStore Object
     """


    def __init__(self, data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class

    def execute(self):

        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
                combine = reduce(lambda x, y: x * y, val_list)
                if combine == 0:
                    combine = sys.float_info.min
                pre_log_score = combine**(1/float(len(val_list)))
            score = -math.log(pre_log_score)

            protein.score = score
            all_scores.append(protein)

        # Higher score is better as a smaller q or pep in a -log will give a larger value
        all_scores = sorted(all_scores, key=lambda k: k.score, reverse=True)

        self.data_class.score_method = 'geometric_mean_log'
        self.data_class.short_score_method = 'gm'
        self.data_class.scored_proteins = all_scores


class IterativeDownweightingV2(Score):
    """
    The following class is an experimental class essentially used for future development of potential scoring schemes
    """
    def __init__(self, data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class

    def execute(self):

        all_scores = []
        print('Scoring Proteins...')
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            # Here take each score and raise it to the power of (1/(1+index_number)).
            # This downweights each successive score by reducing its weight in a decreasing fashion
            # Basically, each score for a protein will provide less and less weight iteratively
            val_list = [val_list[x] ** (1 / float(1 + x)) for x in range(len(val_list))]
            # val_list = [val_list[x]**(1/float(1+(float(x)/10))) for x in range(len(val_list))]
            score = -math.log(reduce(lambda x, y: x * y, val_list))

            protein.score = score
            all_scores.append(protein)

        # Higher score is better as a smaller q or pep in a -log will give a larger value
        all_scores = sorted(all_scores, key=lambda k: k.score, reverse=True)

        self.data_class.score_method = 'iterative_downweighting2'
        self.data_class.short_score_method = 'idw2'
        self.data_class.scored_proteins = all_scores