#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:25:00 2017

@author: hinklet
"""

import math
import numpy
import sys

class Score(object):
    
    def __init__(self):
        None
        
class BestPeptidePerProtein(Score):
     ###Input is DataStore object
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError('scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):
    
        
        all_scores = []

        print 'Scoring Proteins...'
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
        self.data_class.scored_proteins = all_scores
    
class FishersMethod(Score):
     ###Input is output from execute in reader class
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):


        all_scores = []
        print 'Scoring Proteins...'
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            score = -2*sum([math.log(x) for x in val_list])

            protein.score = score

            all_scores.append(protein)
            
        #Here reverse the sorting to descending because a higher score is better
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)
            
        self.data_class.score_method = 'fishers_method'
        self.data_class.scored_proteins = all_scores
            
class MultiplicativeLog(Score):
     ###Input is output from execute in reader class
     ###This scoring treats a lower q value or lower pep as GOOD
     ###We multiply all q values or all pep values for 1 protein together
     ###Then take -log(x) of the multiplied values
     ###Higher scores are better than lower scores
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    def execute(self):            

        all_scores = []
        print 'Scoring Proteins...'
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            
            combine = reduce(lambda x, y: x*y, val_list)
            if combine==0:
                combine=sys.float_info.min
            score = -math.log(combine)
            protein.score = score

            all_scores.append(protein)


        #Higher score is better as a smaller q or pep in a -log will give a larger value
        all_scores = sorted(all_scores,key = lambda k: k.score, reverse=True)
        
        self.data_class.score_method = 'multiplicative_log'
        self.data_class.scored_proteins = all_scores
        

    
class DownweightedMultiplicativeLog(Score):
     ###Input is output from execute in reader class
     ###This scoring treats a lower q value or lower pep as GOOD
     ###We multiply all q values or all pep values for 1 protein together
     ###Then this number is divided by the set q value mean raised to the number of peptides for that protein
     ###Then take -log(x) of the multiplied values
     ###Higher scores are better than lower scores
     
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
        print 'Scoring Proteins...'
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
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
        self.data_class.scored_proteins = all_scores
    
    
class TopTwoCombined(Score):
     ###Input is output from execute in reader class
     ###This scoring treats a lower q value or lower pep as GOOD
     ###Higher scores are better than lower scores
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):
            

        all_scores = []
        print 'Scoring Proteins...'
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
        self.data_class.scored_proteins = all_scores
    
class DownweightedVersion2(Score):
     ###Input is output from execute in reader class
     ###This scoring treats a lower q value or lower pep as GOOD
     ###Higher scores are better than lower scores
     
    def __init__(self,data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class
    
    def execute(self):


        all_scores = []
        print 'Scoring Proteins...'
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
        self.data_class.scored_proteins = all_scores


class IterativeDownweightedLog(Score):
    ###Input is output from execute in reader class
    ###This scoring treats a lower q value or lower pep as GOOD
    ###We multiply all q values or all pep values for 1 protein together
    ###Then this number is divided by the set q value mean raised to the number of peptides for that protein
    ###Then take -log(x) of the multiplied values
    ###Higher scores are better than lower scores

    def __init__(self, data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class

    def execute(self):

        all_scores = []
        print 'Scoring Proteins...'
        for protein in self.pre_score_data:
            val_list = []
            for vals in protein.psm_score_dictionary:
                val_list.append(float(vals['score']))
            mean = numpy.mean(val_list)
            # Here take each score and raise it to the power of (1/(1+index_number)).
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
        self.data_class.scored_proteins = all_scores


class GeometricMeanLog(Score):
    ###Input is output from execute in reader class
    ###This scoring treats a lower q value or lower pep as GOOD
    ###We multiply all q values or all pep values for 1 protein together
    ###Then this number is then raised to the 1/x where x is the number of values (peptide scores) for the specific protein
    ###Then take -log(x) of the multiplied values
    ###Higher scores are better than lower scores

    def __init__(self, data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class

    def execute(self):

        all_scores = []
        print 'Scoring Proteins...'
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
        self.data_class.scored_proteins = all_scores


class IterativeDownweightingV2(Score):
    ###Input is output from execute in reader class
    ###This scoring treats a lower q value or lower pep as GOOD
    ###Higher scores are better than lower scores

    def __init__(self, data_class):
        if data_class.scoring_input:
            self.pre_score_data = data_class.scoring_input
        else:
            raise ValueError(
                'scoring input not found in data class - Please run PreScoreQValue of PreScorePepValue from DataStore to run any scoring type')
        self.data_class = data_class

    def execute(self):

        all_scores = []
        print 'Scoring Proteins...'
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

        self.data_class.score_method = 'iterative_downweighting'
        self.data_class.scored_proteins = all_scores