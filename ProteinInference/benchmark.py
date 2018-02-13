#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 10:15:19 2017

@author: hinklet
"""

import ProteinInference
import matplotlib
matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import numpy
import math
from matplotlib.backends.backend_pdf import PdfPages

class Benchmark(object):
    """
    Generic Benchmark parent class
    """

    def __init__(self):
        None

class RocPlot(Benchmark):
    """
    Class for producing a Roc plot.
    Requires a DataStore object

    Example: rp = ProteinInference.benchmark.RocPlot(data_class = data)
             rp.execute(pdf=True)


    Where data is a DataStore Object
    """

    def __init__(self,data_class):
        fdr01 = ProteinInference.fdrcalc.SetBasedFdr(data_class=data_class, false_discovery_rate=.01)
        fdr01.execute()
        one_percent_data = data_class.fdr_restricted_grouped_scored_proteins
        self.one_percent_length = len(one_percent_data)
        fdr1 = ProteinInference.fdrcalc.SetBasedFdr(data_class=data_class, false_discovery_rate=1)
        fdr1.execute()
        self.all_data = data_class.fdr_restricted_grouped_scored_proteins
        self.data_class = data_class

    def execute(self,pdf=False):


        protein_data = [x[0].identifier for x in self.all_data]


        spec_list = []
        sens_list = []
        numlist = []
        posprots = []
        reg_spec = []
        f1_scores = []
        prec_list = []

        for i in range(len(protein_data)):
            negprots = []
            posprots.append(protein_data[i])
            numlist.append(i)
            positive_ind = set(range(len(protein_data)))-set(numlist)
            for j in positive_ind:
                negprots.append(protein_data[j])
        #    fp = []
        #    tn = []
            numlistneg = [ 1 if '#' in elem else 0 for elem in negprots]
            numlistpos = [ 1 if '#' in elem else 0 for elem in posprots]
            fp = sum(numlistpos)
            tp = len(numlistpos)-fp
            tn = sum(numlistneg)
            fn = len(numlistneg)-tn

            sensitivity = float(tp)/(tp+fn)
            specificity = 1-(float(tn)/(tn+fp))
            precision = float(fp)/(tp+fp)

            try:
                f1 = 2/((1/sensitivity)+(1/precision))
            except ZeroDivisionError:
                f1 = 0

            sens_list.append(sensitivity)
            spec_list.append(specificity)
            prec_list.append(precision)
            reg_spec.append(float(tn)/(tn+fp))
            f1_scores.append(f1)


        #    for hits in negprots:
        #        if '#' in hits[0]:
        #            tn.append(hits)
        #        else:
        #            fp.append(hits)
        #    try:
        #        fpr = float(fp)/(tn+fp)
        #    except ZeroDivisionError:
        #        fpr = float(0)

        #    tp = []
        #    fn = []

        #    for hits in posprots:
        #        if '#' in hits[0]:
        #            fn.append(hits)
        #        else:
        #            tp.append(hits)
        #    tpr = float(tp)/(tp+fn)
        #    fpr_list.append(fpr)
        #    tpr_list.append(tpr)


        maxf1 = max(f1_scores)

        print 'F1 max = '+str(maxf1)

        area = numpy.trapz(sens_list,spec_list)

        index_of_union = [abs(sens_list[x]-area)+abs(reg_spec[x]-area) for x in range(len(reg_spec))]

        min_iu = min(index_of_union)

        index_of_union_index = index_of_union.index(min_iu)

        print 'Index of Union 1-specificty = '+str(spec_list[index_of_union_index])
        print 'Index of Union sensitivity = '+str(sens_list[index_of_union_index])

        concord_prob = [reg_spec[x]*sens_list[x] for x in range(len(reg_spec))]

        concord_max = max(concord_prob)

        concord_ind = concord_prob.index(concord_max)

        print 'Concordance Probability 1-specificty = '+str(spec_list[concord_ind])
        print 'Concordance Probability sensitivity = '+str(sens_list[concord_ind])

        youdens = [reg_spec[x]+sens_list[x]-1 for x in range(len(reg_spec))]

        youdens_ind = max(youdens)

        youind = youdens.index(youdens_ind)

        print 'Youdens 1-specificty = '+str(spec_list[youind])
        print 'Youdens sensitivity = '+str(sens_list[youind])

        coordinate = [[sens_list[x],spec_list[x]] for x in range(len(sens_list))]

        dist_list = [math.hypot(0 - x[1], 1 - x[0]) for x in coordinate]

        min_dist = min(dist_list)

        min_ind = dist_list.index(min_dist)

        best_coord = coordinate[min_ind]

        print 'dist = '+str(min_dist)
        print 'Euc dist sensitivity,1-specificty = '+str(best_coord)


        print "trap area = "+str(area)

        import matplotlib.pyplot as plt
        plt.plot([0,1])
        #plt.scatter(fpr_list,tpr_list,color='r')
        plt.scatter(spec_list,sens_list,color='r')
        plt.plot([spec_list[youind]], [sens_list[youind]], marker='o', markersize=9, color="green")
        plt.plot([spec_list[concord_ind]], [sens_list[concord_ind]], marker='o', markersize=9, color="cyan")
        plt.plot([spec_list[index_of_union_index]], [sens_list[index_of_union_index]], marker='o', markersize=9, color="yellow")
        plt.plot([best_coord[1]], [best_coord[0]], marker='o', markersize=9, color="blue")
        plt.plot([0,best_coord[1]],[1,best_coord[0]])
        plt.plot([spec_list[self.one_percent_length]],[sens_list[self.one_percent_length]], marker='o', markersize=9, color="black")
        plt.ylabel('sensitivity')
        plt.xlabel('1-specificity')
        plt.title('Roc Curve Attempt ProteinInference')
        plt.legend(('Random Guess', 'Max Youdens Index', 'Max Concordance Prob', 'Min Index of Union', 'Min Euclidian Dist', 'Min Dist Line','1 Percent FDR','Sens vs 1-Spec'),
           shadow=True, loc=(.6, .05))
        if not pdf:
            plt.show()
        if pdf:
            pdf.savefig()
            # pp.close()


        youden_data = [self.all_data[x] for x in range(youind)]
        self.data_class.max_youdens_data = youden_data

        decoys = [x[0].identifier for x in youden_data if '#' in x[0].identifier]
        targets = [x[0].identifier for x in youden_data if '#' not in x[0].identifier]
        youdens_fdr = (len(decoys)*2)/float(len(targets))
        print 'FDR from Youdens Data = '+str(youdens_fdr)