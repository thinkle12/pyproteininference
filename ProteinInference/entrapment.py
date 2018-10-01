import ProteinInference
import matplotlib
#
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import math
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
import scipy.stats
import sklearn.metrics

class Entrapment(object):

    def __init__(self):
        None


class GeneratePlot(Entrapment):
    """
    This class is for experienced users only. The purpose of this class is to generate plots of Entrapment FDR vs Decoy FDR
    in an attempt to test to see if entrapment vs decoy fdr is 1 to 1.

    The class requires a DataStore object, an entrapment database (proteins known not to be in the sample),
    a True database (proteins known to be in the sample), a search_id for tagging purposes, a pdf option to save plots as a pdf,
    a picked/non picked boolean, q value restriction, and p value restriction (for labeling the plot),
    and finally an "other_database" if neccessary to add to the entrapment database

    This class produces plots of Entrapment FDR vs Decoy FDR for the purpose of benchmarking the Protein Inference tool

    This class requires the search being used to be from a benchmark dataset of KNOWN protein content.
    """

    def __init__(self, data_class, entrapment_db, true_db, search_id, pdf=None, picked=None, qvr=None, pvr=None, other_database = None, search_database = None):
        self.data_class = data_class
        self.entrapment_db = entrapment_db
        self.true_db = true_db
        self.grouped_scored_data = data_class.grouped_scored_proteins
        self.search_id = search_id
        if picked:
            self.pick = 'Yes Picker'
        else:
            self.pick = 'No Picker'
        self.qvr = qvr
        self.pvr = pvr

        if other_database:
            self.other_database = other_database

        if not other_database:
            self.other_database = None


        if search_database:
            self.search_database = self.search_database

        if not search_database:
            self.search_database = None


        self.average_entrapment_fdr = None
        self.average_decoy_fdr = None
        self.pdf = pdf

        fdr1 = ProteinInference.fdrcalc.PureDecoyFdr(data_class=data_class, false_discovery_rate=0)
        fdr1.execute()
        self.decoy_fdr_calc_list = fdr1.fdr_list
        efdr1 = ProteinInference.fdrcalc.EntrapFdr(data_class=data_class, entrapment_database=self.entrapment_db,
                                                   true_database=self.true_db,other_database=self.other_database, false_discovery_rate=0)
        efdr1.execute()
        self.entrapment_fdr_calc_list = efdr1.fdr_list

    def execute(self):



        true_handle = SeqIO.parse(self.true_db, 'fasta')

        # Get the true proteins minus the decoys...
        true_proteins = []
        for records in true_handle:
            if '#' not in records.id:
                true_proteins.append(records.id)

        # search_handle = SeqIO.parse(self.search_database, 'fasta')
        #
        # search_proteins = []
        # for records in search_handle:
        #     search_proteins.append(records.id)
        #
        # decoy_search_proteins = [x for x in search_proteins if "#" in x]

        entrapment_handle = SeqIO.parse(self.entrapment_db, 'fasta')

        # Get the entrapment proteins minus the decoys...
        entrapment_proteins = []
        for records in entrapment_handle:
            if '#' not in records.id:
                entrapment_proteins.append(records.id)

        # Get the proteins from the search...
        protein_data = [x[0].identifier for x in self.grouped_scored_data]

        # Add in proteins found in the search that arent decoys and arent already in entrapment proteins...
        for prots in protein_data:
            if '#' not in prots and prots not in entrapment_proteins:
                entrapment_proteins.append(prots)

        print 'Number of Grouped Proteins = '+str(len(protein_data))


        # This is with NO Decoys...
        true_to_entrap_ratio = float(len(true_proteins)+1000)/float(len(entrapment_proteins))
        # decoy_to_entrap_ratio = float(len(entrapment_proteins))/float(len(decoy_search_proteins))
        print str(true_to_entrap_ratio)
        # print str(decoy_to_entrap_ratio)

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

            dfdr = float(len(decoys))/float(len(protein_data))
            efdr = (float(len(false_true_positives))/float(len(protein_data)))

            decoy_fdr.append(dfdr)
            entrapment_fdr.append(efdr)

        self.decoy_fdr = decoy_fdr
        self.entrapment_fdr = entrapment_fdr


        # self.decoy_fdr_list.append(decoy_fdr)
        # self.entrapment_fdr_list.append(entrapment_fdr)

        print 'Number of entraped proteins = '+str(len(false_true_positives))
        print 'Number of decoy proteins = '+str(len(decoys))

        # #Plot 1
        # plt.plot([0, 1])
        # plt.plot([0, 0.6666666666666666],'--')
        # plt.plot([0, 1.5],'--')
        # plt.scatter(decoy_fdr,entrapment_fdr)
        # plt.xlim([0, .1])
        # plt.ylim([0, .1])
        # plt.xlabel('Decoy FDR')
        # plt.ylabel('Entrapment FDR')
        # plt.title(self.data_class.score_method+'\n'+'Protein set '+str(db)+' used'+' with SearchID: '+self.search_id)
        # plt.show()
        # plt.close()

        efdr05 = ProteinInference.fdrcalc.EntrapFdr(data_class=self.data_class, entrapment_database=self.entrapment_db,
                                                    true_database=self.true_db,other_database=self.other_database,
                                                    false_discovery_rate=.05)
        efdr05.execute()

        efdr05.restricted_proteins

        self.prots_that_pass = efdr05.restricted_proteins

        # Plot 2
        plt.plot([0, 1])
        plt.plot([0, 0.6666666666666666], '--')
        plt.plot([0, 1.5], '--')
        plt.plot(self.decoy_fdr_calc_list, self.entrapment_fdr_calc_list,'-o')
        plt.xlim([0.00001, .5])
        plt.ylim([0.00001, .5])
        plt.xlabel('Decoy FDR')
        plt.ylabel('Entrapment FDR')
        plt.title(self.data_class.score_method + ' with '+str(len(efdr05.restricted_proteins))+' passing proteins' + '\n' + 'SearchID: ' + self.search_id + ', ' +self.pick + ', qvr = '+str(self.qvr)+ ', pvr = '+str(self.pvr))
        if self.pdf:
            self.pdf.savefig()
        plt.show()
        plt.close()

        full_vertical_distance_list = []
        for j in range(len(self.decoy_fdr_calc_list)):
            full_vertical_distance = self.entrapment_fdr_calc_list[j]-self.decoy_fdr_calc_list[j]
            full_vertical_distance_list.append(abs(full_vertical_distance))

        full_vertical_distance_mean = numpy.mean(full_vertical_distance_list)
        self.full_vertical_distance_mean = full_vertical_distance_mean

        full_perp_dist_list = []
        p1 = numpy.array([0, 0])
        p2 = numpy.array([1, 1])

        for i in range(len(self.decoy_fdr_calc_list)):
            p3 = numpy.array([self.decoy_fdr_calc_list[i], self.entrapment_fdr_calc_list[i]])

            full_perp_distance = numpy.linalg.norm(numpy.cross(p2-p1, p1-p3))/numpy.linalg.norm(p2-p1)
            full_perp_dist_list.append(full_perp_distance)

        full_perp_dist_mean = numpy.mean(full_perp_dist_list)
        self.full_perp_dist_mean = full_perp_dist_mean

        vertical_distance_list = []
        for j in range(len(self.decoy_fdr_calc_list)):
            if self.decoy_fdr_calc_list[j]<.05:
                vertical_distance = self.entrapment_fdr_calc_list[j] - self.decoy_fdr_calc_list[j]
                vertical_distance_list.append(abs(vertical_distance))

        vertical_distance_mean = numpy.mean(vertical_distance_list)
        self.vertical_distance_mean = vertical_distance_mean

        perp_dist_list = []
        p1 = numpy.array([0, 0])
        p2 = numpy.array([1, 1])

        for i in range(len(self.decoy_fdr_calc_list)):
            if self.decoy_fdr_calc_list[i]<.05:
                p3 = numpy.array([self.decoy_fdr_calc_list[i], self.entrapment_fdr_calc_list[i]])

                perp_distance = numpy.linalg.norm(numpy.cross(p2 - p1, p1 - p3)) / numpy.linalg.norm(p2 - p1)
                perp_dist_list.append(perp_distance)

        perp_dist_mean = numpy.mean(perp_dist_list)
        self.perp_dist_mean = perp_dist_mean



        xlist = [x for x in self.decoy_fdr_calc_list if x<.05]

        data_array = numpy.array([self.entrapment_fdr_calc_list[x] for x in range(len(self.decoy_fdr_calc_list)) if self.decoy_fdr_calc_list[x]<.05])

        test_array = numpy.array([x for x in xlist])

        self.data_array = data_array
        self.test_array = test_array


        pr = scipy.stats.pearsonr(test_array,data_array)

        pr = pr[0]

        self.pr = pr
        #
        # pearson_list = []
        # for k in range(len(data_array)):
        #     pearson = scipy.stats.pearsonr(data_array[k],test_array[k])
        #     correl = pearson[0]
        #     pearson_list.append(correl)
        #
        # self.pearson = pearson_list

        # if len(self.decoy_fdr_list) > 1:
        #     max_len = max([len(x) for x in self.decoy_fdr_list])
        #     for kk in range(len(self.decoy_fdr_list)):
        #         if len(self.decoy_fdr_list[kk]) < max_len:
        #             for j in range(10000):
        #                 self.decoy_fdr_list[kk].append(numpy.nan)
        #                 if len(self.decoy_fdr_list[kk]) >= max_len:
        #                     break
        #
        #     decoy_fdr_array = numpy.array(self.decoy_fdr_list)
        #     avg_decoy_fdr = numpy.nanmean(decoy_fdr_array,axis=0).tolist
        #
        # if len(self.entrapment_fdr_list) > 1:
        #     max_len = max([len(x) for x in self.entrapment_fdr_list])
        #     for kk in range(len(self.entrapment_fdr_list)):
        #         if len(self.entrapment_fdr_list[kk]) < max_len:
        #             for j in range(10000):
        #                 self.entrapment_fdr_list[kk].append(numpy.nan)
        #                 if len(self.entrapment_fdr_list[kk]) >= max_len:
        #                     break
        #
        #     entrap_fdr_array = numpy.array(self.entrapment_fdr_list)
        #     avg_entrap_fdr = numpy.nanmean(entrap_fdr_array,axis=0).tolist
        #
        #
        #
        # if len(self.entrapment_fdr_list) > 1 and len(self.decoy_fdr_list) > 1:
        #     self.average_decoy_fdr = avg_decoy_fdr
        #     self.average_entrapment_fdr = avg_entrap_fdr
        #
        #     plt.plot([0, 1])
        #     plt.plot([0, 0.6666666666666666], '--')
        #     plt.plot([0, 1.5], '--')
        #     plt.scatter([float(x) for x in avg_decoy_fdr()], [float(x) for x in avg_entrap_fdr()])
        #     plt.xlim([0, .1])
        #     plt.ylim([0, .1])
        #     plt.xlabel('Decoy FDR')
        #     plt.ylabel('Entrapment FDR')
        #     plt.title(self.data_class.score_method + '\n' + 'Protein set ' + str(db) + ' used '+'averaged '+str(len(self.decoy_fdr_list))+' times')
        #     if self.pdf:
        #         self.pdf.savefig()
        #     # plt.show()
        #     # plt.close()
        #     plt.close()
        #
        #
        #
