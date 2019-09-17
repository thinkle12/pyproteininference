import csv
import uuid
import tempfile
import os
import subprocess

class Qvality(object):

    def __init__(self):
        None

class CalculateQandPepValues(Qvality):
    """
    Class calculates Q values and Pep Values using the command line too Qvality, which is distributed with Percolator

    Input is a DataStore object

    Example: protein_inference.qvality.CalculateQandPepValues(data_class = data)

    """
    def __init__(self, data_class):
        self.grouped_scored_data = data_class.grouped_scored_proteins
        self.data_class = data_class
        self.uuid = uuid.uuid4()
        self.uuid_tag = str(self.uuid)
        self.temp_dir = tempfile.gettempdir()
        self.target_score_file = os.path.join(self.temp_dir,'target_scores-'+ self.uuid_tag)
        self.decoy_score_file = os.path.join(self.temp_dir,'decoy_scores-'+ self.uuid_tag)
        self.qvality_output_filename = os.path.join(self.temp_dir,'qvality_output-'+ self.uuid_tag)

    def execute(self):
        # Need to generate UUIDs for the filenames... because we are going to delete all of them...
        targets = []
        decoys = []
        for prots in self.grouped_scored_data:
            lead_prot = prots[0]
            if '#' in lead_prot.identifier:
                decoys.append(str(lead_prot.score))
            else:
                targets.append(str(lead_prot.score))


        with open(self.target_score_file, "wb") as f:
            for t in targets:
                f.write(t+'\n')
        with open(self.decoy_score_file, "wb") as f:
            for d in decoys:
                f.write(d+'\n')
        p = subprocess.Popen('qvality '+self.target_score_file+' '+self.decoy_score_file+' -o '+self.qvality_output_filename,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        output = p.communicate()

        print('Start Command line Stdout')
        print(output[0])
        print('End Command line Stdout')
        print('Start Command line Stderr')
        print(output[1])
        print('End Command line Stderr')
        # Create a file that has all target scores, and a file that has all decoy scores....

        qo = open(self.qvality_output_filename, 'r')
        qo = qo.read()
        qo = qo.split('\n')
        qvality_output = [x.split('\t') for x in qo]
        self.qvality_output = qvality_output
        self.data_class.qvality_output = qvality_output