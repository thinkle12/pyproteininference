#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:25:00 2017

@author: hinklet
"""

import collections
from Bio import SeqIO
from pyteomics import fasta, parser
from logging import getLogger
import re


class Digest(object):
    LIST_OF_DIGEST_TYPES = ["trypsin", "lysc"]
    AA_LIST = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]
    UNIPROT_STRS = "sp\||tr\|"
    UNIPROT_STR_REGEX = re.compile(UNIPROT_STRS)
    SP_STRING = "sp|"
    METHIONINE = "M"

    def __init__(self):
        pass


class InSilicoDigest(Digest):
    """
    The following class creates protein to peptide, peptide to protein, and swissprot protein mappings.
    These mappings are essential for GlpkGrouper as an InSilicoDigest object is input for GlpkGrouper

    The input is a fasta database, number of missed cleavages, as well as a digestion type ("trypsin").

    Further digestion types need to be added in the future other than just trypsin

    Exmample: Digest.insilicodigest.InSilicoDigest(database_path = "example_human_db.fasta", num_miss_cleavs=2, digest_type='trypsin')
    """

    def __init__(self, database_path, parameter_file_object, id_splitting=True):
        self.peptide_to_protein_dictionary = {}
        self.protein_to_peptide_dictionary = {}
        self.swiss_prot_protein_set = set()
        self.database_path = database_path
        self.num_miss_cleavs = parameter_file_object.missed_cleavages
        self.id_splitting = id_splitting
        self.reviewed_identifier_symbol = (
            parameter_file_object.reviewed_identifier_symbol
        )
        if parameter_file_object.digest_type in self.LIST_OF_DIGEST_TYPES:
            self.digest_type = parameter_file_object.digest_type
        else:
            raise ValueError(
                "digest_type must be equal to one of the following"
                + str(self.LIST_OF_DIGEST_TYPES)
                + " or... (List more digest types here in the future...)"
            )
        self.parameter_file_object = parameter_file_object
        self.logger = getLogger("protein_inference.in_silico_digest.InSilicoDigest")
        self.max_peptide_length = parameter_file_object.restrict_peptide_length

    def digest(self, proseq, miss_cleavage):
        peptides = []
        cut_sites = [0]
        if self.digest_type == "trypsin":
            for i in range(0, len(proseq) - 1):
                if proseq[i] == "K" and proseq[i + 1] != "P":
                    cut_sites.append(i + 1)
                elif proseq[i] == "R" and proseq[i + 1] != "P":
                    cut_sites.append(i + 1)

        if self.digest_type == "lysc":
            for i in range(0, len(proseq) - 1):
                if proseq[i] == "K" and proseq[i + 1] != "P":
                    cut_sites.append(i + 1)
            # Here write more code for other types of digest types....

        if cut_sites[-1] != len(proseq):
            cut_sites.append(len(proseq))

        if len(cut_sites) > 2:
            if miss_cleavage == 0:
                for j in range(0, len(cut_sites) - 1):
                    no_miss_cleave_pep = proseq[cut_sites[j] : cut_sites[j + 1]]
                    peptides.append(no_miss_cleave_pep)
                    # Account for N terminal Methionine Potential Cleavage
                    if j == 0 and proseq[cut_sites[0]] == "M":
                        peptides.append(no_miss_cleave_pep[1:])
                    else:
                        pass

            elif miss_cleavage == 1:
                for j in range(0, len(cut_sites) - 2):
                    peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]])
                    peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]])
                    # Account for N terminal Methionine Potential Cleavage
                    if j == 0 and proseq[cut_sites[0]] == "M":
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]][1:])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]][1:])

                # Account for N terminal Methionine Potential Cleavage
                if cut_sites[-2] == 0 and proseq[cut_sites[-2]] == "M":
                    peptides.append(proseq[cut_sites[-2] : cut_sites[-1]][1:])

                peptides.append(proseq[cut_sites[-2] : cut_sites[-1]])

            elif miss_cleavage == 2:
                for j in range(0, len(cut_sites) - 3):
                    peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]])
                    peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]])
                    peptides.append(proseq[cut_sites[j] : cut_sites[j + 3]])
                    # Account for N terminal Methionine Potential Cleavage
                    if j == 0 and proseq[cut_sites[0]] == "M":
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]][1:])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]][1:])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 3]][1:])

                # Account for N terminal Methionine Potential Cleavage
                if cut_sites[-3] == 0 and proseq[cut_sites[-3]] == "M":
                    peptides.append(proseq[cut_sites[-3] : cut_sites[-2]][1:])
                    peptides.append(proseq[cut_sites[-3] : cut_sites[-1]][1:])
                if cut_sites[-2] == 0 and proseq[cut_sites[-2]] == "M":
                    peptides.append(proseq[cut_sites[-2] : cut_sites[-1]][1:])

                peptides.append(proseq[cut_sites[-3] : cut_sites[-2]])
                peptides.append(proseq[cut_sites[-3] : cut_sites[-1]])
                peptides.append(proseq[cut_sites[-2] : cut_sites[-1]])

            elif miss_cleavage == 3:
                # If len cut sites is greater than 3... then we can do 3 missed cleavages
                if len(cut_sites) > 3:
                    for j in range(0, len(cut_sites) - 4):
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 3]])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 4]])
                        # Account for N terminal Methionine Potential Cleavage
                        if j == 0 and proseq[cut_sites[0]] == "M":
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]][1:])
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]][1:])
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 3]][1:])
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 4]][1:])

                    # Account for N terminal Methionine Potential Cleavage
                    if cut_sites[-3] == 0 and proseq[cut_sites[-3]] == "M":
                        peptides.append(proseq[cut_sites[-3] : cut_sites[-2]][1:])
                        peptides.append(proseq[cut_sites[-3] : cut_sites[-1]][1:])
                    if cut_sites[-2] == 0 and proseq[cut_sites[-2]] == "M":
                        peptides.append(proseq[cut_sites[-2] : cut_sites[-1]][1:])
                    if cut_sites[-4] == 0 and proseq[cut_sites[-4]] == "M":
                        peptides.append(proseq[cut_sites[-4] : cut_sites[-3]][1:])
                        peptides.append(proseq[cut_sites[-4] : cut_sites[-2]][1:])
                        peptides.append(proseq[cut_sites[-4] : cut_sites[-1]][1:])

                    peptides.append(proseq[cut_sites[-3] : cut_sites[-2]])
                    peptides.append(proseq[cut_sites[-3] : cut_sites[-1]])
                    peptides.append(proseq[cut_sites[-2] : cut_sites[-1]])
                    peptides.append(proseq[cut_sites[-4] : cut_sites[-1]])
                    peptides.append(proseq[cut_sites[-4] : cut_sites[-2]])
                    peptides.append(proseq[cut_sites[-4] : cut_sites[-3]])

                else:
                    # If len cut sites not greater than 3... then we do 2 MC
                    for j in range(0, len(cut_sites) - 3):
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]])
                        peptides.append(proseq[cut_sites[j] : cut_sites[j + 3]])
                        # Account for N terminal Methionine Potential Cleavage
                        if j == 0 and proseq[cut_sites[0]] == "M":
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 1]][1:])
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 2]][1:])
                            peptides.append(proseq[cut_sites[j] : cut_sites[j + 3]][1:])

                    # Account for N terminal Methionine Potential Cleavage
                    if cut_sites[-3] == 0 and proseq[cut_sites[-3]] == "M":
                        peptides.append(proseq[cut_sites[-3] : cut_sites[-2]][1:])
                        peptides.append(proseq[cut_sites[-3] : cut_sites[-1]][1:])
                    if cut_sites[-2] == 0 and proseq[cut_sites[-2]] == "M":
                        peptides.append(proseq[cut_sites[-2] : cut_sites[-1]][1:])

                    peptides.append(proseq[cut_sites[-3] : cut_sites[-2]])
                    peptides.append(proseq[cut_sites[-3] : cut_sites[-1]])
                    peptides.append(proseq[cut_sites[-2] : cut_sites[-1]])

        else:  # there is no tryptic site in the protein sequence
            peptides.append(proseq)
            if proseq[0] == "M":
                peptides.append(proseq[1:])

        new_peptides = []
        if any("X" in x for x in peptides):
            for peps in peptides:
                if "X" in peps:
                    x_index = peps.index("X")
                    for aa in self.AA_LIST:
                        peptide_list = list(peps)
                        peptide_list[x_index] = aa
                        new_pep = "".join(peptide_list)
                        new_peptides.append(new_pep)

        peptides = peptides + new_peptides
        # peptides = [x for x in peptides if len(x)>=self.max_peptide_length]
        return peptides

    def digest_fasta_database(self):
        handle = SeqIO.parse(self.database_path, "fasta")

        self.logger.info("Starting Standard Digest...")
        pep_dict = {}
        prot_dict = {}
        sp_set = set()

        for record in handle:
            if self.id_splitting == True:
                identifier_stripped = self.UNIPROT_STR_REGEX.sub("", record.id)
            if self.id_splitting == False:
                identifier_stripped = record.id

            if record.id.startswith(self.SP_STRING):
                sp_set.add(identifier_stripped)

            proseq = str(record.seq)
            peptide_list = self.digest(proseq, self.num_miss_cleavs)
            prot_dict[identifier_stripped] = set(peptide_list)
            for peptide in peptide_list:
                pep_dict.setdefault(peptide, set()).add(identifier_stripped)

        self.logger.info(
            "Digest finished, peptide and protein dictionaries created based on the provided database"
        )

        self.swiss_prot_protein_set = sp_set
        self.peptide_to_protein_dictionary = pep_dict
        self.protein_to_peptide_dictionary = prot_dict

        self.logger.info("Standard Digest Finished...")


class PyteomicsDigest(Digest):
    """
    The following class creates protein to peptide, peptide to protein, and swissprot protein mappings.
    These mappings are essential for GlpkGrouper as an InSilicoDigest object is input for GlpkGrouper

    The input is a fasta database, number of missed cleavages, as well as a digestion type ("trypsin").

    Further digestion types need to be added in the future other than just trypsin

    """

    def __init__(self, database_path, parameter_file_object, id_splitting=True):
        self.peptide_to_protein_dictionary = {}
        self.protein_to_peptide_dictionary = {}
        self.swiss_prot_protein_set = set()
        self.database_path = database_path
        self.num_miss_cleavs = parameter_file_object.missed_cleavages
        self.id_splitting = id_splitting
        self.reviewed_identifier_symbol = (
            parameter_file_object.reviewed_identifier_symbol
        )
        if parameter_file_object.digest_type in self.LIST_OF_DIGEST_TYPES:
            self.digest_type = parameter_file_object.digest_type
        else:
            raise ValueError(
                "digest_type must be equal to one of the following"
                + str(self.LIST_OF_DIGEST_TYPES)
                + " or... (List more digest types here in the future...)"
            )
        self.logger = getLogger("protein_inference.in_silico_digest.PyteomicsDigest")
        self.max_peptide_length = parameter_file_object.restrict_peptide_length

    def digest_fasta_database(self):
        self.logger.info("Starting Pyteomics Digest...")
        pep_dict = {}
        prot_dict = {}
        sp_set = set()

        for description, sequence in fasta.read(self.database_path):
            new_peptides = parser.cleave(
                sequence,
                parser.expasy_rules[self.digest_type],
                self.num_miss_cleavs,
                min_length=self.max_peptide_length,
            )

            # Hopefully this splitting works...
            # IDK how robust this is...
            identifier = description.split(" ")[0]

            # Handle ID Splitting...
            if self.id_splitting == True:
                identifier_stripped = self.UNIPROT_STR_REGEX.sub("", identifier)
            if self.id_splitting == False:
                identifier_stripped = identifier

            # TODO instead of using [:3] we should replace sp| with nothing
            # If SP add to
            if identifier.startswith(self.SP_STRING):
                sp_set.add(identifier_stripped)
            prot_dict[identifier_stripped] = new_peptides
            met_cleaved_peps = set()
            for peptide in new_peptides:
                pep_dict.setdefault(peptide, []).append(identifier_stripped)
                # Need to account for potential N-term Methionine Cleavage
                if sequence.startswith(peptide) and peptide.startswith(self.METHIONINE):
                    # If our sequence starts with the current peptide... and our current peptide starts with methionine...
                    # Then we remove the methionine from the peptide and add it to our dicts...
                    methionine_cleaved_peptide = peptide[1:]
                    met_cleaved_peps.add(methionine_cleaved_peptide)
            for met_peps in met_cleaved_peps:
                pep_dict.setdefault(met_peps, []).append(identifier_stripped)
                prot_dict[identifier_stripped].add(met_peps)

        self.swiss_prot_protein_set = sp_set
        self.peptide_to_protein_dictionary = pep_dict
        self.protein_to_peptide_dictionary = prot_dict

        self.logger.info("Pyteomics Digest Finished...")
