#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 12:52:14 2017

@author: hinklet
"""

import os
from protein_inference.physical import Psm
from protein_inference.datastore import DataStore
import csv
import itertools
from logging import getLogger


class Reader(object):
    """
    Main Reader Class which is parent to all reader subclasses
    """

    MAX_ALLOWED_ALTERNATIVE_PROTEINS = 50

    def __init__(
        self, target_file=None, decoy_file=None, combined_files=None, directory=None
    ):
        self.target_file = target_file
        self.decoy_file = decoy_file
        self.combined_files = combined_files
        self.directory = directory

    def remap(self, fieldnames):
        price_count = itertools.count(1)
        return [
            "alternative_protein_{}".format(next(price_count)) if f == "" else f
            for f in fieldnames
        ]

    def _validate_input(self):
        if (
            not self.target_file
            and not self.decoy_file
            and not self.combined_files
            and not self.directory
        ):
            raise ValueError(
                "No input provided, please supply either target and decoy files, combined files, or a directory of combined target/decoy files"
            )
        else:
            pass

    @classmethod
    def _fix_alternative_proteins(
        cls,
        append_alt_from_db,
        identifiers_sorted,
        max_proteins,
        psm,
        parameter_file_object,
    ):
        # If we are appending alternative proteins from the db
        if append_alt_from_db:
            # Loop over the Identifiers from the DB These are identifiers that contain the current peptide
            for alt_proteins in identifiers_sorted[:max_proteins]:
                # If the identifier is not already in possible proteins and if then len of poss prot is less than the max...
                # Then append
                if (
                    alt_proteins not in psm.possible_proteins
                    and len(psm.possible_proteins) < max_proteins
                ):
                    psm.possible_proteins.append(alt_proteins)
        # Next if the len of possible proteins is greater than max then restrict the list length...
        if len(psm.possible_proteins) > max_proteins:
            psm.possible_proteins = [
                psm.possible_proteins[x] for x in range(max_proteins)
            ]
        else:
            pass

        # If no inference only select first poss protein
        if parameter_file_object.inference_type == "none":
            psm.possible_proteins = [psm.possible_proteins[0]]

        return psm


class PercolatorReader(Reader):
    """
    The following class takes a percolator target file and a percolator decoy file and creates standard protein_inference.physical.Psm() objects.
    This reader class is used as input for the DataStore class which is used in all further protein_inference Classes

    Example: protein_inference.reader.PercolatorReader(target_file = "example_target.txt", decoy_file = "example_decoy.txt")

    Percolator Output is formatted as follows:
    with each entry being tabbed delimited
    PSMId	 score	q-value	posterior_error_prob	peptide	proteinIds
    116108.15139.15139.6.dta	 3.44016	 0.000479928	7.60258e-10	K.MVVSMTLGLHPWIANIDDTQYLAAK.R	CNDP1_HUMAN|Q96KN2	B4E180_HUMAN|B4E180	A8K1K1_HUMAN|A8K1K1	J3KRP0_HUMAN|J3KRP0

    """

    def __init__(
        self,
        digest_class,
        parameter_file_object,
        append_alt_from_db=False,
        target_file=None,
        decoy_file=None,
        combined_files=None,
        directory=None,
    ):
        self.target_file = target_file
        self.decoy_file = decoy_file
        self.combined_files = combined_files
        self.directory = directory
        self._validate_input()
        # Define Indicies based on input

        self.psmid_index = 0
        self.perc_score_index = 1
        self.q_value_index = 2
        self.posterior_error_prob_index = 3
        self.peptide_index = 4
        self.proteinIDs_index = 5
        self.psms = None
        self.search_id = None
        self.digest_class = digest_class
        self.append_alt_from_db = append_alt_from_db

        self.parameter_file_object = parameter_file_object
        self.logger = getLogger("protein_inference.reader.PercolatorReader.read_psms")

    def read_psms(self):
        # Read in and split by line
        if self.target_file and self.decoy_file:
            # If target_file is a list... read them all in and concatenate...
            if isinstance(self.target_file, (list,)):
                all_target = []
                for t_files in self.target_file:
                    self.logger.info(t_files)
                    ptarg = []
                    with open(t_files, "r") as perc_target_file:
                        spamreader = csv.reader(perc_target_file, delimiter="\t")
                        for row in spamreader:
                            ptarg.append(row)
                    del ptarg[0]
                    all_target = all_target + ptarg
            elif self.target_file:
                # If not just read the file...
                ptarg = []
                with open(self.target_file, "r") as perc_target_file:
                    spamreader = csv.reader(perc_target_file, delimiter="\t")
                    for row in spamreader:
                        ptarg.append(row)
                del ptarg[0]
                all_target = ptarg

            # Repeat for decoy file
            if isinstance(self.decoy_file, (list,)):
                all_decoy = []
                for d_files in self.decoy_file:
                    self.logger.info(d_files)
                    pdec = []
                    with open(d_files, "r") as perc_decoy_file:
                        spamreader = csv.reader(perc_decoy_file, delimiter="\t")
                        for row in spamreader:
                            pdec.append(row)
                    del pdec[0]
                    all_decoy = all_decoy + pdec
            elif self.decoy_file:
                pdec = []
                with open(self.decoy_file, "r") as perc_decoy_file:
                    spamreader = csv.reader(perc_decoy_file, delimiter="\t")
                    for row in spamreader:
                        pdec.append(row)
                del pdec[0]
                all_decoy = pdec

            # Combine the lists
            perc_all = all_target + all_decoy

        elif self.combined_files:
            if isinstance(self.combined_files, (list,)):
                all = []
                for f in self.combined_files:
                    self.logger.info(f)
                    p = []
                    with open(f, "r") as perc_files:
                        spamreader = csv.reader(perc_files, delimiter="\t")
                        for row in spamreader:
                            p.append(row)
                    del p[0]
                    all = all + p
            elif self.combined_files:
                # If not just read the file...
                p = []
                with open(self.combined_files, "r") as perc_files:
                    spamreader = csv.reader(perc_files, delimiter="\t")
                    for row in spamreader:
                        p.append(row)
                del p[0]
                all = p
            perc_all = all

        elif self.directory:

            all_files = os.listdir(self.directory)
            all = []
            for files in all_files:
                self.logger.info(files)
                p = []
                with open(files, "r") as perc_file:
                    spamreader = csv.reader(perc_file, delimiter="\t")
                    for row in spamreader:
                        p.append(row)
                del p[0]
                all = all + p
            perc_all = all

        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary

        perc_all_filtered = []
        for psms in perc_all:
            try:
                float(psms[self.posterior_error_prob_index])
                perc_all_filtered.append(psms)
            except ValueError as e:
                pass

        # Filter by pep
        perc_all = sorted(
            perc_all_filtered,
            key=lambda x: float(x[self.posterior_error_prob_index]),
            reverse=False,
        )

        # TODO
        # TRY TO GET PERC_ALL AS A GENERATOR
        # Can do this... just give the option to feed a combined file... and loop over all the files...
        # Hmm... problem is it still needs to be sorted by perc score....

        list_of_psm_objects = []
        peptide_tracker = set()
        all_sp_proteins = set(self.digest_class.swiss_prot_protein_set)
        # We only want to get unique peptides... using all messes up scoring...
        # Create Psm objects with the identifier, percscore, qvalue, pepvalue, and possible proteins...

        # TODO
        # make this for loop a generator...

        self.logger.info("Length of PSM Data: {}".format(len(perc_all)))
        for psm_info in perc_all:
            current_peptide = psm_info[self.peptide_index]
            # Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier=current_peptide)
                # Add all the attributes
                p.percscore = float(psm_info[self.perc_score_index])
                p.qvalue = float(psm_info[self.q_value_index])
                p.pepvalue = float(psm_info[self.posterior_error_prob_index])
                if self.parameter_file_object.inference_type == "none":
                    poss_proteins = [psm_info[self.proteinIDs_index]]
                else:
                    poss_proteins = list(
                        set(
                            psm_info[
                                self.proteinIDs_index : self.proteinIDs_index
                                + self.MAX_ALLOWED_ALTERNATIVE_PROTEINS
                            ]
                        )
                    )
                p.possible_proteins = (
                    poss_proteins  # Restrict to 50 total possible proteins...
                )
                p.psm_id = psm_info[self.psmid_index]

                # Split peptide if flanking
                current_peptide = Psm.split_peptide(peptide_string=current_peptide)

                if not current_peptide.isupper() or not current_peptide.isalpha():
                    # If we have mods remove them...
                    peptide_string = current_peptide.upper()
                    stripped_peptide = Psm.remove_peptide_mods(peptide_string)
                    current_peptide = stripped_peptide

                # Add the other possible_proteins from insilicodigest here...
                try:
                    current_alt_proteins = list(
                        peptide_to_protein_dictionary[current_peptide]
                    )  # This peptide needs to be scrubbed of Mods...
                except KeyError:
                    current_alt_proteins = []
                    self.logger.warning(
                        "Peptide {} was not found in the supplied DB with the following proteins {}".format(
                            current_peptide, ";".join(p.possible_proteins)
                        )
                    )
                    for poss_prot in p.possible_proteins:
                        self.digest_class.peptide_to_protein_dictionary.setdefault(
                            current_peptide, set()
                        ).add(poss_prot)
                        self.digest_class.protein_to_peptide_dictionary.setdefault(
                            poss_prot, set()
                        ).add(current_peptide)
                        self.logger.info(
                            "Adding Peptide {} and Protein {} to Digest dictionaries".format(
                                current_peptide, poss_prot
                            )
                        )

                # Sort Alt Proteins by Swissprot then Trembl...
                identifiers_sorted = DataStore.sort_protein_strings(
                    protein_string_list=current_alt_proteins,
                    sp_proteins=all_sp_proteins,
                    decoy_symbol=self.parameter_file_object.decoy_symbol,
                )

                # Restrict to 50 possible proteins
                p = self._fix_alternative_proteins(
                    append_alt_from_db=self.append_alt_from_db,
                    identifiers_sorted=identifiers_sorted,
                    max_proteins=self.MAX_ALLOWED_ALTERNATIVE_PROTEINS,
                    psm=p,
                    parameter_file_object=self.parameter_file_object,
                )

                # Remove blank alt proteins
                p.possible_proteins = [x for x in p.possible_proteins if x != ""]

                list_of_psm_objects.append(p)
                peptide_tracker.add(current_peptide)

        self.psms = list_of_psm_objects

        self.logger.info("Length of PSM Data: {}".format(len(self.psms)))

        # return perc


class ProteologicPostSearchReader(Reader):
    """
    Potential Future class to read in from another source.
    Potentially from a database, another Psm scoring source, or potentially from Percolator XML source.
    This Method will be used to read from post processing proteologic logical object... which can either be LDA or Percolator Results...
    Essentially it will just be a searchID and an LDA/Percolator ID
    """

    def __init__(
        self,
        proteologic_object,
        search_id,
        postsearch_id,
        digest_class,
        parameter_file_object,
        append_alt_from_db=False,
    ):
        self.proteologic_object = proteologic_object
        self.search_id = search_id
        self.postsearch_id = postsearch_id

        self.psms = None
        self.digest_class = digest_class
        self.append_alt_from_db = append_alt_from_db

        self.parameter_file_object = parameter_file_object
        self.logger = getLogger(
            "protein_inference.reader.ProteologicPostSearchReader.read_psms"
        )

    def read_psms(self):
        self.logger.info("Reading in data from Proteologic...")
        if isinstance(self.proteologic_object, (list,)):
            list_of_psms = []
            for p_objs in self.proteologic_object:
                for psms in p_objs.physical_object.psm_sets:
                    list_of_psms.append(psms)
        else:
            list_of_psms = self.proteologic_object.physical_object.psm_sets

        # Sort this by posterior error prob...
        list_of_psms = sorted(list_of_psms, key=lambda x: float(x.psm_filter.pepvalue))

        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary

        list_of_psm_objects = []
        peptide_tracker = set()
        all_sp_proteins = set(self.digest_class.swiss_prot_protein_set)
        # Peptide tracker is used because we only want UNIQUE peptides...
        # The data is sorted by percolator score... or at least it should be...
        # Or sorted by posterior error probability

        # TODO, We may need to try the quant thing on ALL PSMs, not just unique Peptides...
        for peps in list_of_psms:
            self.peps = peps
            current_peptide = peps.peptide.sequence
            # Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier=current_peptide)
                # Add all the attributes
                p.percscore = float(0)  # Will be stored in table in future I think...
                p.qvalue = float(peps.psm_filter.q_value)
                p.pepvalue = float(peps.psm_filter.pepvalue)
                if peps.peptide.protein not in peps.alternative_proteins:
                    p.possible_proteins = [
                        peps.peptide.protein
                    ] + peps.alternative_proteins
                else:
                    p.possible_proteins = peps.alternative_proteins

                p.possible_proteins = list(filter(None, p.possible_proteins))
                p.psm_id = peps.spectrum.spectrum_identifier

                # Split peptide if flanking
                current_peptide = Psm.split_peptide(peptide_string=current_peptide)

                if not current_peptide.isupper() or not current_peptide.isalpha():
                    # If we have mods remove them...
                    peptide_string = current_peptide.upper()
                    stripped_peptide = Psm.remove_peptide_mods(peptide_string)
                    current_peptide = stripped_peptide

                # Add the other possible_proteins from insilicodigest here...
                try:
                    current_alt_proteins = list(
                        peptide_to_protein_dictionary[current_peptide]
                    )  # This peptide needs to be scrubbed of Mods...
                except KeyError:
                    current_alt_proteins = []
                    self.logger.warning(
                        "Peptide {} was not found in the supplied DB with the following proteins {}".format(
                            current_peptide, ";".join(p.possible_proteins)
                        )
                    )
                    for poss_prot in p.possible_proteins:
                        self.digest_class.peptide_to_protein_dictionary.setdefault(
                            current_peptide, set()
                        ).add(poss_prot)
                        self.digest_class.protein_to_peptide_dictionary.setdefault(
                            poss_prot, set()
                        ).add(current_peptide)
                        self.logger.info(
                            "Adding Peptide {} and Protein {} to Digest dictionaries".format(
                                current_peptide, poss_prot
                            )
                        )

                # Sort Alt Proteins by Swissprot then Trembl...
                identifiers_sorted = DataStore.sort_protein_strings(
                    protein_string_list=current_alt_proteins,
                    sp_proteins=all_sp_proteins,
                    decoy_symbol=self.parameter_file_object.decoy_symbol,
                )

                # Restrict to 50 possible proteins... and append alt proteins from db
                p = self._fix_alternative_proteins(
                    append_alt_from_db=self.append_alt_from_db,
                    identifiers_sorted=identifiers_sorted,
                    max_proteins=self.MAX_ALLOWED_ALTERNATIVE_PROTEINS,
                    psm=p,
                    parameter_file_object=self.parameter_file_object,
                )

                list_of_psm_objects.append(p)
                peptide_tracker.add(current_peptide)

            # TODO, Here we will Keep track of a large set of all PSMs..
            # We will keep track of this new PSMid -> search_id.peptide_id

        self.psms = list_of_psm_objects
        self.logger.info("Finished reading in data from Proteologic...")


class GenericReader(Reader):
    """
    The following class takes a percolator like target file and a percolator like decoy file and creates standard protein_inference.physical.Psm() objects.
    This reader class is used as input for the DataStore class which is used in all further protein_inference Classes

    Example: protein_inference.reader.GenericReader(target_file = "example_target.txt", decoy_file = "example_decoy.txt")

    Percolator Output is formatted as follows:
    with each entry being tabbed delimited
    PSMId	 score	q-value	posterior_error_prob	peptide	proteinIds
    116108.15139.15139.6.dta	 3.44016	 0.000479928	7.60258e-10	K.MVVSMTLGLHPWIANIDDTQYLAAK.R	CNDP1_HUMAN|Q96KN2	B4E180_HUMAN|B4E180	A8K1K1_HUMAN|A8K1K1	J3KRP0_HUMAN|J3KRP0

    """

    PSMID = "PSMId"
    SCORE = "score"
    Q_VALUE = "q-value"
    POSTERIOR_ERROR_PROB = "posterior_error_prob"
    PEPTIDE = "peptide"
    PROTEIN_IDS = "proteinIds"
    EXTRA_PROTEIN_IDS = "alternative_protein_{}"

    def __init__(
        self,
        digest_class,
        parameter_file_object,
        append_alt_from_db=False,
        target_file=None,
        decoy_file=None,
        combined_files=None,
        directory=None,
    ):
        self.target_file = target_file
        self.decoy_file = decoy_file
        self.combined_files = combined_files
        self.directory = directory
        self._validate_input()

        self.psms = None
        self.search_id = None
        self.digest_class = digest_class
        self.load_custom_score = False

        self.append_alt_from_db = append_alt_from_db

        self.parameter_file_object = parameter_file_object
        self.scoring_variable = parameter_file_object.score

        self.logger = getLogger("protein_inference.reader.GenericReader.read_psms")

        if (
            self.scoring_variable != self.Q_VALUE
            and self.scoring_variable != self.POSTERIOR_ERROR_PROB
        ):
            self.load_custom_score = True
            self.logger.info(
                "Pulling custom column based on parameter file input for score, Column: {}".format(
                    self.scoring_variable
                )
            )
        else:
            self.logger.info(
                "Pulling no custom columns based on parameter file input for score, using standard Column: {}".format(
                    self.scoring_variable
                )
            )

        self.MAX_ALTERNATIVE_PROTEIN_COLUMN_NAMES = [
            self.EXTRA_PROTEIN_IDS.format(x)
            for x in range(1, self.MAX_ALLOWED_ALTERNATIVE_PROTEINS + 1)
        ]

        # If we select to not run inference at all
        if self.parameter_file_object.inference_type == "none":
            # Only allow 1 Protein per PSM
            self.MAX_ALLOWED_ALTERNATIVE_PROTEINS = 1

    def read_psms(self):
        self.logger.info("Reading in Input Files using Generic Reader...")
        # Read in and split by line
        # If target_file is a list... read them all in and concatenate...
        if self.target_file and self.decoy_file:
            if isinstance(self.target_file, (list,)):
                all_target = []
                for t_files in self.target_file:
                    ptarg = []
                    with open(t_files, "r") as psm_target_file:
                        self.logger.info(t_files)
                        spamreader = csv.reader(psm_target_file, delimiter="\t")
                        fieldnames = self.remap(next(spamreader))
                        for row in spamreader:
                            ptarg.append(dict(zip(fieldnames, row)))
                    all_target = all_target + ptarg
            else:
                # If not just read the file...
                ptarg = []
                with open(self.target_file, "r") as psm_target_file:
                    self.logger.info(self.target_file)
                    spamreader = csv.reader(psm_target_file, delimiter="\t")
                    fieldnames = self.remap(next(spamreader))
                    for row in spamreader:
                        ptarg.append(dict(zip(fieldnames, row)))
                all_target = ptarg

            # Repeat for decoy file
            if isinstance(self.decoy_file, (list,)):
                all_decoy = []
                for d_files in self.decoy_file:
                    pdec = []
                    with open(d_files, "r") as psm_decoy_file:
                        self.logger.info(d_files)
                        spamreader = csv.reader(psm_decoy_file, delimiter="\t")
                        fieldnames = self.remap(next(spamreader))
                        for row in spamreader:
                            pdec.append(dict(zip(fieldnames, row)))
                    all_decoy = all_decoy + pdec
            else:
                pdec = []
                with open(self.decoy_file, "r") as psm_decoy_file:
                    self.logger.info(self.decoy_file)
                    spamreader = csv.reader(psm_decoy_file, delimiter="\t")
                    fieldnames = self.remap(next(spamreader))
                    for row in spamreader:
                        pdec.append(dict(zip(fieldnames, row)))
                all_decoy = pdec

            # Combine the lists
            all_psms = all_target + all_decoy

        elif self.combined_files:
            if isinstance(self.combined_files, (list,)):
                all = []
                for c_files in self.combined_files:
                    c_all = []
                    with open(c_files, "r") as psm_file:
                        self.logger.info(c_files)
                        spamreader = csv.reader(psm_file, delimiter="\t")
                        fieldnames = self.remap(next(spamreader))
                        for row in spamreader:
                            c_all.append(dict(zip(fieldnames, row)))
                    all = all + c_all
            else:
                c_all = []
                with open(self.combined_files, "r") as psm_file:
                    self.logger.info(self.combined_files)
                    spamreader = csv.reader(psm_file, delimiter="\t")
                    fieldnames = self.remap(next(spamreader))
                    for row in spamreader:
                        c_all.append(dict(zip(fieldnames, row)))
                all = c_all
            all_psms = all

        elif self.directory:
            all_files = os.listdir(self.directory)
            all = []
            for files in all_files:
                p = []
                with open(files, "r") as psm_file:
                    self.logger.info(self.target_file)
                    spamreader = csv.reader(psm_file, delimiter="\t")
                    fieldnames = self.remap(next(spamreader))
                    for row in spamreader:
                        p.append(dict(zip(fieldnames, row)))
                all = all + p
            all_psms = all

        psms_all_filtered = []
        for psms in all_psms:
            if self.POSTERIOR_ERROR_PROB in psms.keys():
                try:
                    float(psms[self.POSTERIOR_ERROR_PROB])
                    psms_all_filtered.append(psms)
                except ValueError as e:
                    pass
            else:
                try:
                    float(psms[self.scoring_variable])
                    psms_all_filtered.append(psms)
                except ValueError as e:
                    pass

        # Filter by pep
        try:
            self.logger.info("Sorting by {}".format(self.POSTERIOR_ERROR_PROB))
            all_psms = sorted(
                psms_all_filtered,
                key=lambda x: float(x[self.POSTERIOR_ERROR_PROB]),
                reverse=False,
            )
        except KeyError:
            self.logger.info(
                "Cannot Sort by {} the values do not exist".format(
                    self.POSTERIOR_ERROR_PROB
                )
            )
            self.logger.info("Sorting by {}".format(self.scoring_variable))
            if self.parameter_file_object.score_type == "additive":
                all_psms = sorted(
                    psms_all_filtered,
                    key=lambda x: float(x[self.scoring_variable]),
                    reverse=True,
                )
            if self.parameter_file_object.score_type == "multiplicative":
                all_psms = sorted(
                    psms_all_filtered,
                    key=lambda x: float(x[self.scoring_variable]),
                    reverse=False,
                )

        # TODO
        # TRY TO GET PERC_ALL AS A GENERATOR
        # Can do this... just give the option to feed a combined file... and loop over all the files...
        # Hmm... problem is it still needs to be sorted by perc score....

        list_of_psm_objects = []
        peptide_tracker = set()
        all_sp_proteins = set(self.digest_class.swiss_prot_protein_set)
        # We only want to get unique peptides... using all messes up scoring...
        # Create Psm objects with the identifier, percscore, qvalue, pepvalue, and possible proteins...

        peptide_to_protein_dictionary = self.digest_class.peptide_to_protein_dictionary

        # TODO
        # make this for loop a generator...
        self.logger.info("Length of PSM Data: {}".format(len(all_psms)))
        for psm_info in all_psms:
            current_peptide = psm_info[self.PEPTIDE]
            # Define the Psm...
            if current_peptide not in peptide_tracker:
                p = Psm(identifier=current_peptide)
                # Attempt to add variables from PSM info...
                # If they do not exist in the psm info then we skip...
                try:
                    p.percscore = float(psm_info[self.SCORE])
                except KeyError:
                    pass
                try:
                    p.qvalue = float(psm_info[self.Q_VALUE])
                except KeyError:
                    pass
                try:
                    p.pepvalue = float(psm_info[self.POSTERIOR_ERROR_PROB])
                except KeyError:
                    pass
                # If user has a custom score IE not q-value or pep_value...
                if self.load_custom_score:
                    # Then we look for it...
                    p.custom_score = float(psm_info[self.scoring_variable])
                p.possible_proteins = []
                p.possible_proteins.append(psm_info[self.PROTEIN_IDS])
                for (
                    alternative_protein_keys
                ) in self.MAX_ALTERNATIVE_PROTEIN_COLUMN_NAMES:
                    try:
                        if psm_info[alternative_protein_keys]:
                            p.possible_proteins.append(
                                psm_info[alternative_protein_keys]
                            )
                    except KeyError:
                        break
                # Remove potential Repeats
                if self.parameter_file_object.inference_type != "none":
                    p.possible_proteins = list(set(p.possible_proteins))

                # Get PSM ID
                p.psm_id = psm_info[self.PSMID]

                # Split peptide if flanking
                current_peptide = Psm.split_peptide(peptide_string=current_peptide)

                if not current_peptide.isupper() or not current_peptide.isalpha():
                    # If we have mods remove them...
                    peptide_string = current_peptide.upper()
                    stripped_peptide = Psm.remove_peptide_mods(peptide_string)
                    current_peptide = stripped_peptide
                # Add the other possible_proteins from insilicodigest here...

                try:
                    current_alt_proteins = list(
                        peptide_to_protein_dictionary[current_peptide]
                    )  # This peptide needs to be scrubbed of Mods...
                except KeyError:
                    current_alt_proteins = []
                    self.logger.warning(
                        "Peptide {} was not found in the supplied DB for Proteins {}".format(
                            current_peptide, ";".join(p.possible_proteins)
                        )
                    )
                    for poss_prot in p.possible_proteins:
                        self.digest_class.peptide_to_protein_dictionary.setdefault(
                            current_peptide, set()
                        ).add(poss_prot)
                        self.digest_class.protein_to_peptide_dictionary.setdefault(
                            poss_prot, set()
                        ).add(current_peptide)
                        self.logger.info(
                            "Adding Peptide {} and Protein {} to Digest dictionaries".format(
                                current_peptide, poss_prot
                            )
                        )

                # Sort Alt Proteins by Swissprot then Trembl...
                identifiers_sorted = DataStore.sort_protein_strings(
                    protein_string_list=current_alt_proteins,
                    sp_proteins=all_sp_proteins,
                    decoy_symbol=self.parameter_file_object.decoy_symbol,
                )

                # Restrict to 50 possible proteins
                p = self._fix_alternative_proteins(
                    append_alt_from_db=self.append_alt_from_db,
                    identifiers_sorted=identifiers_sorted,
                    max_proteins=self.MAX_ALLOWED_ALTERNATIVE_PROTEINS,
                    psm=p,
                    parameter_file_object=self.parameter_file_object,
                )

                list_of_psm_objects.append(p)
                peptide_tracker.add(current_peptide)

        self.psms = list_of_psm_objects

        self.logger.info("Length of PSM Data: {}".format(len(self.psms)))

        self.logger.info("Finished GenericReader.read_psms...")
