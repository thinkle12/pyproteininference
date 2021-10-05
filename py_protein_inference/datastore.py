import collections
import logging
import sys

from Bio import SeqIO

from py_protein_inference.inference import Inference
from py_protein_inference.physical import Protein, Psm
from py_protein_inference.scoring import Score

logger = logging.getLogger(__name__)

# set up our logger
logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


class DataStore(object):
    """
    The following Class serves as the data storage object for a protein inference analysis
    The class serves as a central point that is accessed at virtually every PI processing step


    Attributes:
        main_data_form (list): List of unrestricted Psm objects
        parameter_file_object (py_protein_inference.parameters.ProteinInferenceParameter): protein inference parameter
            object
        restricted_peptides (list): List of non flaking peptide strings present in the current analysis
        main_data_restricted (list): List of restricted :py:class:`py_protein_inference.physical.Psm` objects.
            Restriction is based on the parameter_file_object and the object is created by function
                :py:meth:`py_protein_inference.datastore.DataStore.restrict_psm_data`
        scored_proteins (list): List of scored :py:class:`py_protein_inference.physical.Protein` objects.
            Output from scoring methods from :py:mod:`py_protein_inference.scoring`
        grouped_scored_proteins (list): List of scored :py:class:`py_protein_inference.physical.Protein`
            objects that have been grouped and sorted. Output from
                :py:meth:`py_protein_inference.inference.Inference.run_inference` method
        scoring_input (list): List of non-scored :py:class:`py_protein_inference.physical.Protein` objects.
            Output from :py:meth:`py_protein_inference.datastore.DataStore.create_scoring_input`
        picked_proteins_scored (list): List of :py:class:`py_protein_inference.physical.Protein` objects that pass
            the protein picker algorithm (:py:meth:`py_protein_inference.datastore.DataStore.protein_picker`)
        picked_proteins_removed (list): List of :py:class:`py_protein_inference.physical.Protein` objects that do not
            pass the protein picker algorithm (:py:meth:`py_protein_inference.datastore.DataStore.protein_picker`)
        protein_peptide_dictionary (collections.defaultdict): Dictionary of protein strings (keys) that map to sets
            of peptide strings based on the peptides and proteins found in the search. Protein -> set(Peptides)
        peptide_protein_dictionary (collections.defaultdict): Dictionary of peptide strings (keys) that map to sets
            of protein strings based on the peptides and proteins found in the search. Peptide -> set(Proteins)
        high_low_better (str): Variable that indicates whether a higher or a lower protein score is better.
            This is necessary to sort Protein objects by score properly. Can either be "higher" or "lower"
        psm_score (str): Variable that indicates the :py:class:`py_protein_inference.physical.Psm`
            score being used in the analysis to generate :py:class:`py_protein_inference.physical.Protein` scores
        protein_score (str): String to indicate the protein score method used
        short_protein_score (str): Short String to indicate the protein score method used
        protein_group_objects (list): List of scored :py:class:`py_protein_inference.physical.ProteinGroup`
            objects that have been grouped and sorted. Output from
             :py:meth:`py_protein_inference.inference.Inference.run_inference` method
        decoy_symbol (str): String that is used to differentiate between decoy proteins and target proteins. Ex: "##"
        digest (py_protein_inference.in_silico_digest.Digest): Digest object
            :py:class:`py_protein_inference.in_silico_digest.Digest`
        SCORE_MAPPER (dict): Dictionary that maps potential scores in input files to internal score names
        CUSTOM_SCORE_KEY (str): String that indicates a custom score is being used

    """

    SCORE_MAPPER = {
        "q_value": "qvalue",
        "pep_value": "pepvalue",
        "perc_score": "percscore",
        "score": "percscore",
        "q-value": "qvalue",
        "posterior_error_prob": "pepvalue",
        "posterior_error_probability": "pepvalue",
    }

    CUSTOM_SCORE_KEY = "custom_score"

    HIGHER_PSM_SCORE = "higher"
    LOWER_PSM_SCORE = "lower"

    def __init__(self, reader, digest, validate=True):
        """

        Args:
            reader (py_protein_inference.reader.Reader): Reader object :py:class:`protein_infernece.reader.Reader`
            digest (py_protein_inference.in_silico_digest.Digest): Digest object
                :py:class:`protein_infernece.in_silico_digest.Digest`
            validate (bool): True/False to indicate if the input data should be validated

        Example:
            >>> py_protein_inference.datastore.DataStore(reader = reader, digest=digest)


        """
        # If the reader class is from a percolator.psms then define main_data_form as reader.psms
        # main_data_form is the starting point for all other analyses
        self._init_validate(reader=reader)

        self.parameter_file_object = reader.parameter_file_object  # Parameter object
        self.main_data_restricted = None  # PSM data post restriction
        self.scored_proteins = []  # List of scored Protein objects
        self.grouped_scored_proteins = []  # List of sorted scored Protein objects
        self.scoring_input = None  # List of non scored Protein objects
        self.picked_proteins_scored = None  # List of Protein objects after picker algorithm
        self.picked_proteins_removed = None  # Protein objects removed via picker
        self.protein_peptide_dictionary = None
        self.peptide_protein_dictionary = None
        self.high_low_better = None  # Variable that indicates whether a higher or lower protein score is better
        self.psm_score = None  # PSM Score used
        self.protein_score = None
        self.short_protein_score = None
        self.protein_group_objects = []  # List of sorted protein group objects
        self.decoy_symbol = self.parameter_file_object.decoy_symbol  # Decoy symbol from parameter file
        self.digest = digest  # Digest object

        # Run Checks and Validations
        if validate:
            self.validate_psm_data()
            self.validate_digest()
            self.check_data_consistency()

        # Run method to fix our parameter object if necessary
        self.parameter_file_object.fix_parameters_from_datastore(data=self)

    def get_sorted_identifiers(self, scored=True):
        """
        Retrieves a sorted list of protein strings present in the analysis

        Args:
            scored (bool): True/False to indicate if we should return scored or non-scored identifiers

        Returns:
            list: List of sorted protein identifier strings

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> sorted_proteins = data.get_sorted_identifiers(scored=True)
        """

        if scored:
            self._validate_scored_proteins()
            if self.picked_proteins_scored:
                proteins = set([x.identifier for x in self.picked_proteins_scored])
            else:
                proteins = set([x.identifier for x in self.scored_proteins])
        else:
            self._validate_scoring_input()
            proteins = [x.identifier for x in self.scoring_input]

        all_sp_proteins = set(self.digest.swiss_prot_protein_set)

        our_target_sp_proteins = sorted([x for x in proteins if x in all_sp_proteins and self.decoy_symbol not in x])
        our_decoy_sp_proteins = sorted([x for x in proteins if x in all_sp_proteins and self.decoy_symbol in x])

        our_target_tr_proteins = sorted(
            [x for x in proteins if x not in all_sp_proteins and self.decoy_symbol not in x]
        )
        our_decoy_tr_proteins = sorted([x for x in proteins if x not in all_sp_proteins and self.decoy_symbol in x])

        our_proteins_sorted = (
            our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins
        )

        return our_proteins_sorted

    @classmethod
    def sort_protein_group_objects(cls, protein_group_objects, higher_or_lower):
        """
        Class Method to sort a list of :py:class:`protein_inferenece.physical.ProteinGroup` objects by
        score and number of peptides

        Args:
            protein_group_objects (list): list of :py:class:`protein_inferenece.physical.ProteinGroup` objects
            higher_or_lower (str): String to indicate if a "higher" or "lower" protein score is "better"

        Returns:
            list: list of sorted :py:class:`protein_inferenece.physical.ProteinGroup` objects

        Example:
            >>> list_of_group_objects = py_protein_inference.datastore.DataStore.sort_protein_group_objects(
            >>>     protein_group_objects=list_of_group_objects, higher_or_lower="higher"
            >>> )
        """
        if higher_or_lower == cls.LOWER_PSM_SCORE:

            protein_group_objects = sorted(
                protein_group_objects,
                key=lambda k: (
                    k.proteins[0].score,
                    -k.proteins[0].num_peptides,
                ),
                reverse=False,
            )
        elif higher_or_lower == cls.HIGHER_PSM_SCORE:

            protein_group_objects = sorted(
                protein_group_objects,
                key=lambda k: (
                    k.proteins[0].score,
                    k.proteins[0].num_peptides,
                ),
                reverse=True,
            )

        return protein_group_objects

    @classmethod
    def sort_protein_objects(cls, grouped_protein_objects, higher_or_lower):
        """
        Class Method to sort a list of :py:class:`protein_inferenece.physical.Protein` objects by score and number of
        peptides

        Args:
            grouped_protein_objects (list): list of :py:class:`protein_inferenece.physical.Protein` objects
            higher_or_lower (str): String to indicate if a "higher" or "lower" protein score is "better"

        Returns:
            list: list of sorted :py:class:`protein_inferenece.physical.Protein` objects

        Example:
            >>> scores_grouped = py_protein_inference.datastore.DataStore.sort_protein_objects(
            >>>     grouped_protein_objects=scores_grouped, higher_or_lower="higher"
            >>> )
        """
        if higher_or_lower == cls.LOWER_PSM_SCORE:
            grouped_protein_objects = sorted(
                grouped_protein_objects,
                key=lambda k: (k[0].score, -k[0].num_peptides),
                reverse=False,
            )
        if higher_or_lower == cls.HIGHER_PSM_SCORE:
            grouped_protein_objects = sorted(
                grouped_protein_objects,
                key=lambda k: (k[0].score, k[0].num_peptides),
                reverse=True,
            )
        return grouped_protein_objects

    @classmethod
    def sort_protein_sub_groups(cls, protein_list, higher_or_lower):
        """
        Method to sort protein sub lists

        Args:
            protein_list (list): List of :py:class:`protein_inferenece.physical.Protein` objects to be sorted
            higher_or_lower (str): String to indicate if a "higher" or "lower" protein score is "better"

        Returns:
            list: List of :py:class:`protein_inferenece.physical.Protein` objects to be sorted by score and number of
            peptides

        """

        # Sort the groups based on higher or lower indication, secondarily sort the groups based on number of unique
        # peptides
        # We use the index [1:] as we do not wish to sort the lead protein...
        if higher_or_lower == cls.LOWER_PSM_SCORE:
            protein_list[1:] = sorted(
                protein_list[1:],
                key=lambda k: (float(k.score), -float(k.num_peptides)),
                reverse=False,
            )
        if higher_or_lower == cls.HIGHER_PSM_SCORE:
            protein_list[1:] = sorted(
                protein_list[1:],
                key=lambda k: (float(k.score), float(k.num_peptides)),
                reverse=True,
            )

        return protein_list

    def get_psm_data(self):
        """
        Method to retrieve a list of :py:class:`py_protein_inference.physical.Psm` objects.
        Retrieves restricted data if the data has been restricted or all of the data if the data has not been restricted

        Returns:
            list: list of :py:class:`py_protein_inference.physical.Psm` objects

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> psm_data = data.get_psm_data()
        """
        if not self.main_data_restricted and not self.main_data_form:
            raise ValueError(
                "Both main_data_restricted and main_data_form variables are empty. Please re-load the DataStore "
                "object with a properly loaded Reader object."
            )

        if self.main_data_restricted:
            psm_data = self.main_data_restricted
        else:
            psm_data = self.main_data_form

        return psm_data

    def get_protein_data(self):
        """
        Method to retrieve a list of :py:class:`py_protein_inference.physical.Protein` objects.
        Retrieves picked and scored data if the data has been picked and scored or just the scored data if the data has
         not been picked.

        Returns:
            list: list of :py:class:`py_protein_inference.physical.Protein` objects

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> # Data must ben ran through a py_protein_inference.scoring.Score method
            >>> protein_data = data.get_protein_data()
        """

        if self.picked_proteins_scored:
            scored_proteins = self.picked_proteins_scored
        else:
            scored_proteins = self.scored_proteins

        return scored_proteins

    def get_protein_identifiers_from_psm_data(self):
        """
        Method to retrieve a list of lists of all possible protein identifiers from the psm data

        Returns:
            list: list of lists of protein strings

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> protein_strings = data.get_protein_identifiers_from_psm_data()
        """
        psm_data = self.get_psm_data()

        proteins = [x.possible_proteins for x in psm_data]

        return proteins

    def get_q_values(self):
        """
        Method to retrieve a list of all q values for all PSMs

        Returns:
            list: list of floats (q values)

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> q = data.get_q_values()
        """
        psm_data = self.get_psm_data()

        q_values = [x.qvalue for x in psm_data]

        return q_values

    def get_pep_values(self):
        """
        Method to retrieve a list of all posterior error probabilities for all PSMs

        Returns:
            list: list of floats (pep values)

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> pep = data.get_pep_values()
        """
        psm_data = self.get_psm_data()

        pep_values = [x.pepvalue for x in psm_data]

        return pep_values

    def get_protein_information_dictionary(self):
        """
        Method to retrieve a dictionary of scores for each peptide

        Returns:
            dict: dictionary of scores for each protein

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> protein_dict = data.get_protein_information_dictionary()
        """
        psm_data = self.get_psm_data()

        protein_psm_score_dictionary = collections.defaultdict(list)

        # Loop through all Psms
        for psms in psm_data:
            # Loop through all proteins
            for prots in psms.possible_proteins:
                protein_psm_score_dictionary[prots].append(
                    {
                        "peptide": psms.identifier,
                        "Qvalue": psms.qvalue,
                        "PosteriorErrorProbability": psms.pepvalue,
                        "Percscore": psms.percscore,
                    }
                )

        return protein_psm_score_dictionary

    def restrict_psm_data(self, remove1pep=True):
        """
        Method to restrict the input of PSM data (:py:class:`py_protein_inference.physical.Psm`) objects.
        This method is central to the py_protein_inference module and is able to restrict the Psm data by:
        Q value, Pep Value, Percolator Score, Peptide Length, and Custom Score Input.
        Restriction values are pulled from the :py:class:`py_protein_inference.parameters.ProteinInferenceParameter`
        object

        This method sets the :attr:`main_data_restricted` and :attr:`restricted_peptides` Attributes for the
         DataStore object

        Args:
            remove1pep (bool): True/False on whether or not to remove PEP values that equal 1 even if other restrictions
                are set to not restrict.

        Returns:
            None

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> data.restrict_psm_data(remove1pep=True)
        """

        # Validate that we have the main data variable
        self._validate_main_data_form()

        logger.info("Restricting PSM data")

        peptide_length = self.parameter_file_object.restrict_peptide_length
        posterior_error_prob_threshold = self.parameter_file_object.restrict_pep
        q_value_threshold = self.parameter_file_object.restrict_q
        custom_threshold = self.parameter_file_object.restrict_custom

        main_psm_data = self.main_data_form
        logger.info("Length of main data: {}".format(len(self.main_data_form)))
        # If restrict_main_data is called, we automatically discard everything that has a PEP of 1
        if remove1pep and posterior_error_prob_threshold:
            main_psm_data = [x for x in main_psm_data if x.pepvalue != 1]

        # Restrict peptide length and posterior error probability
        if peptide_length and posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.stripped_peptide) >= peptide_length and psms.pepvalue < float(
                    posterior_error_prob_threshold
                ):
                    restricted_data.append(psms)

        # Restrict peptide length only
        if peptide_length and not posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.stripped_peptide) >= peptide_length:
                    restricted_data.append(psms)

        # Restrict peptide length, posterior error probability, and qvalue
        if peptide_length and posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if (
                    len(psms.stripped_peptide) >= peptide_length
                    and psms.pepvalue < float(posterior_error_prob_threshold)
                    and psms.qvalue < float(q_value_threshold)
                ):
                    restricted_data.append(psms)

        # Restrict peptide length and qvalue
        if peptide_length and not posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if len(psms.stripped_peptide) >= peptide_length and psms.qvalue < float(q_value_threshold):
                    restricted_data.append(psms)

        # Restrict posterior error probability and q value
        if not peptide_length and posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if psms.pepvalue < float(posterior_error_prob_threshold) and psms.qvalue < float(q_value_threshold):
                    restricted_data.append(psms)

        # Restrict qvalue only
        if not peptide_length and not posterior_error_prob_threshold and q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if psms.qvalue < float(q_value_threshold):
                    restricted_data.append(psms)

        # Restrict posterior error probability only
        if not peptide_length and posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = []
            for psms in main_psm_data:
                if psms.pepvalue < float(posterior_error_prob_threshold):
                    restricted_data.append(psms)

        # Restrict nothing... (only PEP gets restricted - takes everything less than 1)
        if not peptide_length and not posterior_error_prob_threshold and not q_value_threshold:
            restricted_data = main_psm_data

        if custom_threshold:
            custom_restricted = []
            if self.parameter_file_object.psm_score_type == Score.MULTIPLICATIVE_SCORE_TYPE:
                for psms in restricted_data:
                    if psms.custom_score <= custom_threshold:
                        custom_restricted.append(psms)

            if self.parameter_file_object.psm_score_type == Score.ADDITIVE_SCORE_TYPE:
                for psms in restricted_data:
                    if psms.custom_score >= custom_threshold:
                        custom_restricted.append(psms)

            restricted_data = custom_restricted

        self.main_data_restricted = restricted_data

        logger.info("Length of restricted data: {}".format(len(restricted_data)))

        self.restricted_peptides = [x.non_flanking_peptide for x in restricted_data]

    def create_scoring_input(self):
        """
        Method to create the scoring input.
        This method initializes a list of :py:class:`py_protein_inference.physical.Protein` objects to get them ready
        to be scored by :py:mod:`py_protein_inference.scoring.Score` methods
        This method also takes into account the inference type and aggregates peptides -> proteins accordingly.

        This method sets the :attr:`scoring_input` and :attr:`score` Attributes for the DataStore object

        The score selected comes from the protein inference parameter object

        Returns:
            None

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> data.create_scoring_input()
        """

        logger.info("Creating Scoring Input")

        psm_data = self.get_psm_data()

        protein_psm_dict = collections.defaultdict(list)

        try:
            score_key = self.SCORE_MAPPER[self.parameter_file_object.psm_score]
        except KeyError:
            score_key = self.CUSTOM_SCORE_KEY

        if self.parameter_file_object.inference_type != Inference.PEPTIDE_CENTRIC:
            # Loop through all Psms
            for psms in psm_data:
                psms.assign_main_score(score=score_key)
                # Loop through all proteins
                for prots in psms.possible_proteins:
                    protein_psm_dict[prots].append(psms)

        else:
            self.peptide_to_protein_dictionary()
            sp_proteins = self.digest.swiss_prot_protein_set
            for psms in psm_data:

                # Assign main score
                psms.assign_main_score(score=score_key)
                protein_set = self.peptide_protein_dictionary[psms.non_flanking_peptide]
                # Sort protein_set by sp-alpha, decoy-sp-alpha, tr-alpha, decoy-tr-alpha
                sorted_protein_list = self.sort_protein_strings(
                    protein_string_list=protein_set,
                    sp_proteins=sp_proteins,
                    decoy_symbol=self.parameter_file_object.decoy_symbol,
                )
                # Restrict the number of identifiers by the value in param file max_identifiers_peptide_centric
                sorted_protein_list = sorted_protein_list[: self.parameter_file_object.max_identifiers_peptide_centric]
                protein_name = ";".join(sorted_protein_list)
                protein_psm_dict[protein_name].append(psms)

        protein_list = []
        for pkey in sorted(protein_psm_dict.keys()):
            protein_object = Protein(identifier=pkey)
            protein_object.psms = protein_psm_dict[pkey]
            protein_object.raw_peptides = set([x.identifier for x in protein_psm_dict[pkey]])
            protein_list.append(protein_object)

        self.psm_score = self.parameter_file_object.psm_score
        self.scoring_input = protein_list

    def protein_to_peptide_dictionary(self):
        """
        Method that returns a map of protein strings to sets of peptide strings and is essentially half
         of a BiPartite graph
        This method sets the :attr:`protein_peptide_dictionary` Attribute for the
        :py:class:`py_protein_inference.datastore.DataStore` object

        Returns:
            collections.defaultdict: Dictionary of protein strings (keys) that map to sets of peptide strings based
            on the peptides and proteins found in the search. Protein -> set(Peptides)

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> protein_peptide_dict = data.protein_to_peptide_dictionary()
        """
        psm_data = self.get_psm_data()

        res_pep_set = set(self.restricted_peptides)
        default_dict_proteins = collections.defaultdict(set)
        for peptide_objects in psm_data:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = peptide_objects.non_flanking_peptide
                if cur_peptide in res_pep_set:
                    default_dict_proteins[prots].add(cur_peptide)

        self.protein_peptide_dictionary = default_dict_proteins

        return default_dict_proteins

    def peptide_to_protein_dictionary(self):
        """
        Method that returns a map of peptide strings to sets of protein strings and is essentially half of a
        BiPartite graph
        This method sets the :attr:`peptide_protein_dictionary` Attribute for the
         :py:class:`py_protein_inference.datastore.DataStore` object

        Returns:
            collections.defaultdict: Dictionary of peptide strings (keys) that map to sets of protein strings based
                on the peptides and proteins found in the search. Peptide -> set(Proteins)

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> peptide_protein_dict = data.peptide_to_protein_dictionary()
        """
        psm_data = self.get_psm_data()

        res_pep_set = set(self.restricted_peptides)
        default_dict_peptides = collections.defaultdict(set)
        for peptide_objects in psm_data:
            for prots in peptide_objects.possible_proteins:
                cur_peptide = peptide_objects.non_flanking_peptide
                if cur_peptide in res_pep_set:
                    default_dict_peptides[cur_peptide].add(prots)
                else:
                    pass

        self.peptide_protein_dictionary = default_dict_peptides

        return default_dict_peptides

    def unique_to_leads_peptides(self):
        """
        Method to retrieve peptides that are unique based on the data from the searches
        (Not based on the database digestion)

        Returns:
            set

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> unique_peps = data.unique_to_leads_peptides()
        """
        if self.grouped_scored_proteins:
            lead_peptides = [list(x[0].peptides) for x in self.grouped_scored_proteins]
            flat_peptides = [item for sublist in lead_peptides for item in sublist]
            counted_peps = collections.Counter(flat_peptides)
            unique_to_leads_peptides = set([x for x in counted_peps if counted_peps[x] == 1])
        else:
            unique_to_leads_peptides = set()

        return unique_to_leads_peptides

    def higher_or_lower(self):
        """
        Method to determine if a higher or lower score is better for a given combination of score input and score type

        This method sets the :attr:`high_low_better` Attribute for the DataStore object.

        This method depends on the output from the Score class to be sorted properly from best to worst score

        Returns:
            str: String indicating "higher" or "lower" depending on if a higher or lower score is a better protein score

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> high_low = data.higher_or_lower()
        """

        if not self.high_low_better:
            logger.info("Determining If a higher or lower score is better based on scored proteins")
            worst_score = self.scored_proteins[-1].score
            best_score = self.scored_proteins[0].score

            if float(best_score) > float(worst_score):
                higher_or_lower = self.HIGHER_PSM_SCORE

            if float(best_score) < float(worst_score):
                higher_or_lower = self.LOWER_PSM_SCORE

            logger.info("best score = {}".format(best_score))
            logger.info("worst score = {}".format(worst_score))

            if best_score == worst_score:
                raise ValueError(
                    "Best and Worst scores were identical, equal to {}. Score type {} produced the error, "
                    "please change psm_score type.".format(best_score, self.psm_score)
                )

            self.high_low_better = higher_or_lower

        else:
            higher_or_lower = self.high_low_better

        return higher_or_lower

    def get_protein_identifiers(self, data_form):
        """
        Method to retrieve the protein string identifiers

        Args:
            data_form (str): Can be one of the following: "main", "restricted", "picked", "picked_removed"

        Returns:
            list: list of protein identifier strings

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> protein_strings = data.get_protein_identifiers(data_form="main")
        """
        if data_form == "main":
            # All the data (unrestricted)
            data_to_select = self.main_data_form
            prots = [[x.possible_proteins] for x in data_to_select]
            proteins = prots

        if data_form == "restricted":
            # Proteins that pass certain restriction criteria (peptide length, pep, qvalue)
            data_to_select = self.main_data_restricted
            prots = [[x.possible_proteins] for x in data_to_select]
            proteins = prots

        if data_form == "picked":
            # Here we look at proteins that are 'picked' (aka the proteins that beat out their matching target/decoy)
            data_to_select = self.picked_proteins_scored
            prots = [x.identifier for x in data_to_select]
            proteins = prots

        if data_form == "picked_removed":
            # Here we look at the proteins that were removed due to picking (aka the proteins that
            # have a worse score than their target/decoy counterpart)
            data_to_select = self.picked_proteins_removed
            prots = [x.identifier for x in data_to_select]
            proteins = prots

        return proteins

    def get_protein_information(self, protein_string):
        """
        Method to retrieve attributes for a specific scored protein

        Args:
            protein_string (str): Protein Identifier String

        Returns:
            list: list of protein attributes

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> protein_attr = data.get_protein_information(protein_string="RAF1_HUMAN|P04049")
        """
        all_scored_protein_data = self.scored_proteins
        identifiers = [x.identifier for x in all_scored_protein_data]
        protein_scores = [x.score for x in all_scored_protein_data]
        groups = [x.group_identification for x in all_scored_protein_data]
        reviewed = [x.reviewed for x in all_scored_protein_data]
        peptides = [x.peptides for x in all_scored_protein_data]
        # Peptide scores currently broken...
        peptide_scores = [x.peptide_scores for x in all_scored_protein_data]
        picked = [x.picked for x in all_scored_protein_data]
        num_peptides = [x.num_peptides for x in all_scored_protein_data]

        main_index = identifiers.index(protein_string)

        list_structure = [
            [
                "identifier",
                "protein_score",
                "groups",
                "reviewed",
                "peptides",
                "peptide_scores",
                "picked",
                "num_peptides",
            ]
        ]
        list_structure.append([protein_string])
        list_structure[-1].append(protein_scores[main_index])
        list_structure[-1].append(groups[main_index])
        list_structure[-1].append(reviewed[main_index])
        list_structure[-1].append(peptides[main_index])
        list_structure[-1].append(peptide_scores[main_index])
        list_structure[-1].append(picked[main_index])
        list_structure[-1].append(num_peptides[main_index])

        return list_structure

    def exclude_non_distinguishing_peptides(self, protein_subset_type="hard"):
        """
        Method to Exclude peptides that are not distinguishing on either the search or database level

        The method sets the :attr:`scoring_input` and :attr:`restricted_peptides` variables for the
        :py:class:`py_protein_inference.datastore.DataStore` object

        Args:
            protein_subset_type (str): Either "hard" or "soft". Hard will select distinguishing peptides based on
                the database digestion. "soft" will only use peptides identified in the search.

        Returns:
            None

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> data.exclude_non_distinguishing_peptides(protein_subset_type="hard")
        """

        logger.info("Applying Exclusion Model")

        our_proteins_sorted = self.get_sorted_identifiers(scored=False)

        if protein_subset_type == "hard":
            # Hard protein subsetting defines protein subsets on the digest level (Entire protein is used)
            # This is how Percolator PI does subsetting
            peptides = [self.digest.protein_to_peptide_dictionary[x] for x in our_proteins_sorted]
        elif protein_subset_type == "soft":
            # Soft protein subsetting defines protein subsets on the Peptides identified from the search
            peptides = [set(x.raw_peptides) for x in self.scoring_input]
        else:
            # If neither is dfined we do "hard" exclusion
            peptides = [self.digest.protein_to_peptide_dictionary[x] for x in our_proteins_sorted]

        # Get frozen set of peptides....
        # We will also have a corresponding list of proteins...
        # They will have the same index...
        peptide_sets = [frozenset(e) for e in peptides]
        # Find a way to sort this list of sets...
        # We can sort the sets if we sort proteins from above...
        logger.info("{} number of peptide sets".format(len(peptide_sets)))
        non_subset_peptide_sets = set()
        i = 0
        # Get all peptide sets that are not a subset...
        while peptide_sets:
            i = i + 1
            peptide_set = peptide_sets.pop()
            if any(peptide_set.issubset(s) for s in peptide_sets) or any(
                peptide_set.issubset(s) for s in non_subset_peptide_sets
            ):
                continue
            else:
                non_subset_peptide_sets.add(peptide_set)
            if i % 10000 == 0:
                logger.info("Parsed {} Peptide Sets".format(i))

        logger.info("Parsed {} Peptide Sets".format(i))

        # Get their index from peptides which is the initial list of sets...
        list_of_indeces = []
        for pep_sets in non_subset_peptide_sets:
            ind = peptides.index(pep_sets)
            list_of_indeces.append(ind)

        non_subset_proteins = set([our_proteins_sorted[x] for x in list_of_indeces])

        logger.info("Removing direct subset Proteins from the data")
        # Remove all proteins from scoring input that are a subset of another protein...
        self.scoring_input = [x for x in self.scoring_input if x.identifier in non_subset_proteins]

        logger.info("{} proteins in scoring input after removing subset proteins".format(len(self.scoring_input)))

        # For all the proteins that are not a complete subset of another protein...
        # Get the raw peptides...
        raw_peps = [x.raw_peptides for x in self.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        flat_peptides = [Psm.split_peptide(peptide_string=item) for sublist in raw_peps for item in sublist]

        # Count the number of peptides in this list...
        # This is the number of proteins this peptide maps to....
        counted_peptides = collections.Counter(flat_peptides)

        # If the count is greater than 1... exclude the protein entirely from scoring input... :)
        raw_peps_good = set([x for x in counted_peptides.keys() if counted_peptides[x] <= 1])

        # Alter self.scoring_input by removing psms and peptides that are not in raw_peps_good
        current_score_input = list(self.scoring_input)
        for j in range(len(current_score_input)):
            k = j + 1
            psm_list = []
            new_raw_peptides = []
            current_psms = current_score_input[j].psms
            current_raw_peptides = current_score_input[j].raw_peptides

            for psm_scores in current_psms:
                if psm_scores.non_flanking_peptide in raw_peps_good:
                    psm_list.append(psm_scores)

            for rp in current_raw_peptides:
                if Psm.split_peptide(peptide_string=rp) in raw_peps_good:
                    new_raw_peptides.append(rp)

            current_score_input[j].psms = psm_list
            current_score_input[j].raw_peptides = new_raw_peptides

            if k % 10000 == 0:
                logger.info("Redefined {} Peptide Sets".format(k))

        logger.info("Redefined {} Peptide Sets".format(j))

        filtered_score_input = [x for x in current_score_input if x.psms]

        self.scoring_input = filtered_score_input

        # Recompute the flat peptides
        raw_peps = [x.raw_peptides for x in self.scoring_input if x.identifier in non_subset_proteins]

        # Make the raw peptides a flat list
        new_flat_peptides = set([Psm.split_peptide(peptide_string=item) for sublist in raw_peps for item in sublist])

        self.scoring_input = [x for x in self.scoring_input if x.psms]

        self.restricted_peptides = [x for x in self.restricted_peptides if x in new_flat_peptides]

    def protein_picker(self):
        """
        Method to run the protein picker algorithm.

        Proteins must be scored first with :py:meth:`py_protein_inference.scoring.Score.score_psms`

        The algorithm will match target and decoy proteins identified from the PSMs from the search.
        If a target and matching decoy is found then target/decoy competition is performed.
        In the Target/Decoy pair the protein with the better score is kept and the one with the worse score is
        discarded from the analysis

        The method sets the :attr:`picked_proteins_scored` and :attr:`picked_proteins_removed` variables for
        the :py:class:`py_protein_inference.datastore.DataStore` object

        Returns:
            None

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> data.protein_picker()
        """

        self._validate_scored_proteins()

        logger.info("Running Protein Picker")

        # Use higher or lower class to determine if a higher protein score or lower protein score is better
        # based on the scoring method used
        higher_or_lower = self.higher_or_lower()
        # Here we determine if a lower or higher score is better
        # Since all input is ordered from best to worst we can do the following

        index_to_remove = []
        # data.scored_proteins is simply a list of Protein objects...
        # Create list of all decoy proteins
        decoy_proteins = [x.identifier for x in self.scored_proteins if self.decoy_symbol in x.identifier]
        # Create a list of all potential matching targets (some of these may not exist in the search)
        matching_targets = [x.replace(self.decoy_symbol, "") for x in decoy_proteins]

        # Create a list of all the proteins from the scored data
        all_proteins = [x.identifier for x in self.scored_proteins]
        logger.info("{} proteins scored".format(len(all_proteins)))

        total_targets = []
        total_decoys = []
        decoys_removed = []
        targets_removed = []
        # Loop over all decoys identified in the search
        logger.info("Picking Proteins...")
        for i in range(len(decoy_proteins)):
            cur_decoy_index = all_proteins.index(decoy_proteins[i])
            cur_decoy_protein_object = self.scored_proteins[cur_decoy_index]
            total_decoys.append(cur_decoy_protein_object.identifier)

            # Try, Except here because the matching target to the decoy may not be a result from the search
            try:
                cur_target_index = all_proteins.index(matching_targets[i])
                cur_target_protein_object = self.scored_proteins[cur_target_index]
                total_targets.append(cur_target_protein_object.identifier)

                if higher_or_lower == self.HIGHER_PSM_SCORE:
                    if cur_target_protein_object.score > cur_decoy_protein_object.score:
                        index_to_remove.append(cur_decoy_index)
                        decoys_removed.append(cur_decoy_index)
                        cur_target_protein_object.picked = True
                        cur_decoy_protein_object.picked = False
                    else:
                        index_to_remove.append(cur_target_index)
                        targets_removed.append(cur_target_index)
                        cur_decoy_protein_object.picked = True
                        cur_target_protein_object.picked = False

                if higher_or_lower == self.LOWER_PSM_SCORE:
                    if cur_target_protein_object.score < cur_decoy_protein_object.score:
                        index_to_remove.append(cur_decoy_index)
                        decoys_removed.append(cur_decoy_index)
                        cur_target_protein_object.picked = True
                        cur_decoy_protein_object.picked = False
                    else:
                        index_to_remove.append(cur_target_index)
                        targets_removed.append(cur_target_index)
                        cur_decoy_protein_object.picked = True
                        cur_target_protein_object.picked = False
            except ValueError:
                pass

        logger.info("{} total decoy proteins".format(len(total_decoys)))
        logger.info("{} matching target proteins also found in search".format(len(total_targets)))
        logger.info("{} decoy proteins to be removed".format(len(decoys_removed)))
        logger.info("{} target proteins to be removed".format(len(targets_removed)))

        logger.info("Removing Lower Scoring Proteins...")
        picked_list = []
        removed_proteins = []
        for protein_objects in self.scored_proteins:
            if protein_objects.picked:
                picked_list.append(protein_objects)
            else:
                removed_proteins.append(protein_objects)
        self.picked_proteins_scored = picked_list
        self.picked_proteins_removed = removed_proteins
        logger.info("Finished Removing Proteins")

    def calculate_q_values(self, regular=True):
        """
        Method calculates Q values FDR on the lead protein in the group on the :attr:`protein_group_objects`
        instance variable
        FDR is calculated As (2*decoys)/total if regular is set to True and is (decoys)/total if regular is set to False

        This method updates the :attr:`protein_group_objects` for the
        :py:class:`py_protein_inference.datastore.DataStore` object by updating the q_value variable of the
         :py:class:`py_protein_inference.physical.ProteinGroup` objects

        Returns:
            None

        Example:
            >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> # Data must be scored first
            >>> data.calculate_q_values()
        """

        self._validate_protein_group_objects()

        logger.info("Calculating Q values from the protein group objects")

        # pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x.proteins[0] for x in self.protein_group_objects]
        # Now pick out only the lead protein identifiers
        lead_proteins = [x.identifier for x in lead_score]

        lead_proteins.reverse()

        logger.info("Calculating FDRs")
        fdr_list = []
        for i in range(len(lead_proteins)):
            binary_decoy_target_list = [1 if elem.startswith(self.decoy_symbol) else 0 for elem in lead_proteins]
            total = len(lead_proteins)
            decoys = sum(binary_decoy_target_list)
            # Calculate FDR at every step starting with the entire list...
            # Delete first entry (worst score) every time we go through a cycle
            if regular:
                fdr = (2 * decoys) / (float(total))
            else:
                fdr = (decoys) / (float(total))
            fdr_list.append(fdr)
            del lead_proteins[0]

        qvalue_list = []
        new_fdr_list = []
        logger.info("Calculating Q Values")
        for fdrs in fdr_list:
            new_fdr_list.append(fdrs)
            qvalue = min(new_fdr_list)
            # qvalue = fdrs
            qvalue_list.append(qvalue)

        qvalue_list.reverse()

        logger.info("Assigning Q Values")
        for k in range(len(self.protein_group_objects)):
            self.protein_group_objects[k].q_value = qvalue_list[k]

        fdr_restricted = [x for x in self.protein_group_objects if x.q_value <= self.parameter_file_object.fdr]

        fdr_restricted_set = [self.grouped_scored_proteins[x] for x in range(len(fdr_restricted))]

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        logger.info(
            "Protein Group leads that pass with more than 1 PSM with a {} FDR = {}".format(
                self.parameter_file_object.fdr,
                str(len(fdr_restricted_set) - len(onehitwonders)),
            )
        )
        logger.info(
            "Protein Group lead One hit Wonders that pass {} FDR = {}".format(
                self.parameter_file_object.fdr, len(onehitwonders)
            )
        )

        logger.info(
            "Number of Protein groups that pass a {} percent FDR: {}".format(
                str(self.parameter_file_object.fdr * 100), len(fdr_restricted_set)
            )
        )

        logger.info("Finished Q value Calculation")

    def entrapment_fdr(self, true_database, false_discovery_rate=0.05):
        """
        Calculates Entrapment FDR on the lead protein in the groups.
        Input is a DataStore object, an entrapment database, and a false discovery rate

        Example:
             >>> data = py_protein_inference.datastore.DataStore(reader = reader, digest=digest)
            >>> # Data must be scored first
            >>> data.entrapment_fdr(true_database = "example_entrap.fasta",false_discovery_rate=0.1)

        FDR values are calculated As (entrapped proteins)/total

        This method is useful for calculating an entrapment FDR IF using a benchmark dataset with KNOWN protein content.
        Entrapment DB would be target proteins known to NOT be in the sample. However, these entrapped proteins
        should be in the main database that
        the search was searched against via comet/mascot
        """

        true_handle = SeqIO.parse(true_database, "fasta")
        true_proteins = []
        for records in true_handle:
            if not records.id.startsiwith(
                self.decoy_symbol
            ):  # TODO this should not be an in... should be a startswith and use params
                true_proteins.append(records.id)

        protein_data = [x[0].identifier for x in self.grouped_scored_proteins]
        false_true_positives = []
        decoys = []
        for i in range(len(protein_data)):
            # Get a list of false true positives from the protein data... basically if its not a decoy and if its not
            # in true_proteins
            if not protein_data[i].startswith(self.decoy_symbol) and protein_data[i] not in true_proteins:
                false_true_positives.append(protein_data[i])
            if protein_data[i].startswith(self.decoy_symbol):
                decoys.append(protein_data[i])

        entrapment_proteins_set = list(false_true_positives)
        entrapment_proteins_set = set(entrapment_proteins_set)

        # pick out the lead scoring protein for each group... lead score is at 0 position
        lead_score = [x[0] for x in self.grouped_scored_proteins]
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
            if fdr < false_discovery_rate:
                break
            else:
                # Here we delete the worst score every time... thus making our list smaller and smaller
                del lead_proteins[0]

        lead_proteins.reverse()

        fdr_restricted_set = [self.grouped_scored_proteins[x] for x in range(len(lead_proteins))]

        onehitwonders = []
        for groups in fdr_restricted_set:
            if int(groups[0].num_peptides) == 1:
                onehitwonders.append(groups[0])

        logger.info(
            "Protein Group leads that pass with more than 1 PSM with a {}  Entrapment FDR =  {}".format(
                false_discovery_rate, str(len(fdr_restricted_set) - len(onehitwonders))
            )
        )
        logger.info(
            "Protein Group lead One hit Wonders that pass {} Entrapment FDR = {}".format(
                false_discovery_rate, len(onehitwonders)
            )
        )

        logger.info(
            "Number of Protein groups that pass a {} Entrapment FDR: {}".format(
                str(false_discovery_rate * 100), len(fdr_restricted_set)
            )
        )

        return fdr_restricted_set

    def validate_psm_data(self):
        """
        Method that validates the PSM data
        """
        self._validate_decoys_from_data()
        self._validate_isoform_from_data()

    def validate_digest(self):
        """
        Method that validates the :py:class:`py_protein_inference.in_silico_digest.Digest` object
        """
        self._validate_reviewed_v_unreviewed()
        self._check_target_decoy_split()

    def check_data_consistency(self):
        """
        Method that checks for data consistency
        """
        self._check_data_digest_overlap_psms()
        self._check_data_digest_overlap_proteins()

    def _check_data_digest_overlap_psms(self):
        """
        Method that logs the overlap between the digested fasta file and the input files on the PSM level
        """
        peptides = [x.stripped_peptide for x in self.main_data_form]
        peptides_in_digest = set(self.digest.peptide_to_protein_dictionary.keys())
        peptides_from_search_in_digest = [x for x in peptides if x in peptides_in_digest]
        percentage = float(len(set(peptides))) / float(len(set(peptides_from_search_in_digest)))
        logger.info("{} PSMs identified from input files".format(len(peptides)))
        logger.info(
            "{} PSMs identified from input files that are also present in database digestion".format(
                len(peptides_from_search_in_digest)
            )
        )
        logger.info(
            "{}; ratio of PSMs identified from input files to those that are present in the search"
            " and in the database digestion".format(percentage)
        )

    def _check_data_digest_overlap_proteins(self):
        """
        Method that logs the overlap between the digested fasta file and the input files on the Protein level
        """
        proteins = [x.possible_proteins for x in self.main_data_form]
        flat_proteins = set([item for sublist in proteins for item in sublist])
        proteins_in_digest = set(self.digest.protein_to_peptide_dictionary.keys())
        proteins_from_search_in_digest = [x for x in flat_proteins if x in proteins_in_digest]
        percentage = float(len(flat_proteins)) / float(len(proteins_from_search_in_digest))
        logger.info("{} proteins identified from input files".format(len(flat_proteins)))
        logger.info(
            "{} proteins identified from input files that are also present in database digestion".format(
                len(proteins_from_search_in_digest)
            )
        )
        logger.info(
            "{}; ratio of proteins identified from input files that are also present in database digestion".format(
                percentage
            )
        )

    def _check_target_decoy_split(self):
        """
        Method that logs the number of target and decoy proteins from the digest
        """
        # Check the number of targets vs the number of decoys from the digest
        targets = [
            x
            for x in self.digest.protein_to_peptide_dictionary.keys()
            if self.parameter_file_object.decoy_symbol not in x
        ]
        decoys = [
            x for x in self.digest.protein_to_peptide_dictionary.keys() if self.parameter_file_object.decoy_symbol in x
        ]
        ratio = float(len(targets)) / float(len(decoys))
        logger.info("Number of Target Proteins in Digest: {}".format(len(targets)))
        logger.info("Number of Decoy Proteins in Digest: {}".format(len(decoys)))
        logger.info("Ratio of Targets Proteins to Decoy Proteins: {}".format(ratio))

    def _validate_score(self):
        # Make sure the score is actually in our data file header... Not sure if we can do this...?
        pass

    def _validate_decoys_from_data(self):
        """
        Method that checks to make sure that target and decoy proteins exist in the data files
        """
        # Check to see if we find decoys from our input files
        proteins = [x.possible_proteins for x in self.main_data_form]
        flat_proteins = set([item for sublist in proteins for item in sublist])
        targets = [x for x in flat_proteins if self.parameter_file_object.decoy_symbol not in x]
        decoys = [x for x in flat_proteins if self.parameter_file_object.decoy_symbol in x]
        logger.info("Number of Target Proteins in Data Files: {}".format(len(targets)))
        logger.info("Number of Decoy Proteins in Data Files: {}".format(len(decoys)))

    def _validate_isoform_from_data(self):
        """
        Method that validates whether or not isoforms are able to be identified in the data files
        """
        # Check to see if we find any proteins with isoform info in name in our input files
        proteins = [x.possible_proteins for x in self.main_data_form]
        flat_proteins = set([item for sublist in proteins for item in sublist])
        if self.parameter_file_object.isoform_symbol:
            non_iso = [x for x in flat_proteins if self.parameter_file_object.isoform_symbol not in x]

        else:
            non_iso = [x for x in flat_proteins]

        if self.parameter_file_object.isoform_symbol:
            iso = [x for x in flat_proteins if self.parameter_file_object.isoform_symbol in x]

        else:
            iso = []
        logger.info("Number of Non Isoform Labeled Proteins in Data Files: {}".format(len(non_iso)))
        logger.info("Number of Isoform Labeled Proteins in Data Files: {}".format(len(iso)))

    def _validate_reviewed_v_unreviewed(self):
        """
        Method that logs whether or not we can distinguish from reviewed and unreviewd protein identifiers in the digest
        """
        # Check to see if we get reviewed prots in digest...
        reviewed_proteins = len(self.digest.swiss_prot_protein_set)
        proteins_in_digest = len(set(self.digest.protein_to_peptide_dictionary.keys()))
        unreviewed_proteins = proteins_in_digest - reviewed_proteins
        logger.info("Number of Total Proteins in from Digest: {}".format(proteins_in_digest))
        logger.info("Number of Reviewed Proteins in from Digest: {}".format(reviewed_proteins))
        logger.info("Number of Unreviewed Proteins in from Digest: {}".format(unreviewed_proteins))

    @classmethod
    def sort_protein_strings(cls, protein_string_list, sp_proteins, decoy_symbol):
        """
        Method that sorts protein strings in the following order: Target Reviewed, Decoy Reviewed, Target Unreviewed,
         Decoy Unreviewed

        Args:
            protein_string_list (list): List of Protein Strings
            sp_proteins (set): Set of Reviewed Protein Strings
            decoy_symbol (str): Symbol to denote a decoy protein identifier IE "##"

        Returns:
            list: List of sorted protein strings

        Example:
            >>> list_of_group_objects = datastore.DataStore.sort_protein_strings(
            >>>     protein_string_list=protein_string_list, sp_proteins=sp_proteins, decoy_symbol="##"
            >>> )
        """

        our_target_sp_proteins = sorted([x for x in protein_string_list if x in sp_proteins and decoy_symbol not in x])
        our_decoy_sp_proteins = sorted([x for x in protein_string_list if x in sp_proteins and decoy_symbol in x])

        our_target_tr_proteins = sorted(
            [x for x in protein_string_list if x not in sp_proteins and decoy_symbol not in x]
        )
        our_decoy_tr_proteins = sorted([x for x in protein_string_list if x not in sp_proteins and decoy_symbol in x])

        identifiers_sorted = (
            our_target_sp_proteins + our_decoy_sp_proteins + our_target_tr_proteins + our_decoy_tr_proteins
        )

        return identifiers_sorted

    def input_has_q(self):
        """
        Method that checks to see if the input data has q values
        """
        len_q = len([x.qvalue for x in self.main_data_form if x.qvalue])
        len_all = len(self.main_data_form)
        if len_q == len_all:
            status = True
            logger.info("Input has Q value; Can restrict by Q value")
        else:
            status = False
            logger.warning("Input does not have Q value; Cannot restrict by Q value")

        return status

    def input_has_pep(self):
        """
        Method that checks to see if the input data has pep values
        """
        len_pep = len([x.pepvalue for x in self.main_data_form if x.pepvalue])
        len_all = len(self.main_data_form)
        if len_pep == len_all:
            status = True
            logger.info("Input has Pep value; Can restrict by Pep value")
        else:
            status = False
            logger.warning("Input does not have Pep value; Cannot restrict by Pep value")

        return status

    def input_has_custom(self):
        """
        Method that checks to see if the input data has custom score values
        """
        len_c = len([x.custom_score for x in self.main_data_form if x.custom_score])
        len_all = len(self.main_data_form)
        if len_c == len_all:
            status = True
            logger.info("Input has Custom value; Can restrict by Custom value")

        else:
            status = False
            logger.warning("Input does not have Custom value; Cannot restrict by Custom value")

        return status

    def get_protein_objects(self, fdr_restricted=False):
        if fdr_restricted:
            protein_objects = [
                x.proteins for x in self.protein_group_objects if x.q_value <= self.parameter_file_object.fdr
            ]
        else:
            protein_objects = self.grouped_scored_proteins

        return protein_objects

    def _init_validate(self, reader):
        if reader.psms:
            self.main_data_form = reader.psms  # Unrestricted PSM data
            self.restricted_peptides = [x.non_flanking_peptide for x in self.main_data_form]
        else:
            raise ValueError(
                "Psms variable from Reader object is either empty or does not exist. "
                "Make sure your files contain proper data and that you run the 'read_psms' "
                "method on your Reader object."
            )

    def _validate_main_data_form(self):
        if self.main_data_form:
            pass
        else:
            raise ValueError(
                "Main Data is not defined, thus method cannot be ran. Please make sure PSM data is properly"
                " loaded from the Reader object"
            )

    def _validate_main_data_restricted(self):
        if self.main_data_restricted:
            pass
        else:
            raise ValueError(
                "Main Data Restricted is not defined, thus method cannot be ran. Please make sure PSM data is properly"
                " loaded from the Reader object and make sure to run DataStore method 'restrict_psm_data'."
            )

    def _validate_scored_proteins(self):
        if self.picked_proteins_scored or self.scored_proteins:
            pass
        else:
            raise ValueError(
                "Proteins have not been scored, Please initialize a Score object and run a score method with"
                " 'score_psms' instance method."
            )

    def _validate_scoring_input(self):
        if self.scoring_input:
            pass
        else:
            raise ValueError(
                "Scoring input has not been created, Please run 'create_scoring_input' method from the DataStore "
                "object to continue."
            )

    def _validate_protein_group_objects(self):
        if self.protein_group_objects and self.grouped_scored_proteins:
            pass
        else:
            raise ValueError(
                "Either 'protein_group_objects' or 'grouped_scored_proteins' or both DataStore variables are undefined."
                " Please make sure you run an inference method from the Inference class before proceeding."
            )

    def generate_fdr_vs_target_hits(self, fdr_max=0.2):

        fdr_vs_count = []
        count_list = []
        for pg in self.protein_group_objects:
            if not pg.proteins[0].identifier.startswith(self.decoy_symbol):
                count_list.append(pg)
            fdr_vs_count.append([pg.q_value, len(count_list)])

        fdr_vs_count = [x for x in fdr_vs_count if x[0] < fdr_max]

        return fdr_vs_count
