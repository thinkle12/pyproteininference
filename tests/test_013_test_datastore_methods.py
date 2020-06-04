import tempfile
from unittest import TestCase
from pkg_resources import resource_filename

import protein_inference

TEST_DATABASE = resource_filename('protein_inference', '../tests/data/test_database.fasta')
PARAMETER_FILE = resource_filename('protein_inference', '../tests/data/test_params_inclusion.yaml')
OUTPUT_DIR = tempfile.gettempdir()
# OUTPUT_DIR = resource_filename('protein_inference', '../tests/output/')

TARGET_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_target.txt')
DECOY_FILE = resource_filename('protein_inference', '../tests/data/test_perc_data_decoy.txt')

temp_dir = tempfile.gettempdir()

class TestDataStoreMethods(TestCase):

    def test_datastore(self):

        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        ### STEP 1: Load parameter file ###
        protein_inference_parameters = protein_inference.parameters.ProteinInferenceParameter(yaml_param_filepath=PARAMETER_FILE)

        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        ### STEP 2: Start with running an In Silico Digestion ###
        digest = protein_inference.in_silico_digest.InSilicoDigest(database_path=TEST_DATABASE,
                                                 parameter_file_object=protein_inference_parameters, id_splitting=True)
        digest.digest_fasta_database()

        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        ### STEP 3: Read PSM Data ###
        pep_and_prot_data = protein_inference.reader.GenericReader(target_file=TARGET_FILE,
                                                                   decoy_file=DECOY_FILE,
                                                                   combined_files=None,
                                                                   parameter_file_object=protein_inference_parameters,
                                                                   digest_class=digest)
        pep_and_prot_data.read_psms()

        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        ### STEP 4: Initiate the datastore class ###
        data = protein_inference.datastore.DataStore(pep_and_prot_data, digest_class=digest)

        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        ### Step 5: Restrict the PSM data
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters)

        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        ### Step 6: Generate protein scoring input
        data.create_scoring_input(score_input=protein_inference_parameters.score)

        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        ### Step 7: Remove non unique peptides if running exclusion
        if protein_inference_parameters.inference_type == "exclusion":
            # This gets ran if we run exclusion...
            data.exclude_non_distinguishing_peptides(digest_class=digest)

        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        ### STEP 8: Score our PSMs given a score method
        score = protein_inference.scoring.Score(data_class=data)
        score.score_psms(score_method=protein_inference_parameters.score_method)

        ### STEP 9: Run protein picker on the data
        ### STEP 9: Run protein picker on the data
        ### STEP 9: Run protein picker on the data
        if protein_inference_parameters.picker:
            data.protein_picker()
        else:
            pass

        ### STEP 10: Apply Inference
        ### STEP 10: Apply Inference
        ### STEP 10: Apply Inference
        inference_type = protein_inference_parameters.inference_type

        # For parsimony... Run GLPK setup, runner, grouper...
        if inference_type == 'parsimony':
            group = protein_inference.inference.Parsimony(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "inclusion":
            group = protein_inference.inference.Inclusion(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "exclusion":
            group = protein_inference.inference.Exclusion(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "none":
            group = protein_inference.inference.FirstProtein(data_class=data, digest_class=digest)
            group.infer_proteins()

        if inference_type == "peptide_centric":
            group = protein_inference.inference.PeptideCentric(data_class=data, digest_class=digest)
            group.infer_proteins()

        ### STEP 11: Run FDR and Q value Calculations
        ### STEP 11: Run FDR and Q value Calculations
        ### STEP 11: Run FDR and Q value Calculations
        data.set_based_fdr(false_discovery_rate=float(protein_inference_parameters.fdr))
        data.calculate_q_values()

        # Start datastore tests

        scored_identifiers = data.get_sorted_identifiers(digest_class=digest,scored=True)
        self.assertListEqual(scored_identifiers,['ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056', 'HNRPU_HUMAN|Q00839', 'RAF1_HUMAN|P04049', 'RPOC_SHIF8|Q0SY12', 'TCAF1_HUMAN|Q9Y4C2', '##TCAF2_HUMAN|##A6NFQ2', 'B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7'])

        identifiers = data.get_sorted_identifiers(digest_class=digest,scored=False)
        self.assertListEqual(identifiers, ['ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056', 'HNRPU_HUMAN|Q00839', 'RAF1_HUMAN|P04049', 'RPOC_SHIF8|Q0SY12', 'TCAF1_HUMAN|Q9Y4C2', '##TCAF1_HUMAN|##Q9Y4C2', '##TCAF2_HUMAN|##A6NFQ2', 'B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7'])

        psm_data = data.get_psm_data()
        self.assertEqual(len(psm_data), 26)

        protein_data = data.get_protein_data()
        self.assertEqual(len(protein_data), 9)

        proteins_from_psms = data.get_protein_identifiers_from_psm_data()
        self.assertEqual(len(proteins_from_psms), 26)

        q_values = data.get_q_values()
        self.assertEqual(len(q_values), 26)

        pep_values = data.get_pep_values()
        self.assertEqual(len(pep_values), 26)

        protein_info_dict = data.get_protein_information_dictionary()
        self.assertDictEqual(protein_info_dict['RAF1_HUMAN|P04049'][0], {'peptide': 'R.CQTCGYKFHEHCSTK.V', 'Qvalue': 0.00032, 'PosteriorErrorProbability': 3.5e-06, 'Percscore': 7.0})

        # Restrict with all filteres
        protein_inference_parameters.restrict_pep=.05
        protein_inference_parameters.restrict_q=.001
        protein_inference_parameters.restrict_peptide_length=9

        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 8)

        protein_inference_parameters.restrict_pep=None
        # Restrict length and q
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 8)

        protein_inference_parameters.restrict_q=None
        # Restrict length only
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 23)

        protein_inference_parameters.restrict_peptide_length=None
        # No Restriction
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 27)

        # No restriction and no removing 1 pep values
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=False)
        self.assertEqual(len(data.main_data_restricted), 27)

        protein_inference_parameters.restrict_pep=0.05
        protein_inference_parameters.restrict_q=0.05
        # Restrict pep and q without length
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 24)

        protein_inference_parameters.restrict_q=.001
        protein_inference_parameters.restrict_pep=None
        protein_inference_parameters.restrict_peptide_length=None
        # restrict pep and length without q
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 9)

        protein_inference_parameters.restrict_q=None
        protein_inference_parameters.restrict_pep=0.003
        protein_inference_parameters.restrict_peptide_length=None
        # restrict pep and length without q
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 24)

        protein_inference_parameters.restrict_q = None
        protein_inference_parameters.restrict_pep=0.0004
        protein_inference_parameters.restrict_peptide_length=9
        # restrict pep and length without q
        data.restrict_psm_data(parameter_file_object=protein_inference_parameters, remove1pep=True)
        self.assertEqual(len(data.main_data_restricted), 10)

        prot_dict = data.protein_to_peptide_dictionary()
        self.assertSetEqual(prot_dict['RAF1_HUMAN|P04049'], set(['FQMFQLIDIAR', 'WHGDVAVKILK', 'SASEPSLHR', 'CQTCGYKFHEHCSTK', 'QTAQGMDYLHAK']))
        pep_dict = data.peptide_to_protein_dictionary()
        self.assertSetEqual(pep_dict['QTAQGMDYLHAK'],set(['ARAF_HUMAN|P10398', 'RAF1_HUMAN|P04049', 'BRAF_HUMAN|P15056']))

        utl = data.unique_to_leads_peptides()
        self.assertSetEqual(utl, set(['FQMFQLIDIAR', 'TFFSLAFCDFCLK', 'WHGDVAVKILK', 'CGVEVTQTK', 'RVDYSGR', 'VTAEDVLKPGTADILVPR', 'IALASPDMIR', 'CQTCGYKFHEHCSTK', 'YCWMSTGLYIPGR', 'NTLLHEQWCDLLEENSVDAVK', 'VADLFEAR', 'VIDIWAAANDR', 'MGAEAIQALLK', 'IPQESGGTK', 'LIPAGTGYAYHQDR', 'MEPTPVPFCGAK', 'FATSDLNDLYR', 'LYLLTQMPH', 'EGLNVLQY#FISTHGAR']))
        hl = data.higher_or_lower()
        self.assertEqual(hl,"higher")

        df_main = data.get_protein_identifiers(data_form='main')
        self.assertListEqual(sorted(df_main[2][0]),sorted(['ARAF_HUMAN|P10398', 'RAF1_HUMAN|P04049', 'BRAF_HUMAN|P15056']))
        df_restricted = data.get_protein_identifiers(data_form='restricted')
        self.assertEqual(len(df_restricted),10)
        df_picked = data.get_protein_identifiers(data_form='picked')
        self.assertListEqual(df_picked,['RPOC_SHIF8|Q0SY12', 'RAF1_HUMAN|P04049', 'ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056', 'TCAF1_HUMAN|Q9Y4C2', 'B3KX72_HUMAN|B3KX72', 'HNRPU_HUMAN|Q00839', 'Q96BA7_HUMAN|Q96BA7', '##TCAF2_HUMAN|##A6NFQ2'])
        df_picked_removed = data.get_protein_identifiers(data_form='picked_removed')
        self.assertListEqual(df_picked_removed,['##TCAF1_HUMAN|##Q9Y4C2'])
        df_fdr_restricted = data.get_protein_identifiers(data_form='fdr_restricted')
        self.assertListEqual(df_fdr_restricted,['RPOC_SHIF8|Q0SY12', 'RAF1_HUMAN|P04049', 'ARAF_HUMAN|P10398', 'BRAF_HUMAN|P15056', 'TCAF1_HUMAN|Q9Y4C2', 'HNRPU_HUMAN|Q00839', 'B3KX72_HUMAN|B3KX72', 'Q96BA7_HUMAN|Q96BA7'])


        p1 = data.get_protein_information(protein_string="RAF1_HUMAN|P04049")
        self.assertListEqual(p1[1],['RAF1_HUMAN|P04049', 70.7434325345954, {2}, True, {'SASEPSLHR', 'QTAQGMDYLHAK', 'VFLPNKQR', 'WHGDVAVKILK', 'FQMFQLIDIAR', 'CQTCGYKFHEHCSTK'}, None, True, 6])



