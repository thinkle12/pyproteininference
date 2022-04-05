from unittest import TestCase

from pyproteininference.physical import Psm


class TestLoadPeptideParsing(TestCase):
    def test_peptide_parsing(self):

        # GFY Sequence Flanking
        peptide1 = "R.ALASLM*TYK.-"
        peptide1_sf = Psm.split_peptide(peptide_string=peptide1)
        self.assertEqual("ALASLM*TYK", peptide1_sf)

        peptide1_no_mods = Psm.remove_peptide_mods(peptide_string=peptide1_sf)
        self.assertEqual("ALASLMTYK", peptide1_no_mods)

        # GFY Sequence no Flanking
        peptide2 = "ALASLM*TYK"
        peptide2_sf = Psm.split_peptide(peptide_string=peptide2)
        self.assertEqual("ALASLM*TYK", peptide2_sf)

        peptide2_no_mods = Psm.remove_peptide_mods(peptide_string=peptide2_sf)
        self.assertEqual("ALASLMTYK", peptide2_no_mods)

        # Max Quant Sequence Flanking _ and string mod
        peptide3 = "_ALASLM(ox)TYK_"
        peptide3_sf = Psm.split_peptide(peptide_string=peptide3, delimiter="_")
        self.assertEqual("ALASLM(ox)TYK", peptide3_sf)

        peptide3_no_mods = Psm.remove_peptide_mods(peptide_string=peptide3_sf)
        self.assertEqual("ALASLMTYK", peptide3_no_mods)

        # Max Quant with float mass mod
        peptide4 = "ALASLM(137.9)TYK"
        peptide4_sf = Psm.split_peptide(peptide_string=peptide4)
        self.assertEqual("ALASLM(137.9)TYK", peptide4_sf)

        peptide4_no_mods = Psm.remove_peptide_mods(peptide_string=peptide4_sf)
        self.assertEqual("ALASLMTYK", peptide4_no_mods)

        # Max Quant sequence with int mass mod
        peptide5 = "ALASLM(1)TYK"
        peptide5_sf = Psm.split_peptide(peptide_string=peptide5)
        self.assertEqual("ALASLM(1)TYK", peptide5_sf)

        peptide5_no_mods = Psm.remove_peptide_mods(peptide_string=peptide5_sf)
        self.assertEqual("ALASLMTYK", peptide5_no_mods)

        # Max Quant with float mass and flanking with "."
        peptide6 = "R.ALASLM(137.9)TYK.-"
        peptide6_sf = Psm.split_peptide(peptide_string=peptide6)
        self.assertEqual("ALASLM(137.9)TYK", peptide6_sf)

        peptide6_no_mods = Psm.remove_peptide_mods(peptide_string=peptide6_sf)
        self.assertEqual("ALASLMTYK", peptide6_no_mods)

        # Comet Pepxml with brackets and int mod and n term
        peptide7 = "n[230]M[147]PVLSAR"
        peptide7_sf = Psm.split_peptide(peptide_string=peptide7)
        self.assertEqual("n[230]M[147]PVLSAR", peptide7_sf)

        peptide7_no_mods = Psm.remove_peptide_mods(peptide_string=peptide7_sf)
        self.assertEqual("MPVLSAR", peptide7_no_mods)

        # Blank sequence
        peptide8 = "MPVLSAR"
        peptide8_sf = Psm.split_peptide(peptide_string=peptide8)
        self.assertEqual("MPVLSAR", peptide8_sf)

        peptide8_no_mods = Psm.remove_peptide_mods(peptide_string=peptide8_sf)
        self.assertEqual("MPVLSAR", peptide8_no_mods)

        # Blank sequence Mascot
        peptide9 = "MLSVGLGFLR"
        peptide9_sf = Psm.split_peptide(peptide_string=peptide9)
        self.assertEqual("MLSVGLGFLR", peptide9_sf)

        peptide9_no_mods = Psm.remove_peptide_mods(peptide_string=peptide9_sf)
        self.assertEqual("MLSVGLGFLR", peptide9_no_mods)

        # Peaks, no flanking and float mod within parentheses
        peptide10 = "VM(+15.99)YKFLTV"
        peptide10_sf = Psm.split_peptide(peptide_string=peptide10)
        self.assertEqual("VM(+15.99)YKFLTV", peptide10_sf)

        peptide10_no_mods = Psm.remove_peptide_mods(peptide_string=peptide10_sf)
        self.assertEqual("VMYKFLTV", peptide10_no_mods)

        # Peaks, with flanking and float mod within parentheses
        peptide11 = "R.VM(+15.99)YKFLTV.-"
        peptide11_sf = Psm.split_peptide(peptide_string=peptide11)
        self.assertEqual("VM(+15.99)YKFLTV", peptide11_sf)

        peptide11_no_mods = Psm.remove_peptide_mods(peptide_string=peptide11_sf)
        self.assertEqual("VMYKFLTV", peptide11_no_mods)

        # MSGF+ sequence with no flanking and no parentheses or brackets
        peptide12 = "NLANPTSVILASIQM+15.995LEYLGMADK"
        peptide12_sf = Psm.split_peptide(peptide_string=peptide12)
        self.assertEqual("NLANPTSVILASIQM+15.995LEYLGMADK", peptide12_sf)

        peptide12_no_mods = Psm.remove_peptide_mods(peptide_string=peptide12_sf)
        self.assertEqual("NLANPTSVILASIQMLEYLGMADK", peptide12_no_mods)

        # MSGF+ with flanking and parentheses and no brackets
        peptide13 = "K.NLANPTSVILASIQM+15.995LEYLGMADK.A"
        peptide13_sf = Psm.split_peptide(peptide_string=peptide13)
        self.assertEqual("NLANPTSVILASIQM+15.995LEYLGMADK", peptide13_sf)

        peptide13_no_mods = Psm.remove_peptide_mods(peptide_string=peptide13_sf)
        self.assertEqual("NLANPTSVILASIQMLEYLGMADK", peptide13_no_mods)
