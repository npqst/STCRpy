import unittest
import os

import stcrpy
from stcrpy.tcr_processing import TCRParser, abTCR

STRUCTURES = "./test_files/structures"
PARSER_FILES = "./test_files/TCRParser_test_files"


class TestTCRParser(unittest.TestCase):

    def test_imports(self):
        from stcrpy.tcr_processing import TCRParser, abTCR

    def test_get_tcr_structure_class_I(self):
        parser = TCRParser.TCRParser()

        pdb_file = f"{PARSER_FILES}/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["DE", "IJ"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(["FG", "AB"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_antigens()]) == set(
            ["C", "H"]
        )

    def test_get_tcr_structure_class_II(self):
        parser = TCRParser.TCRParser()

        pdb_file = f"{PARSER_FILES}/6r0e.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["DE"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(["AB"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_antigens()]) == set(["C"])

    def test_delta_beta_tcr_parsed_as_abTCR(self):
        tcr = stcrpy.load_TCRs(f"{STRUCTURES}/6vrn.cif")
        assert all([isinstance(x, abTCR) for x in tcr])

    def test_save(self):
        tcr = stcrpy.load_TCRs(f"{STRUCTURES}/4nhu.cif")

        from stcrpy.tcr_processing.TCRIO import TCRIO

        io = TCRIO()

        if os.path.exists("./test_files/out/tcr_processing/test_BA_TCR_only.pdb"):
            os.remove("./test_files/out/tcr_processing/test_BA_TCR_only.pdb")
        if os.path.exists("./test_files/out/tcr_processing/test_DC_TCR_only.pdb"):
            os.remove("./test_files/out/tcr_processing/test_DC_TCR_only.pdb")
        if os.path.exists("./test_files/out/tcr_processing/test_BA.pdb"):
            os.remove("./test_files/out/tcr_processing/test_BA.pdb")
        if os.path.exists("./test_files/out/tcr_processing/test_DC.pdb"):
            os.remove("./test_files/out/tcr_processing/test_DC.pdb")

        for x in tcr:
            io.save(
                x, save_as=f"./test_files/out/tcr_processing/test_{x.id}_TCR_only.pdb"
            )

        for x in tcr:
            io.save(
                x,
                tcr_only=True,
                save_as=f"./test_files/out/tcr_processing/test_{x.id}.pdb",
            )

        self.assertTrue(
            os.path.exists("./test_files/out/tcr_processing/test_BA_TCR_only.pdb")
            and os.path.exists("./test_files/out/tcr_processing/test_DC_TCR_only.pdb")
            and os.path.exists("./test_files/out/tcr_processing/test_BA.pdb")
            and os.path.exists("./test_files/out/tcr_processing/test_DC.pdb")
        )

    def test_MR1_parsing(self):
        tcr1, tcr2 = stcrpy.load_TCRs(f"{STRUCTURES}/5d7i.cif")
        tcr1.get_MHC()[0]
        tcr2.get_MHC()[0]
        assert tcr2.get_MHC()[0].get_MHC_type() == "MR1"

        tcr1, tcr2 = stcrpy.load_TCRs(f"{STRUCTURES}/4pjf.cif")
        tcr1.get_MHC()[0]
        tcr2.get_MHC()[0]
        assert tcr2.get_MHC()[0].get_MHC_type() == "MR1"

    def test_scMH2_parsing(self):
        with self.assertWarns(UserWarning):
            tcrs = stcrpy.load_TCRs(f"{STRUCTURES}/6u3n.cif")
        tcrs[0].get_MHC()[0].get_MHC_type()

        with self.assertWarns(UserWarning):
            tcrs = stcrpy.load_TCRs(f"{STRUCTURES}/6mkr.cif")
        tcrs[0].get_MHC()[0].get_MHC_type()


class TestTCR(unittest.TestCase):
    def test_crop_class_I(self):
        tcrs = stcrpy.load_TCR(f"{PARSER_FILES}/5hyj.pdb")
        tcr = tcrs[0]

        tcr.crop()

        self.assertTrue(max([res.id[1] for chain in tcr for res in chain]) <= 128)
        self.assertTrue(max([res.id[1] for chain in tcr.MHC[0] for res in chain if res.id[0] == ' ']) % 1000 <= 92)

    def test_crop_class_I_dont_remove_hetatoms(self):
        tcrs = stcrpy.load_TCR(f"{PARSER_FILES}/5hyj.pdb")
        tcr = tcrs[0]

        tcr.crop(remove_het_atoms=False)

        self.assertTrue(max([res.id[1] for chain in tcr for res in chain if res.id[0] == ' ']) <= 128)
        self.assertFalse(max([res.id[1] for chain in tcr for res in chain]) <= 128)

    def test_crop_class_II(self):
        tcr = stcrpy.load_TCR(f"{PARSER_FILES}/8vcy_class_II.pdb")
        tcr.crop()

        self.assertTrue(max([res.id[1] for chain in tcr for res in chain if res.id[0] == ' ']) <= 128)
        self.assertTrue(max([res.id[1] for chain in tcr.MHC[0] for res in chain if res.id[0] == ' ']) <= 92)
        tcr.get_MHC()[0].get_MHC_type()


class TestabTCR(unittest.TestCase):
    def setUp(self):
        pdb_file = f"{PARSER_FILES}/5hyj.pdb"
        self.tcrs = stcrpy.load_TCR(pdb_file)

    def test_standardise_chain_names(self):
        tcr = self.tcrs[1]

        tcr.standardise_chain_names()

        self.assertEqual(tcr.id, 'ED')
        self.assertEqual(tcr.VA, 'D')
        self.assertEqual(tcr.get_VA().id, 'D')
        self.assertEqual(tcr.VB, 'E')
        self.assertEqual(tcr.get_VB().id, 'E')

        self.assertEqual(tcr.antigen[0].id, 'C')

        self.assertEqual(tcr.MHC[0].id, 'AB')
        self.assertEqual(tcr.MHC[0].get_alpha().id, 'A')
        self.assertEqual(tcr.MHC[0].get_B2M().id, 'B')

    def test_standardise_chain_names_when_already_standard(self):
        tcr = self.tcrs[0]

        tcr.standardise_chain_names()

        self.assertEqual(tcr.id, 'ED')
        self.assertEqual(tcr.VA, 'D')
        self.assertEqual(tcr.get_VA().id, 'D')
        self.assertEqual(tcr.VB, 'E')
        self.assertEqual(tcr.get_VB().id, 'E')

        self.assertEqual(tcr.antigen[0].id, 'C')

        self.assertEqual(tcr.MHC[0].id, 'AB')
        self.assertEqual(tcr.MHC[0].get_alpha().id, 'A')
        self.assertEqual(tcr.MHC[0].get_B2M().id, 'B')

    def test_standardise_chain_names_mhc_chain_e(self):
        tcr, _ = stcrpy.load_TCR(f"{PARSER_FILES}/2e7l.pdb")

        tcr.standardise_chain_names()

        self.assertEqual(tcr.id, 'ED')
        self.assertEqual(tcr.VA, 'D')
        self.assertEqual(tcr.get_VA().id, 'D')
        self.assertEqual(tcr.VB, 'E')
        self.assertEqual(tcr.get_VB().id, 'E')

        self.assertEqual(tcr.antigen[0].id, 'C')

        self.assertEqual(tcr.MHC[0].id, 'A')
        self.assertEqual(tcr.MHC[0].get_alpha().id, 'A')

    def test_swapped_chain_ids(self):
        tcr = stcrpy.load_TCR(f"{PARSER_FILES}/swapped_chain_ids_5hyj_DECAB.pdb")
        tcr.standardise_chain_names()

        self.assertEqual(tcr.id, 'ED')
        self.assertEqual(tcr.VA, 'D')
        self.assertEqual(tcr.get_VA().id, 'D')
        self.assertEqual(tcr.VB, 'E')
        self.assertEqual(tcr.get_VB().id, 'E')

        self.assertEqual(tcr.antigen[0].id, 'C')

        self.assertEqual(tcr.MHC[0].id, 'AB')
        self.assertEqual(tcr.MHC[0].get_alpha().id, 'A')
        self.assertEqual(tcr.MHC[0].get_B2M().id, 'B')


class TestgdTCR(unittest.TestCase):
    def test_standardise_chain_names(self):
        pdb_file = f"{PARSER_FILES}/1hxm_dbTCRs.pdb"
        tcrs = stcrpy.load_TCR(pdb_file)
        tcr = tcrs[0]

        tcr.standardise_chain_names()

        self.assertEqual(tcr.id, 'ED')
        self.assertEqual(tcr.VD, 'D')
        self.assertEqual(tcr.get_VD().id, 'D')
        self.assertEqual(tcr.VG, 'E')
        self.assertEqual(tcr.get_VG().id, 'E')


class TestSymmetryMates(unittest.TestCase):
    def test_symmetry_mates(self):
        try:
            import pymol
        except ImportError:
            return

        tcrs = stcrpy.load_TCRs(f"{STRUCTURES}/6ulr.cif")
        assert len(tcrs[0].get_MHC()) == 1
