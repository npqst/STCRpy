import unittest
import glob

import stcrpy
from stcrpy.tcr_processing import TCRParser, abTCR, TCR, MHCchain, MHC


class TestTCRParser(unittest.TestCase):

    def test_imports(self):
        import stcrpy
        from stcrpy.tcr_processing import TCRParser, abTCR, TCR, MHCchain, MHC

    def test_get_tcr_structure_class_I(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./test_files/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["DE", "IJ"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(["FG", "AB"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_antigens()]) == set(
            ["C", "H"]
        )

    def test_get_tcr_structure_class_II(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./test_files/6r0e.cif"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["DE"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(["AB"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_antigens()]) == set(["C"])

    def test_all_stcrdab(self):
        from tqdm import tqdm

        with open("./test_files/tcr_pdb_codes.txt") as f:
            pdb_codes = f.readlines()
        pdb_codes = [x.strip() for x in pdb_codes]
        badly_parsed_pdb = []
        errors = {}
        pdb_types = {}
        for pdb_code in tqdm(pdb_codes):
            # pdb_id = pdb_file.split("/")[-1].split(".")[0]
            try:
                tcr = stcrpy.fetch_TCR(pdb_code)
            except Exception as e:
                errors[pdb_code] = e
        print(errors)
        assert len(badly_parsed_pdb) == 0

    def test_delta_beta_tcr_parsed_as_abTCR(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./test_files/DB_test_T104_rank_0_model_0_refined.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["AB"])
        assert all([isinstance(x, abTCR) for x in tcr.get_TCRs()])

    def test_save(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./test_files/4nhu.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)

        from stcrpy.tcr_processing.TCRIO import TCRIO

        io = TCRIO()

        for x in tcr.get_TCRs():
            io.save(x, save_as=f"./test_files/test_{x.id}_TCR_only.pdb")

        for x in tcr.get_TCRs():
            io.save(x, tcr_only=True, save_as=f"./test_files/test_{x.id}.pdb")

        pdb_file = (
            "../stcrpy/tcr_geometry/reference_data/dock_reference_1_imgt_numbered.pdb"
        )
        tcr = parser.get_tcr_structure("test", pdb_file)
        for x in tcr.get_TCRs():
            io.save(x, save_as=f"./test_files/test_{x.id}.pdb")

    def test_error_prone_tcrs(self):
        parser = TCRParser.TCRParser()
        pdb_files = glob.glob("./test_files/TCRParser_test_files/*")
        for file in pdb_files:
            pdb_id = file.split("/")[-1].split(".")[0]
            print(pdb_id)
            tcr_structure = parser.get_tcr_structure(pdb_id, file)
            for tcr in tcr_structure.get_TCRs():
                assert isinstance(tcr, TCR)

    def test_MHC_single_chain_handling(self):
        with open("./test_files/tcr_pdb_codes.txt") as f:
            pdb_codes = f.readlines()
        pdb_codes = [x.strip() for x in pdb_codes]

        badly_parsed_pdb = []
        errors = {}
        single_chain_MHC = {}
        apo_TCRs = {}
        for pdb_file in pdb_codes:
            pdb_id = pdb_file.split("/")[-1].split(".")[0]
            try:
                tcr = stcrpy.fetch_TCR(pdb_id)
                if len(list(tcr.get_TCRs())) == 0:
                    badly_parsed_pdb.append(pdb_id)
                else:
                    for t in tcr.get_TCRs():
                        if len(t.get_MHC()) == 0:
                            apo_TCRs[f"{pdb_id}_{t.id}"] = t
                            print(pdb_id, "No MHC found")
                            continue
                        mhc = t.get_MHC()
                        print(pdb_id, mhc)
                        if isinstance(mhc[0], MHCchain):
                            single_chain_MHC[pdb_id] = t
            except Exception as e:
                errors[pdb_id] = e
        print(badly_parsed_pdb)
        print(len(badly_parsed_pdb))
        print(single_chain_MHC)

    def test_MHC_association(self):
        with open("./test_files/tcr_pdb_codes.txt") as f:
            pdb_codes = f.readlines()
        pdb_codes = [x.strip() for x in pdb_codes]

        badly_parsed_pdb = []
        errors = {}
        apo_TCRs = {}
        true_apo = [
            "1hxm",
            "1kb5",
            "1kgc",
            "1nfd",
            "1tcr",
            "2bnu",
            "2cde",
            "2cdf",
            "2cdg",
            "2eyr",
            "2eys",
            "2eyt",
            "2ial",
        ]
        for pdb_file in pdb_codes:
            pdb_id = pdb_file.split("/")[-1].split(".")[0]
            try:
                tcr = stcrpy.fetch_TCR(pdb_id)
                if len(list(tcr.get_TCRs())) == 0:
                    badly_parsed_pdb.append(pdb_id)
                else:
                    for t in tcr.get_TCRs():
                        if len(t.get_MHC()) == 0 and any(
                            [isinstance(x, MHC) for x in tcr.child_list]
                        ):
                            apo_TCRs[f"{pdb_id}_{t.id}"] = t
                            print(pdb_id, "No MHC found")
                            continue
            except Exception as e:
                errors[pdb_id] = e
        print(badly_parsed_pdb)
        print(len(badly_parsed_pdb))

    def test_MR1_parsing(self):
        import stcrpy

        tcr1, tcr2 = stcrpy.fetch_TCR("5d7i")
        tcr1.get_MHC()[0]
        tcr2.get_MHC()[0]

        tcr1, tcr2 = stcrpy.fetch_TCR("4pjf")
        tcr1.get_MHC()[0]
        tcr2.get_MHC()[0]

    def test_scMH2_parsing(self):
        import stcrpy

        with self.assertWarns(UserWarning):
            tcr = stcrpy.fetch_TCR("6u3n")
            # should raise warning saying that other MHC Class II chain is missing
        tcr.get_MHC()[0].get_MHC_type()

        with self.assertWarns(UserWarning):
            tcr = stcrpy.fetch_TCR("6mkr")
            # should raise warning saying that other MHC Class II chain is missing
        tcr.get_MHC()[0].get_MHC_type()
