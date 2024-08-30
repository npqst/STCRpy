import unittest
import glob

from STCRpy.tcr_processing import TCRParser, abTCR, TCR, MHCchain, MHC


class TestTCRParser(unittest.TestCase):
    def test_get_tcr_structure_class_I(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./STCRpy/test/test_files/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["DE", "IJ"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(["FG", "AB"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_antigens()]) == set(
            ["C", "H"]
        )

    def test_get_tcr_structure_class_II(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./STCRpy/test/test_files/6r0e.cif"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["DE"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(["AB"])
        assert set(["".join(sorted(x.id)) for x in tcr.get_antigens()]) == set(["C"])

    def test_all_stcrdab(self):
        import glob

        parser = TCRParser.TCRParser()
        stcrdab_pdb_files = glob.glob(
            "/home/quast/Projects/STCRDab/Data/entries/*/structure/imgt/*.pdb"
        )
        stcrdab_pdb_files.sort()
        badly_parsed_pdb = []
        errors = {}
        pdb_types = {}
        for pdb_file in stcrdab_pdb_files:
            pdb_id = pdb_file.split("/")[-1].split(".")[0]
            try:
                tcr = parser.get_tcr_structure(pdb_id, pdb_file)
                if len(list(tcr.get_TCRs())) == 0:
                    badly_parsed_pdb.append(pdb_id)
                else:
                    pdb_types[pdb_id] = (
                        type(list(tcr.get_TCRs())[0]),
                        list([str(x) for x in tcr.get_TCRs()]),
                    )
            except Exception as e:
                errors[pdb_id] = e
        print(badly_parsed_pdb)
        print(len(badly_parsed_pdb))

    def test_docked_tcr(self):
        import glob

        parser = TCRParser.TCRParser()
        dock_pdb_files = glob.glob(
            "/home/quast/Collaborations/SamuelsLab_Weizmann/docking/results/N17.2/310569-N17-2_NRAS_rank_0/structures/it1/renumbered_complex_*.pdb"
        )
        dock_pdb_files.sort()
        badly_parsed_pdb = []
        for pdb_file in dock_pdb_files:
            pdb_id = pdb_file.split("/")[-1].split(".")[0]
            tcr = parser.get_tcr_structure(pdb_id, pdb_file)
            if len(list(tcr.get_TCRs())) == 0:
                badly_parsed_pdb.append(pdb_id)
        print(badly_parsed_pdb)
        print(len(badly_parsed_pdb))

    def test_antigen_only(self):
        parser = TCRParser.TCRParser()
        pdb_file = "/home/quast/Collaborations/SamuelsLab_Weizmann/data/antigen_xtal_structures/7mle_1_aligned.cif"
        tcr = parser.get_tcr_structure("7mle", pdb_file)
        print(tcr)

    def test_delta_beta_tcr_parsed_as_abTCR(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./STCRpy/test/test_files/DB_test_T104_rank_0_model_0_refined.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        assert set(["".join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(["AB"])
        assert all([isinstance(x, abTCR) for x in tcr.get_TCRs()])

    def test_save(self):
        parser = TCRParser.TCRParser()

        pdb_file = "./STCRpy/test/test_files/4nhu.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)

        from ..STCRpy.tcr_processing.TCRIO import TCRIO

        io = TCRIO()

        for x in tcr.get_TCRs():
            io.save(x, save_as=f"./STCRpy/test/test_files/test_{x.id}_TCR_only.pdb")

        for x in tcr.get_TCRs():
            io.save(
                x, tcr_only=True, save_as=f"./STCRpy/test/test_files/test_{x.id}.pdb"
            )

        pdb_file = "./STCRpy/STCRpy/tcr_geometry/reference_data/dock_reference_1_imgt_numbered.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        for x in tcr.get_TCRs():
            io.save(x, save_as=f"./STCRpy/test/test_files/test_{x.id}.pdb")

    def test_error_prone_tcrs(self):
        parser = TCRParser.TCRParser()
        pdb_files = glob.glob("./STCRpy/test/test_files/TCRParser_test_files/*")
        for file in pdb_files:
            pdb_id = file.split("/")[-1].split(".")[0]
            print(pdb_id)
            tcr_structure = parser.get_tcr_structure(pdb_id, file)
            for tcr in tcr_structure.get_TCRs():
                assert isinstance(tcr, TCR)

    def test_MHC_single_chain_handling(self):
        import glob

        parser = TCRParser.TCRParser()
        stcrdab_pdb_files = glob.glob(
            "/home/quast/Projects/STCRDab/Data/entries/*/structure/imgt/*.pdb"
        )
        stcrdab_pdb_files.sort()
        badly_parsed_pdb = []
        errors = {}
        single_chain_MHC = {}
        apo_TCRs = {}
        for pdb_file in stcrdab_pdb_files:
            pdb_id = pdb_file.split("/")[-1].split(".")[0]
            try:
                tcr = parser.get_tcr_structure(pdb_id, pdb_file)
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
        assert len(single_chain_MHC) == 0

    def test_MHC_association(self):
        import glob

        parser = TCRParser.TCRParser()
        stcrdab_pdb_files = glob.glob(
            "/home/quast/Projects/STCRDab/Data/entries/*/structure/imgt/*.pdb"
        )
        stcrdab_pdb_files.sort()
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
        for pdb_file in stcrdab_pdb_files:
            pdb_id = pdb_file.split("/")[-1].split(".")[0]
            try:
                tcr = parser.get_tcr_structure(pdb_id, pdb_file)
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
        assert len(apo_TCRs) == 0
