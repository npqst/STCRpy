import unittest
import os

from ..TCRpy.tcr_processing import TCRParser


class TestTCRParser(unittest.TestCase):
    def test_get_tcr_structure_class_I(self):
        parser = TCRParser.TCRParser()

        pdb_file = 'TCRpy/test/test_files/5hyj.pdb'
        tcr = parser.get_tcr_structure('test', pdb_file)
        assert set([''.join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(['DE', 'IJ'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(['FG', 'AB'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_antigens()]) == set(['C', 'H'])

    def test_get_tcr_structure_class_II(self):
        parser = TCRParser.TCRParser()

        pdb_file = 'TCRpy/test/test_files/6r0e.cif'
        tcr = parser.get_tcr_structure('test', pdb_file)
        assert set([''.join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(['DE'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(['AB'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_antigens()]) == set(['C'])

    def test_all_stcrdab(self):
        import glob
        parser = TCRParser.TCRParser()
        stcrdab_pdb_files = glob.glob(
            '/home/quast/Projects/STCRDab/Data/entries/*/structure/imgt/*.pdb'
        )
        stcrdab_pdb_files.sort()
        badly_parsed_pdb = []
        errors = {}
        pdb_types = {}
        for pdb_file in stcrdab_pdb_files:
            pdb_id = pdb_file.split('/')[-1].split('.')[0]
            try:
                tcr = parser.get_tcr_structure(pdb_id, pdb_file)
                if len(list(tcr.get_TCRs())) == 0:
                    badly_parsed_pdb.append(pdb_id)
                else:
                    pdb_types[pdb_id] = (type(list(tcr.get_TCRs())[0]), list([str(x) for x in tcr.get_TCRs()]))                
            except Exception as e:
                errors[pdb_id] = e
        print(badly_parsed_pdb)
        print(len(badly_parsed_pdb))

    def test_docked_tcr(self):
        import glob
        parser = TCRParser.TCRParser()
        dock_pdb_files = glob.glob(
            '/home/quast/Collaborations/SamuelsLab_Weizmann/docking/results/N17.2/310569-N17-2_NRAS_rank_0/structures/it1/renumbered_complex_*.pdb'
        )
        dock_pdb_files.sort()
        badly_parsed_pdb = []
        for pdb_file in dock_pdb_files:
            pdb_id = pdb_file.split('/')[-1].split('.')[0]
            tcr = parser.get_tcr_structure(pdb_id, pdb_file)
            if len(list(tcr.get_TCRs())) == 0:
                badly_parsed_pdb.append(pdb_id)
        print(badly_parsed_pdb)
        print(len(badly_parsed_pdb))


    def test_antigen_only(self):
        parser = TCRParser.TCRParser()
        pdb_file = '/home/quast/Collaborations/SamuelsLab_Weizmann/data/antigen_xtal_structures/7mle_1_aligned.cif'
        tcr = parser.get_tcr_structure('7mle', pdb_file)
        print(tcr)


    def test_delta_beta_tcr(self):
        parser = TCRParser.TCRParser()

        pdb_file = 'TCRpy/test/test_files/DB_test_T104_rank_0_model_0_refined.pdb'
        tcr = parser.get_tcr_structure('test', pdb_file)
        assert set([''.join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(['AB'])


    def test_save(self):
        parser = TCRParser.TCRParser()

        pdb_file = 'TCRpy/test/test_files/4nhu.pdb'
        tcr = parser.get_tcr_structure('test', pdb_file)

        from ..TCRpy.tcr_processing.TCRIO import TCRIO

        io = TCRIO()
        for x in tcr.get_TCRs():
            io.save(
                x,
                save_as=f'TCRpy/test/test_files/test_{x.id}_TCR_only.pdb'
                )
        
        for x in tcr.get_TCRs():
            io.save(
                x,
                tcr_only=True,
                save_as=f'TCRpy/test/test_files/test_{x.id}.pdb'
                )
