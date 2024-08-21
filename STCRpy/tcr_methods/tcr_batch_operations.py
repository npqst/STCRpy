import pandas as pd

from ..tcr_processing.TCRParser import TCRParser
from ..tcr_interactions.TCRInteractionProfiler import TCRInteractionProfiler
from ..tcr_geometry.TCRGeom import TCRGeom
from ..tcr_geometry.TCRGeomFiltering import DockingGeometryFilter


class TCRBatchOperator:
    def __init__(self):
        self._tcr_parser = TCRParser()

    def _load_geometry_calculator(self):
        self._geometry_calculator = TCRGeom()

    def _load_geometry_filter(self):
        self._geometry_filter = DockingGeometryFilter()

    def tcrs_from_file_list(self, file_list):
        for file in file_list:
            tcr_id = file.split("/")[-1].split(".")[0]
            for tcr in self.tcr_parser.get_structure(tcr_id, file).get_TCRs():
                yield tcr

    def tcrs_from_file_dict(self, file_dict):
        for tcr_id, file in file_dict.items():
            for tcr in self.tcr_parser.get_structure(tcr_id, file).get_TCRs():
                yield tcr_id, tcr

    def get_TCR_pMHC_interactions(self, tcr_generator, renumber=True, save_as_csv=None):
        self.interaction_profiler = TCRInteractionProfiler()
        interaction_analysis_dict = {}
        for tcr in tcr_generator:
            tcr_id = f"{tcr.parent.parent.id}_{tcr.id}"
            if len(tcr) == 2:  # handle case where tcr is passed as (key, value)
                tcr_id, tcr = tcr
            interaction_analysis_dict[tcr_id] = (
                self.interaction_profiler.parse_tcr_pmhc_complex(tcr, renumber=renumber)
            )
        interactions_df = pd.concat(
            interaction_analysis_dict.values(),
            keys=interaction_analysis_dict.keys(),
            axis=0,
        )

        if save_as_csv is not None:
            interactions_df.to_csv(save_as_csv)

        return interactions_df


def batch_load_TCRs(tcr_files):
    if isinstance(tcr_files, dict):
        return dict(TCRBatchOperator().tcrs_from_file_dict(tcr_files))
    else:
        return list(TCRBatchOperator().tcrs_from_file_list(tcr_files))


def get_TCR_interactions(tcr_files, renumber=True, save_as_csv=None):
    batch_ops = TCRBatchOperator()
    if isinstance(tcr_files, list):
        tcr_generator = batch_ops.tcrs_from_file_list(tcr_files)
    if isinstance(tcr_files, dict):
        tcr_generator = batch_ops.tcrs_from_file_dict(tcr_files)

    return batch_ops.get_TCR_pMHC_interactions(
        tcr_generator, renumber=renumber, save_as_csv=save_as_csv
    )
