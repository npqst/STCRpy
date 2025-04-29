import warnings
import os
import pandas as pd

from ..tcr_processing.TCRParser import TCRParser
from ..tcr_interactions.TCRInteractionProfiler import TCRInteractionProfiler
from ..tcr_geometry.TCRGeom import TCRGeom
from ..tcr_geometry.TCRGeomFiltering import DockingGeometryFilter
from ..tcr_formats.tcr_formats import get_sequences


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
            try:
                for tcr in self._tcr_parser.get_tcr_structure(tcr_id, file).get_TCRs():
                    yield tcr
            except Exception as e:
                warnings.warn(f"Loading {file} failed with error {str(e)}")
                yield None

    def tcrs_from_file_dict(self, file_dict):
        for tcr_id, file in file_dict.items():
            try:
                for tcr in self._tcr_parser.get_tcr_structure(tcr_id, file).get_TCRs():
                    yield tcr_id, tcr
            except Exception as e:
                warnings.warn(f"Loading {tcr_id}: {file} failed with error {str(e)}")
                yield None

    def get_TCR_pMHC_interactions(self, tcr_generator, renumber=True, save_as_csv=None):
        interaction_analysis_dict = {}
        for tcr in tcr_generator:
            if tcr is None:  # handles case where file could not be parsed in generator
                continue
            tcr_id = f"{tcr.parent.parent.id}_{tcr.id}"
            if isinstance(
                tcr, tuple
            ):  # handle case where tcr is passed as (key, value)
                tcr_id, tcr = tcr
            try:
                interaction_analysis_dict[tcr_id] = tcr.profile_peptide_interactions()
            except Exception as e:
                warnings.warn(
                    f"Interactions profile failed for {tcr} with error {str(e)}"
                )
        interactions_df = pd.concat(
            interaction_analysis_dict.values(),
            keys=interaction_analysis_dict.keys(),
            axis=0,
        )

        if save_as_csv is not None:
            interactions_df.to_csv(save_as_csv)

        return interactions_df

    def get_TCR_geometry(self, tcr_generator, mode="rudolph", save_as_csv=None):
        geometries_dict = {}
        for tcr in tcr_generator:
            if tcr is None:  # handles case where file could not be parsed in generator
                continue

            if isinstance(
                tcr, tuple
            ):  # handle case where tcr is passed as (key, value)
                tcr_id, tcr = tcr
            else:
                tcr_id = f"{tcr.parent.parent.id}_{tcr.id}"
            try:
                geometries_dict[tcr_id] = tcr.calculate_docking_geometry(
                    mode=mode, as_df=True
                )
            except Exception as e:
                warnings.warn(
                    f"Geometry calculation failed for {tcr} with error {str(e)}"
                )
        geometries_df = pd.concat(geometries_dict).droplevel(1)

        if save_as_csv is not None:
            geometries_df.to_csv(save_as_csv)

        return geometries_df

    def get_germlines_and_alleles(self, tcr_generator, save_as_csv=None):
        germlines_and_alleles_dict = {}
        for tcr in tcr_generator:
            if tcr is None:  # handles case where file could not be parsed in generator
                continue
            tcr_id = f"{tcr.parent.parent.id}_{tcr.id}"
            if isinstance(
                tcr, tuple
            ):  # handle case where tcr is passed as (key, value)
                tcr_id, tcr = tcr
            germlines_and_alleles_dict[tcr_id] = tcr.get_germlines_and_alleles()

        germlines_and_alleles_df = pd.DataFrame(germlines_and_alleles_dict).T

        if save_as_csv is not None:
            germlines_and_alleles_df.to_csv(save_as_csv)

        return germlines_and_alleles_df

    def full_analysis(self, tcr_generator, geometry_mode="rudolph", save_dir=None):
        from tqdm import tqdm

        germlines_and_alleles_dict = {}
        geometries_dict = {}
        interaction_analysis_dict = {}

        for tcr in tqdm(tcr_generator):
            if tcr is None:  # handles case where file could not be parsed in generator
                continue
            if isinstance(
                tcr, tuple
            ):  # handle case where tcr is passed as (key, value)
                tcr_id, tcr = tcr
            else:
                tcr_id = f"{tcr.parent.parent.id}_{tcr.id}"
            try:
                germlines_and_alleles_dict[tcr_id] = tcr.get_germlines_and_alleles()
            except Exception as e:
                warnings.warn(
                    f"Germline and allele retrieval failed for {tcr} with error {str(e)}"
                )
            try:
                geometries_dict[tcr_id] = tcr.calculate_docking_geometry(
                    mode=geometry_mode, as_df=True
                )
            except Exception as e:
                warnings.warn(
                    f"Geometry calculation failed for {tcr} with error {str(e)}"
                )
            try:
                interaction_analysis_dict[tcr_id] = tcr.profile_peptide_interactions()
            except Exception as e:
                warnings.warn(
                    f"Interaction profiling failed for {tcr} with error {str(e)}"
                )
        germlines_and_alleles_df = pd.DataFrame(germlines_and_alleles_dict).T

        geometries_df = pd.concat(geometries_dict).droplevel(1)

        interactions_df = pd.concat(
            interaction_analysis_dict.values(),
            keys=interaction_analysis_dict.keys(),
            axis=0,
        )

        if save_dir is not None:
            geometries_df.to_csv(os.path.join(save_dir, "geometries.csv"))
            germlines_and_alleles_df.to_csv(
                os.path.join(save_dir, "germlines_and_alleles.csv")
            )
            interactions_df.to_csv(os.path.join(save_dir, "interactions.csv"))

        return germlines_and_alleles_df, geometries_df, interactions_df


def batch_load_TCRs(tcr_files):
    if isinstance(tcr_files, dict):
        return dict(TCRBatchOperator().tcrs_from_file_dict(tcr_files))
    else:
        return list(TCRBatchOperator().tcrs_from_file_list(tcr_files))


def batch_yield_TCRs(tcr_files):
    if isinstance(tcr_files, dict):
        return TCRBatchOperator().tcrs_from_file_dict(tcr_files)
    else:
        return TCRBatchOperator().tcrs_from_file_list(tcr_files)


def get_TCR_interactions(tcr_files, renumber=True, save_as_csv=None):
    batch_ops = TCRBatchOperator()
    if isinstance(tcr_files, list):
        tcr_generator = batch_ops.tcrs_from_file_list(tcr_files)
    if isinstance(tcr_files, dict):
        tcr_generator = batch_ops.tcrs_from_file_dict(tcr_files)

    return batch_ops.get_TCR_pMHC_interactions(
        tcr_generator, renumber=renumber, save_as_csv=save_as_csv
    )


def get_TCR_geometry(tcr_files, mode="rudolph", save_as_csv=None):
    batch_ops = TCRBatchOperator()
    if isinstance(tcr_files, list):
        tcr_generator = batch_ops.tcrs_from_file_list(tcr_files)
    if isinstance(tcr_files, dict):
        tcr_generator = batch_ops.tcrs_from_file_dict(tcr_files)

    return batch_ops.get_TCR_geometry(tcr_generator, mode=mode, save_as_csv=save_as_csv)


def get_germlines_and_alleles(tcr_files, save_as_csv=None):
    batch_ops = TCRBatchOperator()
    if isinstance(tcr_files, list):
        tcr_generator = batch_ops.tcrs_from_file_list(tcr_files)
    if isinstance(tcr_files, dict):
        tcr_generator = batch_ops.tcrs_from_file_dict(tcr_files)

    return batch_ops.get_germlines_and_alleles(tcr_generator, save_as_csv=save_as_csv)


def analyse_tcrs(tcr_files, save_dir=None):
    batch_ops = TCRBatchOperator()
    if isinstance(tcr_files, list):
        tcr_generator = batch_ops.tcrs_from_file_list(tcr_files)
    if isinstance(tcr_files, dict):
        tcr_generator = batch_ops.tcrs_from_file_dict(tcr_files)

    return batch_ops.full_analysis(tcr_generator, save_dir=save_dir)
