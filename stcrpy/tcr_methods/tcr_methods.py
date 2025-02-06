import warnings

from ..tcr_processing.TCRParser import TCRParser
from .tcr_batch_operations import batch_load_TCRs


def load_TCRs(tcr_structure_files, tcr_ids=None):
    tcr_parser = TCRParser()
    if isinstance(tcr_structure_files, str):  # loading single file
        tcr_id = tcr_structure_files.split("/")[-1].split(".")[
            0
        ]  # set tcr_id to file name without extension
        if tcr_ids is not None:
            if not isinstance(tcr_ids, str):
                warnings.warn(f"TCR ID: {tcr_ids} for a single TCR should be type str.")
            tcr_id = tcr_ids

        tcr_structure = tcr_parser.get_tcr_structure(tcr_id, tcr_structure_files)
        return list(tcr_structure.get_TCRs())

    if len(tcr_structure_files) > 10:
        warnings.warn(
            "Loading more than 10 TCR structure objects into memory. Consider applying generator methods to reduce memory load."
        )

    if tcr_ids is not None:
        if len(tcr_structure_files) == len(tcr_ids):
            return batch_load_TCRs(dict(zip(tcr_ids, tcr_structure_files)))
        else:
            warnings.warn(
                f"Length of TCR IDs {len(tcr_ids)} does not match length of files {len(tcr_structure_files)}. TCR IDs reverted to default."
            )
    return batch_load_TCRs(tcr_structure_files)
