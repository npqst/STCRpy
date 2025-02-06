def tcrs_to_AF3_json(tcrs, path=None, **kwargs):
    from ..tcr_formats.tcr_formats import to_AF3_json
    import json

    if isinstance(tcrs[0], str):
        from .tcr_methods import load_TCRs

        tcrs = load_TCRs(tcrs)
    else:
        from ..tcr_processing.TCR import TCR

        assert isinstance(tcrs[0], TCR)
    multiple_job_json = [to_AF3_json(tcr, save=False, **kwargs) for tcr in tcrs]
    path = path if path is not None else "stcrpy_AF3_TCRs.json"
    with open(path, "w") as f:
        json.dump(multiple_job_json, f)
    print(f"{len(tcrs)} saved as AF3 json job: {path}")
    return multiple_job_json
