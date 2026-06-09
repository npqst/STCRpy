"""
Full-dataset integration tests. Marked @pytest.mark.slow — CI runs these on push to
main only (see .github/workflows/unittest-workflow.yml).

Each test draws a fixed random sample from the full STCRDab PDB code list and
exercises one feature area end-to-end on real-world structures.
"""
import random
import unittest
import warnings

import pytest

import stcrpy
from stcrpy.tcr_processing import abTCR, gdTCR

PDB_CODES_FILE = "./test_files/TCRParser_test_files/tcr_pdb_codes.txt"
SAMPLE_SIZE = 30

# V-gene name prefix expected for each TCR domain type
DOMAIN_GENE_PREFIXES = {
    "VA": "TRAV",
    "VB": "TRBV",
    "VG": "TRGV",
    "VD": "TRDV",
}


def _load_pdb_codes():
    with open(PDB_CODES_FILE) as f:
        return [line.strip() for line in f if line.strip()]


def _fetch_sample():
    """Return (pdb_code, tcrs) pairs for a random sample.

    Structures that raise an exception are collected into an errors dict and
    emitted as a UserWarning at the end so they appear in the test output even
    when the test ultimately passes.
    """
    results = []
    errors = {}
    for pdb_code in random.sample(_load_pdb_codes(), SAMPLE_SIZE):
        try:
            tcrs = stcrpy.fetch_TCRs(pdb_code)
            results.append((pdb_code, tcrs))
        except Exception as e:
            errors[pdb_code] = e

    if errors:
        warnings.warn(
            f"The following {len(errors)} PDB codes raised exceptions during fetch/parse "
            f"and were excluded from the sample: {errors}",
            UserWarning,
            stacklevel=2,
        )
    return results, errors


@pytest.mark.slow
class TestFullDatasetParsing(unittest.TestCase):
    """Fetch and parse a random sample; assert each structure produces valid TCR objects."""

    def test_subset_of_stcrdab(self):
        sample, errors = _fetch_sample()

        self.assertEqual(errors, {}, msg=f"Failed to fetch/parse: {errors}")
        self.assertGreater(len(sample), 0)

        for pdb_code, tcrs in sample:
            self.assertGreater(len(tcrs), 0, msg=f"{pdb_code}: no TCRs returned")
            for tcr in tcrs:
                self.assertIsInstance(
                    tcr, (abTCR, gdTCR),
                    msg=f"{pdb_code}: unexpected TCR type {type(tcr)}",
                )
                domain = tcr.get_domain_assignment()
                self.assertIsInstance(domain, dict)
                self.assertGreater(
                    len(domain), 0,
                    msg=f"{pdb_code} {tcr.id}: empty domain assignment",
                )
                # Every assigned domain must have a corresponding child chain
                for domain_name, chain_id in domain.items():
                    self.assertIn(
                        chain_id, tcr.child_dict,
                        msg=f"{pdb_code} {tcr.id}: domain {domain_name} → chain {chain_id} missing",
                    )


@pytest.mark.slow
class TestFullDatasetGeometry(unittest.TestCase):
    """
    Geometry calculations on a random sample.

    - TCRGeom is run for every TCR that has an associated MHC.
    - The docking geometry score (DockingGeometryFilter) is additionally
      calculated for TCRs with MHC class I, where the STCRDab distribution
      model applies.
    """

    def test_geometry_with_mhc(self):
        from stcrpy.tcr_geometry.TCRGeom import TCRGeom

        sample, errors = _fetch_sample()
        self.assertEqual(errors, {}, msg=f"Failed to fetch/parse: {errors}")

        geom_errors = {}
        score_errors = {}
        geom_computed = 0
        score_computed = 0

        for pdb_code, tcrs in sample:
            for tcr in tcrs:
                if not tcr.get_MHC():
                    continue

                try:
                    geom = TCRGeom(tcr)
                    self.assertIsNotNone(
                        geom.scanning_angle,
                        msg=f"{pdb_code} {tcr.id}: scanning_angle is None",
                    )
                    geom_computed += 1
                except Exception as e:
                    geom_errors[f"{pdb_code}:{tcr.id}"] = e

                # Docking geometry score for MHC class I only
                mhc_type = tcr.get_MHC()[0].get_MHC_type()
                if mhc_type == "MH1":
                    try:
                        score = tcr.score_docking_geometry()
                        self.assertIsInstance(
                            score, float,
                            msg=f"{pdb_code} {tcr.id}: score_docking_geometry did not return float",
                        )
                        self.assertFalse(
                            score != score,  # NaN check
                            msg=f"{pdb_code} {tcr.id}: score_docking_geometry returned NaN",
                        )
                        score_computed += 1
                    except Exception as e:
                        score_errors[f"{pdb_code}:{tcr.id}"] = e

        if geom_errors:
            warnings.warn(
                f"{len(geom_errors)} TCRs raised exceptions during TCRGeom calculation: {geom_errors}",
                UserWarning,
            )
        if score_errors:
            warnings.warn(
                f"{len(score_errors)} MHC class I TCRs raised exceptions during score_docking_geometry: {score_errors}",
                UserWarning,
            )

        self.assertEqual(geom_errors, {}, msg=f"TCRGeom failures: {geom_errors}")
        self.assertEqual(score_errors, {}, msg=f"score_docking_geometry failures: {score_errors}")
        self.assertGreater(geom_computed, 0, "No MHC-containing TCRs found in sample")
        self.assertGreater(score_computed, 0, "No MHC class I TCRs found in sample")


@pytest.mark.slow
class TestFullDatasetSequenceOperations(unittest.TestCase):
    """
    Germline gene assignment on a random sample.

    For each TCR chain, the assigned V-gene name must start with the IMGT
    prefix that corresponds to the TCR domain type (TRAV, TRBV, TRGV, TRDV).
    """

    def test_germline_assignments_sample(self):
        sample, errors = _fetch_sample()
        self.assertEqual(errors, {}, msg=f"Failed to fetch/parse: {errors}")

        germline_errors = {}

        for pdb_code, tcrs in sample:
            for tcr in tcrs:
                try:
                    germline_info = tcr.get_germline_assignments()
                    self.assertIsInstance(
                        germline_info, dict,
                        msg=f"{pdb_code} {tcr.id}: get_germline_assignments did not return dict",
                    )
                    self.assertGreater(
                        len(germline_info), 0,
                        msg=f"{pdb_code} {tcr.id}: empty germline assignment",
                    )

                    domain_assignment = tcr.get_domain_assignment()

                    for domain_name, chain_id in domain_assignment.items():
                        if chain_id not in germline_info:
                            continue
                        chain_germline = germline_info[chain_id]
                        if chain_germline is None or "v_gene" not in chain_germline:
                            continue

                        v_gene_entry = chain_germline["v_gene"]
                        if v_gene_entry is None:
                            continue

                        # v_gene is [(species, gene_name), score]
                        gene_name = v_gene_entry[0][1]

                        expected_prefix = DOMAIN_GENE_PREFIXES.get(domain_name)
                        if expected_prefix is not None:
                            self.assertTrue(
                                gene_name.startswith(expected_prefix),
                                msg=(
                                    f"{pdb_code} {tcr.id}: domain {domain_name} (chain {chain_id}) "
                                    f"has v_gene '{gene_name}', expected prefix '{expected_prefix}'"
                                ),
                            )
                except AssertionError:
                    raise
                except Exception as e:
                    germline_errors[f"{pdb_code}:{tcr.id}"] = e

        if germline_errors:
            warnings.warn(
                f"{len(germline_errors)} TCRs raised exceptions during germline assignment: {germline_errors}",
                UserWarning,
            )
        self.assertEqual(
            germline_errors, {},
            msg=f"Germline assignment failures: {germline_errors}",
        )
