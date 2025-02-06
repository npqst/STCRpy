import unittest

import numpy as np

from stcrpy.tcr_geometry.TCRGeomFiltering import (
    DockingGeometryFilter,
    DockingGeometryPDF,
    GaussianMixturePDF,
    GammaPDF,
    ScanningAnglePDF,
    ZDistPDF,
    PitchAnglePDF,
)


class TestTCRGeometryFilters(unittest.TestCase):

    def test_docking_geometry_filter(self):
        dock_geom_filter = DockingGeometryFilter()

        scanning_angles = [47.5, 95.3, -121.3]
        pitch_angles = [5.3, 25.9, 13.0]
        z_dists = [28.31, 25.4, 31.0]

        # test single angles
        score = dock_geom_filter.score_docking_geometry(
            scanning_angles[0], pitch_angles[0], z_dists[0]
        )
        assert isinstance(score, float)

        # test list of angles
        scores = dock_geom_filter.score_docking_geometry(
            scanning_angles, pitch_angles, z_dists
        )
        assert len(scores) == 3
        assert all([isinstance(x, float) for x in scores])
        assert scores[2] < scores[0]

        # test implementation as self call works
        scores_self_call = dock_geom_filter(scanning_angles, pitch_angles, z_dists)
        assert all([s == scores_self_call[i] for i, s in enumerate(scores)])

        # test alternative weights
        scores = dock_geom_filter.score_docking_geometry(
            [45.0, 45.0, 90.0], pitch_angles, z_dists, weights=[1, 0, 0]
        )
        assert scores[0] == scores[1] != scores[-1]
        scores_self_call = dock_geom_filter(
            [45.0, 45.0, 90.0], pitch_angles, z_dists, weights=[1, 0, 0]
        )
        assert all([s == scores_self_call[i] for i, s in enumerate(scores)])

        scores = dock_geom_filter.score_docking_geometry(
            scanning_angles, [13.0, 13.0, 5.0], z_dists, weights=[0, 1, 0]
        )
        assert scores[0] == scores[1] != scores[-1]

        scores = dock_geom_filter.score_docking_geometry(
            scanning_angles, pitch_angles, [29.0, 29.0, 35.0], weights=[0, 0, 1]
        )
        assert scores[0] == scores[1] != scores[-1]

        scores = dock_geom_filter.score_docking_geometry(
            scanning_angles, pitch_angles, z_dists, weights=np.random.rand(3)
        )

    def test_docking_geometry_pdf(self):
        with self.assertRaises(NotImplementedError):
            DockingGeometryPDF()

    def test_gaussian_mixture_pdf(self):
        # test initialisation from parameters
        weights = np.array([0.5, 0.2, 0.3])
        means = np.array([[[10.0]], [[12]], [[-6.5]]])
        covariances = np.array([[2.3], [0.4], [8.3]])
        gmm = GaussianMixturePDF(weights=weights, means=means, covariances=covariances)
        probs = gmm.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialistaion from data
        data = 20 * np.random.rand(100) - 10.0
        gmm = GaussianMixturePDF(data=data)
        assert gmm.means.item() > -5.0 and gmm.means.item() < 5.0
        probs = gmm.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialisation from data with number of components
        data = 20 * np.random.rand(100) - 10.0
        gmm = GaussianMixturePDF(data=data, n_components=4)
        assert gmm.means.shape == (4, 1)
        assert gmm.covariances.shape == (4, 1, 1)
        probs = gmm.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialisation fails if inputs not provided:
        with self.assertRaises(ValueError):
            GaussianMixturePDF(data=None, means=None, covariances=None)

    def test_gamma_pdf(self):
        # test initialisation from parameters
        alpha = 2.43
        loc = 0.64
        scale = 3.78
        gamma = GammaPDF(alpha=alpha, scale=scale, loc=loc)
        probs = gamma.pdf([2.3, 0.01, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialistaion from data
        data = np.random.rand(100) * np.exp(-np.arange(0, 5, 5.0 / 100.0))
        gamma = GammaPDF(data=data)
        assert (
            gamma.alpha.item() > 0.0
            and gamma.loc.item() > 0.0
            and gamma.scale.item() > 0.0
        )
        probs = gamma.pdf([2.3, 0.01, 34.0])
        assert len(probs) == 3
        assert all([isinstance(p, float) for p in probs])
        assert all([p > 0.0 for p in probs])

        # test initialisation fails if inputs not provided:
        with self.assertRaises(AssertionError):
            GammaPDF(data=None, loc=None, alpha=None, scale=None)
        with self.assertRaises(AssertionError):
            GammaPDF(data=None, loc=loc, alpha=alpha, scale=None)
        with self.assertRaises(AssertionError):
            GammaPDF(data=None, loc=loc, alpha=None, scale=scale)
        with self.assertRaises(AssertionError):
            GammaPDF(data=None, loc=None, alpha=alpha, scale=scale)

    def test_scanning_angle_pdf(self):
        # test initialisation from parameters
        means = np.array([[[10.0]], [[12]], [[-6.5]]])
        covariances = np.array([[2.3], [0.4], [8.3]])
        scan_angle_pdf = ScanningAnglePDF(mean=means, variance=covariances)
        assert isinstance(scan_angle_pdf, GaussianMixturePDF)
        probs = scan_angle_pdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialistaion from data
        data = 20 * np.random.rand(100) - 10.0
        scan_angle_pdf = ScanningAnglePDF(scanning_angles=data)
        assert scan_angle_pdf.means.item() > -5.0 and scan_angle_pdf.means.item() < 5.0
        probs = scan_angle_pdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])
        # scan_angle_pdf.plot(# plotting in unittest is causing weird kernel crash, but seems to be ok when not run in test environment?
        #     data=data,
        #     save_as="./test_files/out/geometric_filters/scan_angle_test.png",
        # )

    def test_z_dist_pdf(self):
        # test initialisation from parameters
        weights = np.array([0.5, 0.2, 0.3])
        means = np.array([[[10.0]], [[12]], [[-6.5]]])
        covariances = np.array([[2.3], [0.4], [8.3]])
        zdistpdf = ZDistPDF(weights=weights, mean=means, variance=covariances)
        assert isinstance(zdistpdf, GaussianMixturePDF)
        probs = zdistpdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialistaion from data
        data = 20 * np.random.rand(100) - 10.0
        zdistpdf = ZDistPDF(z_dists=data)
        assert zdistpdf.means.item() > -5.0 and zdistpdf.means.item() < 5.0
        probs = zdistpdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialistaion from data
        data = 20 * np.random.rand(100) - 10.0
        zdistpdf = ZDistPDF(z_dists=data, n_components=3)
        assert (
            zdistpdf.means.shape == (3, 1)
            and zdistpdf.covariances.shape == (3, 1, 1)
            and zdistpdf.weights.shape == (3,)
        )
        probs = zdistpdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])
        # zdistpdf.plot(            # plotting in unittest is causing weird kernel crash, but seems to be ok when not run in test environment?
        #     data=data,
        #     save_as="./test_files/out/geometric_filters/zdist_test.png",
        # )

    def test_pitch_angle_pdf(self):
        # test initialisation from parameters
        alpha = 2.344
        loc = 46.44
        scale = 23.44
        pitch_pdf = PitchAnglePDF(alpha=alpha, loc=loc, scale=scale)
        assert isinstance(pitch_pdf, GammaPDF)
        probs = pitch_pdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([p < 1.0 for p in probs])

        # test initialistaion from data
        data = np.exp(-0.2 * np.arange(-2, 5, 7 / 100)) * np.random.rand(100) - 10.0
        pitch_pdf = PitchAnglePDF(pitch_angles=data)
        assert (
            isinstance(pitch_pdf.alpha, float)
            and isinstance(pitch_pdf.loc, float)
            and isinstance(pitch_pdf.scale, float)
        )
        probs = pitch_pdf.pdf([2.3, 45.0, 34.0])
        assert len(probs) == 3
        assert all([isinstance(p, float) for p in probs])
        # pitch_pdf.plot(  # plotting in unittest is causing weird kernel crash, but seems to be ok when not run in test environment?
        #     data=data,
        #     save_as="./test_files/out/geometric_filters/pitch_angle_test.png",
        # )

    def test_default_STCRDab_parameter_initialisation(self):
        from stcrpy.tcr_geometry import TCRGeomFiltering

        scan_angle_pdf = ScanningAnglePDF(
            mean=TCRGeomFiltering.SCANNING_ANGLE_MEAN,
            variance=TCRGeomFiltering.SCANNING_ANGLE_VARIANCE,
        )
        default_scan_angle_pdf = ScanningAnglePDF()
        probs = scan_angle_pdf.pdf([2.3, 45.0, 34.0])
        default_probs = default_scan_angle_pdf.pdf([2.3, 45.0, 34.0])
        assert tuple(probs) == tuple(default_probs)

        zdistpdf = ZDistPDF(
            weights=TCRGeomFiltering.Z_DIST_WEIGHTS,
            mean=TCRGeomFiltering.Z_DIST_MEAN,
            variance=TCRGeomFiltering.Z_DIST_VARIANCE,
        )
        default_zdistpdf = ZDistPDF()
        probs = zdistpdf.pdf([2.3, 45.0, 34.0])
        default_probs = default_zdistpdf.pdf([2.3, 45.0, 34.0])
        assert tuple(probs) == tuple(default_probs)

        pitch_pdf = PitchAnglePDF(
            alpha=TCRGeomFiltering.PITCH_ANGLE_ALPHA,
            loc=TCRGeomFiltering.PITCH_ANGLE_LOC,
            scale=TCRGeomFiltering.PITCH_ANGLE_SCALE,
        )
        default_pitch_pdf = PitchAnglePDF()
        probs = pitch_pdf.pdf([2.3, 45.0, 34.0])
        default_probs = default_pitch_pdf.pdf([2.3, 45.0, 34.0])
        assert tuple(probs) == tuple(default_probs)
