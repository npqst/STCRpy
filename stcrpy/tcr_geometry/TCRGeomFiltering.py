import warnings
from scipy.stats import norm, gamma
from sklearn.mixture import GaussianMixture
import numpy as np


# DEFAULT PARAMETERS FROM FIT TO STCRDAB (SAMPLED JULY 2024)

# Scanning angle is fit with unimodel gaussian
SCANNING_ANGLE_MEAN = np.array([[67.92200037]])
SCANNING_ANGLE_VARIANCE = np.array([[202.11710875]])

# Pitch angle is fit with gamma distribution
PITCH_ANGLE_ALPHA = 1.1993765213484138
PITCH_ANGLE_LOC = 0.036500005492399054
PITCH_ANGLE_SCALE = 5.781270719307026

# Z distance of TCR CoM is fit with bimodal gaussian mixture
Z_DIST_WEIGHTS = np.array([0.43420871, 0.56579129])
Z_DIST_VARIANCE = np.array([[[0.46845205]], [[0.90790909]]])
Z_DIST_MEAN = np.array([[27.04673494], [28.62071785]])


class DockingGeometryFilter:

    def __init__(self):

        self.scanning_angle_pdf = ScanningAnglePDF(
            mean=SCANNING_ANGLE_MEAN, variance=SCANNING_ANGLE_VARIANCE
        )
        self.pitch_angle_pdf = PitchAnglePDF(
            alpha=PITCH_ANGLE_ALPHA, loc=PITCH_ANGLE_LOC, scale=PITCH_ANGLE_SCALE
        )
        self.z_dist_pdf = ZDistPDF(
            weights=Z_DIST_WEIGHTS, mean=Z_DIST_MEAN, variance=Z_DIST_VARIANCE
        )

    def __call__(self, *args, **kwargs):
        return self.score_docking_geometry(*args, **kwargs)

    def score_docking_geometry(
        self, scanning_angles, pitch_angles, z_dists, weights=np.array([1, 1, 1])
    ):
        scanning_angle_log_probs = self.scanning_angle_pdf.log_probs(scanning_angles)
        pitch_angle_log_probs = self.pitch_angle_pdf.log_probs(pitch_angles)
        z_dist_log_probs = self.z_dist_pdf.log_probs(z_dists)

        scores = (
            weights[0] * scanning_angle_log_probs
            + weights[1] * pitch_angle_log_probs
            + weights[2] * z_dist_log_probs
        )
        if len(scores) == 1:
            return scores.item()
        return scores

    def set_scanning_angle_pdf(self, new_scanning_angle_pdf):
        assert isinstance(
            new_scanning_angle_pdf, DockingGeometryPDF
        ), "Bespoke scanning angle PDF must inherit from DockingGeometryPDF"
        try:
            log_probs = new_scanning_angle_pdf.log_probs(47.2)
            assert (
                isinstance(log_probs, float) and log_probs < 0.0
            ), "Log probs returned by new PDF not correctly configured. Scanning angle PDF has not been set, using default."
            self.scanning_angle_pdf = new_scanning_angle_pdf
        except Exception as e:
            warnings.warn(
                f"Bespoke scanning angle PDF failed test with error {str(e)}.  Scanning angle PDF has not been set, using default."
            )
            return

    def set_zdist_angle_pdf(self, new_zdist_pdf):
        assert isinstance(
            new_zdist_pdf, DockingGeometryPDF
        ), "Bespoke Z distance PDF must inherit from DockingGeometryPDF"
        try:
            log_probs = new_zdist_pdf.log_probs(29.3)
            assert (
                isinstance(log_probs, float) and log_probs < 0.0
            ), "Log probs returned by new PDF not correctly configured. Z distance PDF has not been set, using default."
            self.zdist_pdf = new_zdist_pdf
        except Exception as e:
            warnings.warn(
                f"Bespoke Z distance PDF failed test with error {str(e)}.  Z distance PDF has not been set, using default."
            )
            return

    def set_pitch_angle_pdf(self, new_pitch_angle_pdf):
        assert isinstance(
            new_pitch_angle_pdf, DockingGeometryPDF
        ), "Bespoke pitch angle PDF must inherit from DockingGeometryPDF"
        try:
            log_probs = new_pitch_angle_pdf.log_probs(12.3)
            assert (
                isinstance(log_probs, float) and log_probs < 0.0
            ), "Log probs returned by new PDF not correctly configured. Pitch angle PDF has not been set, using default."
            self.pitch_angle_pdf = new_pitch_angle_pdf
        except Exception as e:
            warnings.warn(
                f"Bespoke pitch angle PDF failed test with error {str(e)}.  Pitch angle PDF has not been set, using default."
            )
            return


class DockingGeometryPDF:
    def __init__(self, data=None):
        if data is None:
            self._pdf_from_parameters()
        else:
            if any([param is not None for param in self.params]):
                warnings.warn(
                    """Some parameters have been provided alongside data. Using data to fit new parameters."""
                )
            self.fit(data)

    def _pdf_from_parameters(self):
        raise NotImplementedError

    def fit(self, data):
        raise NotImplementedError

    def pdf(self, x):
        pass

    def log_probs(self, x):
        return np.log(self.pdf(x))

    def plot(self, data=None, plot_min=0, plot_max=150, save_as=None):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        if data is not None:
            plt.hist(data, bins=50, density=True)
            data_range = max(data) - min(data)
            plot_min = min(data) - (0.1 * data_range)
            plot_max = max(data) + (0.1 * data_range)
        x = np.linspace(plot_min, plot_max, 100)
        plt.plot(x, self.pdf(x))
        if save_as:
            plt.savefig(save_as)
            print(f"PDF distribution plot saved to: {save_as}")
        else:
            plt.show()


class GaussianMixturePDF(DockingGeometryPDF):
    def __init__(
        self,
        data=None,
        n_components=None,
        weights=None,
        means=None,
        covariances=None,
    ):
        if (
            means is not None and covariances is not None
        ):  # check if parameters are provided
            if data is not None:
                warnings.warn(
                    "Data and parameters provided for GMM. Default to fitting GMM from data"
                )
            if data is None and n_components is not None:
                warnings.warn(
                    """Nr of components for Gaussian Mixture Model is only set when fitting to data.
                    When initialising the PDF from parameters the number of components is inferred from the parameters. """
                )

            if weights is None:  # set weights of GMM if not provided
                weights = np.array([1.0 / len(means) for _ in range(len(means))])
            self.n_components = len(means)
        self.params = (weights, means, covariances)

        if (
            data is not None
        ):  # if initialising from data number or components is defined, defaults to 1.
            self.n_components = n_components if n_components is not None else 1

        if data is None and (means is None or covariances is None):
            raise ValueError(
                "GMM needs to be initialised form parameters or data. Check input arguments."
            )

        self.model = GaussianMixture(n_components=self.n_components)
        super().__init__(data)

    def _pdf_from_parameters(self):
        assert all(
            [param is not None for param in self.params]
        ), """Weights, means and covariances must be provided to define PDF from parameters instead of fitting to data.
        Please provide parameters or data to fit."""
        self.weights, self.means, self.covariances = self.params

    def fit(self, data):
        if len(data.shape) == 1:
            data = data.reshape(-1, 1)
        self.model.fit(data)
        self.weights, self.means, self.covariances = (
            self.model.weights_,
            self.model.means_,
            self.model.covariances_,
        )

    def pdf(self, x):
        pdf = np.asarray(
            [
                self.weights[i]
                * norm.pdf(x, self.means[i, 0], np.sqrt(self.covariances[i, 0]))
                for i in range(self.n_components)
            ]
        )
        pdf = np.sum(pdf, axis=0)
        return pdf


class GammaPDF(DockingGeometryPDF):
    def __init__(self, data=None, alpha=None, loc=None, scale=None):
        self.params = (alpha, loc, scale)
        self.model = gamma  # gamma from scipy.stats
        super().__init__(data)

    def _pdf_from_parameters(self):
        assert all(
            [param is not None for param in self.params]
        ), """Alpha, location and scale must be provided to define PDF from parameters instead of fitting to data.
        Please provide parameters or data to fit."""
        self.alpha, self.loc, self.scale = self.params

    def fit(self, data):
        self.alpha, self.loc, self.scale = self.model.fit(data)

    def pdf(self, x):
        pdf = self.model.pdf(x, a=self.alpha, loc=self.loc, scale=self.scale)
        return pdf


class ScanningAnglePDF(GaussianMixturePDF):
    def __init__(
        self,
        scanning_angles=None,
        mean=SCANNING_ANGLE_MEAN,
        variance=SCANNING_ANGLE_VARIANCE,
    ):
        super().__init__(data=scanning_angles, means=mean, covariances=variance)


class ZDistPDF(GaussianMixturePDF):
    def __init__(
        self,
        z_dists=None,
        n_components=None,
        weights=Z_DIST_WEIGHTS,
        mean=Z_DIST_MEAN,
        variance=Z_DIST_VARIANCE,
    ):
        super().__init__(
            z_dists,
            n_components=n_components,
            weights=weights,
            means=mean,
            covariances=variance,
        )


class PitchAnglePDF(GammaPDF):
    def __init__(
        self,
        pitch_angles=None,
        alpha=PITCH_ANGLE_ALPHA,
        loc=PITCH_ANGLE_LOC,
        scale=PITCH_ANGLE_SCALE,
    ):
        super().__init__(data=pitch_angles, alpha=alpha, loc=loc, scale=scale)
