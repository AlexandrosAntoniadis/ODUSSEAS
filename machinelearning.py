import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.linear_model._base import LinearModel
from sklearn.metrics import explained_variance_score, r2_score
from sklearn.model_selection import train_test_split

from spec_utils import Spectrum

error_matrix = {
    "photometry": {
        "115000": (0.10, 65),
        "110000": (0.10, 68),
        "94600": (0.12, 77),
        "75000": (0.13, 78),
        "48000": (0.13, 80),
    },
    "interferometry": {
        "115000": (0.11, 90),
        "110000": (0.11, 92),
        "94600": (0.12, 95),
        "75000": (0.13, 97),
        "48000": (0.13, 99),
    },
}


def make_plot(x, y, extra_coords, title, xlabel, ylabel, figsize, fout):
    set_res = 15
    plt.clf()
    _, ax = plt.subplots(figsize=figsize)
    if title:
        ax.set_title(title, fontsize=set_res * 1.5)
    ax.set_xlabel(xlabel, fontsize=set_res * 1.5)
    ax.set_ylabel(ylabel, fontsize=set_res * 1.5)
    ax.plot(*extra_coords, "--b", lw=2)
    ax.plot(x, y, "ko")
    ax.tick_params(axis="both", labelsize=set_res * 1.5)
    ax.spines["right"].set_visible(True)
    ax.spines["top"].set_visible(True)
    plt.savefig(fout, bbox_inches="tight")
    plt.close()


def mean_abso_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred)), axis=0).values


def get_regression_model(name: str) -> LinearModel:
    match name:
        case "linear":
            return linear_model.LinearRegression()
        case "ridge":
            return linear_model.Ridge(alpha=1)
        case "ridgecv":
            return linear_model.RidgeCV(alphas=[0.1, 1.0, 10.0], cv=5, gcv_mode=None)
        case "multitasklasso":
            return linear_model.MultiTaskLasso(alpha=1)
        case "multitaskelasticnet":
            return linear_model.MultiTaskElasticNet(alpha=1)
    raise ValueError(f"Regression '{name}' is not valid")


@dataclass
class Reference:
    resolution: int
    reference: str

    def __post_init__(self) -> None:
        self.data = []
        self.header = None
        with open(self.fname) as f:
            for i, line in enumerate(f):
                if i > 0:
                    self.data.append(line.strip().split(","))
                else:
                    self.header = line.strip().split(",")
        if not self.header:
            raise ValueError("Header not found")
        self.data = np.array(self.data)

    @property
    def fname(self) -> str:
        return f"reference_data/res{self.resolution}_{self.reference}RefEWPar.csv"

    def subset_with_wavelengths(
        self, wavelengths: np.ndarray
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        res = dict()
        for row, col in zip(self.data.T, self.header):
            if col in ("names", "FeH", "Teff"):
                res[col] = row
                continue
            wavelength = np.round(np.float64(col), 3)
            if not wavelength in wavelengths:
                continue
            res[wavelength] = row
        df = pd.DataFrame(res)
        df_x = df.drop(["names", "FeH", "Teff"], axis=1).astype(float)
        df_y = df[["names", "FeH", "Teff"]]
        df_y.loc[:, "FeH"] = df_y["FeH"].astype(float)
        df_y.loc[:, "Teff"] = df_y["Teff"].astype(float)
        return df_x, df_y


@dataclass
class Measurement:
    fname: str
    resolution: int

    def __post_init__(self) -> None:
        self.spectrum = Spectrum(self.fname, self.resolution)

        self.newdf = pd.read_csv(
            os.path.join("results", f"{self.spectrum.name}.dat"),
            header=None,
            names=["wavelength", "EW"],
        )
        self.newdf.dropna(axis=1, inplace=True)
        self.newdf[self.newdf.EW > 0.00001]
        self.newdf.wavelength = self.newdf.wavelength.round(3)


@dataclass
class Result:
    resolution: str
    reference: str
    spectrum: str
    params_sub_result: Dict[str, List[float]]

    def __str__(self) -> str:
        d = self.params_sub_result
        wide_error_FeH = round(
            (
                (np.std(d["fehs"])) ** 2
                + (error_matrix[self.reference][self.resolution][0]) ** 2
            )
            ** (1 / 2),
            3,
        )
        wide_error_Teff = round(
            (
                (np.std(d["teffs"])) ** 2
                + (error_matrix[self.reference][self.resolution][1]) ** 2
            )
            ** (1 / 2),
            3,
        )
        feh_mean = round(np.mean(d["fehs"]), 3)
        feh_std = round(np.std(d["fehs"]), 3)
        feh_mae = round(np.mean(d["feh_maes"]), 3)
        teff_mean = int(np.mean(d["teffs"]))
        teff_std = int(np.std(d["teffs"]))
        teff_mae = int(np.mean(d["teff_maes"]))
        r2_mean = round(np.mean(d["r2scores"]), 3)
        var_mean = round(np.mean(d["variances"]), 3)

        return f"{self.spectrum} {feh_mean} {feh_std} {feh_mae} {wide_error_FeH} {teff_mean} {teff_std} {teff_mae} {wide_error_Teff} {r2_mean} {var_mean}\n"


def ML(spectra: Dict[str, int], regression: str, reference: str) -> None:
    reg = get_regression_model(regression)
    plots = "Model_Prediction_Plots"
    os.makedirs(plots, exist_ok=True)

    results = []
    for fname, resolution in spectra.items():
        # This is the EW we already measured
        measurement = Measurement(fname, resolution)
        spectrum = measurement.spectrum
        print(f"Getting atmospheric parameters for {spectrum.name} ...")
        newdf = measurement.newdf
        # Get the reference data
        ref = Reference(resolution, reference)
        df_x, df_y = ref.subset_with_wavelengths(newdf.wavelength.values)

        result = Result(
            resolution=str(resolution),
            reference=reference,
            spectrum=spectrum.name,
            params_sub_result={
                "fehs": [],
                "teffs": [],
                "feh_maes": [],
                "teff_maes": [],
                "variances": [],
                "r2scores": [],
            },
        )
        for _ in range(100):
            x_train, x_test, y_train, y_test = train_test_split(
                df_x, df_y, test_size=0.20
            )

            labels_train = y_train[["names"]]
            labels_test = y_test[["names"]]

            y_train.drop(["names"], axis=1, inplace=True)
            y_test.drop(["names"], axis=1, inplace=True)

            reg.fit(x_train.values, y_train)

            y_pred_test = reg.predict(x_test.values)

            variancescore = float(explained_variance_score(y_test, y_pred_test))
            r2score = float(r2_score(y_test, y_pred_test))

            pred_newdf = reg.predict(
                newdf.set_index("wavelength").T
            )  # applying the saved model

            mae_test = mean_abso_error(y_test[:], y_pred_test[:])

            train_givenvalues = pd.concat([labels_train, y_train], axis=1)
            train_givenvalues = train_givenvalues.reset_index(drop=True)

            test_givenvalues = pd.concat([labels_test, y_test], axis=1)
            test_givenvalues = test_givenvalues.reset_index(drop=True)
            new_labeltest = labels_test.reset_index(drop=True)
            test_predvalues = pd.concat(
                [new_labeltest, pd.DataFrame(y_pred_test)], axis=1
            )

            result.params_sub_result["fehs"].append(pred_newdf[0][0])
            result.params_sub_result["teffs"].append(pred_newdf[0][1])
            result.params_sub_result["feh_maes"].append(mae_test[0])
            result.params_sub_result["teff_maes"].append(mae_test[1])
            result.params_sub_result["r2scores"].append(r2score)
            result.params_sub_result["variances"].append(variancescore)

        results.append(result)

        # plots of the FeH test
        set_res = 15
        make_plot(
            x=test_predvalues.values[:, 1],
            y=test_givenvalues.values[:, 1],
            extra_coords=[(-0.8, 0.4), (-0.8, 0.4)],
            title="[Fe/H] model testing",
            xlabel="M.L. [Fe/H] [dex]",
            ylabel="Ref. [Fe/H] [dex]",
            figsize=(set_res * 0.8, set_res * 0.5),
            fout=f"{plots}/{spectrum.name}_FeH_test_comparison.pdf",
        )

        make_plot(
            x=test_predvalues.values[:, 1],
            y=test_predvalues.values[:, 1] - test_givenvalues.values[:, 1],
            extra_coords=[(-0.8, 0.4), (0, 0)],
            title=None,
            xlabel="M.L. [Fe/H] [dex]",
            ylabel=r"$\Delta$[Fe/H] [dex]",
            figsize=(set_res * 0.8, set_res * 0.2),
            fout=f"{plots}/{spectrum.name}_Diff_FeH_test_comparison.pdf",
        )

        make_plot(
            x=test_predvalues.values[:, 2],
            y=test_givenvalues.values[:, 2],
            extra_coords=[(2700, 4000), (2700, 4000)],
            title=r"T$_{\mathrm{eff}}$ model testing",
            xlabel=r"M.L. T$_{\mathrm{eff}}$ [K]",
            ylabel=r"Ref. T$_{\mathrm{eff}}$ [K]",
            figsize=(set_res * 0.8, set_res * 0.5),
            fout=f"{plots}/{spectrum.name}_Teff_test_comparison.pdf",
        )

        make_plot(
            x=test_predvalues.values[:, 2],
            y=test_predvalues.values[:, 2] - test_givenvalues.values[:, 2],
            extra_coords=[(2700, 4000), (0, 0)],
            title=None,
            xlabel=r"M.L. T$_{\mathrm{eff}}$ [K]",
            ylabel=r"$\Delta$T$_{\mathrm{eff}}$ [K]",
            figsize=(set_res * 0.8, set_res * 0.2),
            fout=f"{plots}/{spectrum.name}_Diff_Teff_test_comparison.pdf",
        )

    # Save all the results
    with open("Parameter_Results.dat", "w") as fh:
        fh.write(
            "newstars [Fe/H] STD_FeH MAE_[Fe/H] Wide_Error_[Fe/H] Teff STD_Teff MAE_Teff Wide_Error_Teff R2_score EV_score\n"
        )
        for result in results:
            fh.write(str(result))
    print("Results saved to Parameter_Results.dat")
