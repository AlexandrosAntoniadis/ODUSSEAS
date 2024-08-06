import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.interpolate import interp1d

from odusseas.build_pew import area_between, cut_data

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


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))


def find_rv(
    wavelength,
    flux,
    mask_lines=[
        5371.50,
        5394.64,
        5397.10,
        5405.80,
        5409.80,
        5429.70,
        5432.55,
        5434.52,
        5446.90,
        5535.50,
        6090.20,
        6102.72,
        6119.52,
        6122.22,
    ],
    delta_rv=200,
):
    wave_mask = np.arange(mask_lines[0] - 5.0, mask_lines[-1] + 5.0, 0.01)
    flux_mask = 1.0 + wave_mask * 0.0
    for line in mask_lines:
        flux_mask -= gaussian(wave_mask, line, 0.1)

    ind_wi = np.where(wavelength < wave_mask[0])[0][-1]
    ind_wf = np.where(wavelength > wave_mask[-1])[0][0]

    wave_t = wavelength[ind_wi:ind_wf]
    flux_t = flux[ind_wi:ind_wf]

    rv, cc = pyasl.crosscorrRV(
        wave_t, flux_t, wave_mask, flux_mask, -delta_rv, delta_rv, 0.1, skipedge=900
    )
    maxind = np.argmax(cc)

    return rv[maxind]


class ReferenceEnum(str, Enum):
    photometry = "photometry"
    interferometry = "interferometry"


class RegressionEnum(str, Enum):
    linear = "linear"
    ridge = "ridge"
    ridgecv = "ridgecv"
    multitasklasso = "multitasklasso"
    multitaskelasticnet = "multitaskelasticnet "


@dataclass
class Spectrum:
    fname: str
    resolution: int

    def __post_init__(self):
        self.flux = fits.getdata(self.fname)
        hdr = fits.getheader(self.fname)
        w0, dw, N = hdr["CRVAL1"], hdr["CDELT1"], hdr["NAXIS1"]
        cdelt1 = 0.010
        self.wavelength = w0 + dw * np.arange(N)
        if round(dw, 4) != 0.010:
            f2 = interp1d(self.wavelength, self.flux, kind="linear")
            self.wavelength = np.arange(self.wavelength[0], self.wavelength[-1], cdelt1)
            self.flux = f2(self.wavelength)

        rv = find_rv(self.wavelength, self.flux)
        wave_rv = self.wavelength / (1 + rv / 3.0e5) if abs(rv) > 0 else self.wavelength
        f2 = interp1d(wave_rv, self.flux, kind="linear")
        self.wavelength = np.arange(wave_rv[0], wave_rv[-1] - cdelt1, cdelt1)
        self.flux = f2(self.wavelength)

    @property
    def name(self) -> str:
        return os.path.basename(self.fname).replace(".fits", "")

    def pseudo_EW(
        self, wavelength_range: Tuple[float, float], dw: float = 0.4, verbose=False
    ):
        """
        Calculate the pseudo equivalent width of the spectrum in a given range.
        """
        w1, w2 = wavelength_range
        wavelength, flux = cut_data(self.wavelength, self.flux, w1, w2)
        # Find central wavelength
        idx_fmin = np.argmin(flux)
        critical_wl = wavelength[idx_fmin]
        # Work on left side
        wavelength_left, flux_left = cut_data(
            wavelength, flux, critical_wl - dw, critical_wl
        )
        idx_left = np.argmax(flux_left)
        wl_max_flux_left = wavelength_left[idx_left]
        fl_max_flux_left = flux_left[idx_left]
        # Work on right side
        wavelength_right, flux_right = cut_data(
            wavelength, flux, critical_wl, critical_wl + dw
        )
        idx_right = np.argmax(flux_right)
        wl_max_flux_right = wavelength_right[idx_right]
        fl_max_flux_right = flux_right[idx_right]
        # Set the area
        x1, x2 = wl_max_flux_left, wl_max_flux_right
        y1, y2 = fl_max_flux_left, fl_max_flux_right
        g = np.polyfit([x1, x2], [y1, y2], deg=1)
        N = len(flux_left[idx_left::]) + len(flux_right[0:idx_right])
        x = np.linspace(x1, x2, N)
        idx = (x1 <= wavelength) & (wavelength <= x2)
        f = flux[idx]
        g = np.poly1d(g)(x)
        area = area_between(f, g, dx=x[1] - x[0]) * 1000
        if verbose:
            print(r"Area of line: {:.2f}mÅ".format(area))
        return area

    def save_pseudo_EW(self, output: np.ndarray, wcentral: np.ndarray) -> None:
        os.makedirs("results", exist_ok=True)
        fname = os.path.join(
            "results", os.path.basename(self.fname).replace(".fits", ".dat")
        )

        table = np.zeros((len(output), 2))
        table[:, 0] = wcentral
        table[:, 1] = output
        np.savetxt(fname, table, delimiter=",", fmt="%s")

    def get_wavelength_ranges(self, filepath: str) -> np.ndarray:
        wavelength_ranges = np.loadtxt(filepath, skiprows=2)
        w1_max, w2_min = np.max(wavelength_ranges[:, 0]), np.min(
            wavelength_ranges[:, 1]
        )
        # Adjust wavelength ranges to only include lines in the range of the spectrum
        starwavelength_min, starwavelength_max = self.wavelength[0], self.wavelength[-1]
        wavelength_ranges = wavelength_ranges
        if starwavelength_min > w2_min:
            mask = np.where(wavelength_ranges[:, 1] > starwavelength_min)
            wavelength_ranges = wavelength_ranges[mask[0]]
        if starwavelength_max < w1_max:
            mask = np.where(wavelength_ranges[:, 0] < starwavelength_max)
            wavelength_ranges = wavelength_ranges[mask[0]]
        return wavelength_ranges


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
        return (
            Path(__file__).parent
            / f"reference_data/res{self.resolution}_{self.reference}RefEWPar.csv"
        )

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
