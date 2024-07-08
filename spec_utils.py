import os
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.interpolate import interp1d

from build_pew import area_between, cut_data


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


def rv_correction(fname: str, cdelt1=0.010):
    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr["CRVAL1"], hdr["CDELT1"], hdr["NAXIS1"]
    wave = w0 + dw * np.arange(N)

    rv = find_rv(wave, flux)

    print("RV:", rv)

    wave_rv = wave / (1 + rv / 3.0e5) if abs(rv) > 0 else wave

    f2 = interp1d(wave_rv, flux, kind="linear")
    wave_int = np.arange(wave_rv[0], wave_rv[-1] - cdelt1, cdelt1)
    flux_int = f2(wave_int)
    flux = flux_int
    hdr["CRVAL1"] = wave_rv[0]
    hdr["CDELT1"] = cdelt1
    fits.writeto(fname, flux, hdr, overwrite=True)


@dataclass
class Spectrum:
    fname: str
    resolution: int
    do_rv_cor: bool = False

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

        if self.do_rv_cor:
            rv = find_rv(self.wavelength, self.flux)
            wave_rv = (
                self.wavelength / (1 + rv / 3.0e5) if abs(rv) > 0 else self.wavelength
            )
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


def EWmeasurements(
    spectra: Dict[str, int], do_rv_cor: bool, verbose: bool = False
) -> None:
    dw = 0.4
    for fname, resolution in spectra.items():
        spectrum = Spectrum(fname, resolution, do_rv_cor=do_rv_cor)
        print(f"Calculating EW for {spectrum.name} ...")
        wavelength_ranges = spectrum.get_wavelength_ranges("lines.rdb")

        # Calculate and save the pseudo EWs
        output = np.zeros(len(wavelength_ranges))
        for j, wavelength_range in enumerate(wavelength_ranges):
            output[j] = spectrum.pseudo_EW(
                wavelength_range=wavelength_range,
                dw=dw,
                verbose=verbose,
            )

        wcentral = np.array([(w[0] + w[1]) / 2 for w in wavelength_ranges])
        spectrum.save_pseudo_EW(output, wcentral)
