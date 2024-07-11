from typing import Dict

import numpy as np

from odusseas.utils import Spectrum


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
