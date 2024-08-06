from pathlib import Path
from typing import Dict

import numpy as np

from odusseas.utils import Spectrum


def EWmeasurements(spectra: Dict[str, int], verbose: bool = False) -> None:
    dw = 0.4
    for fname, resolution in spectra.items():
        spectrum = Spectrum(fname, resolution)
        print(f"Calculating EW for {spectrum.name} ...")
        wavelength_ranges = spectrum.get_wavelength_ranges(
            Path(__file__).parent / "lines.rdb"
        )

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
