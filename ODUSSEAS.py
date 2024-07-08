from enum import Enum
from typing import Annotated

import typer

import machinelearning
import spec_utils


class Reference(str, Enum):
    photometry = "photometry"
    interferometry = "interferometry"


class Regression(str, Enum):
    linear = "linear"
    ridge = "ridge"
    ridgecv = "ridgecv"
    multitasklasso = "multitasklasso"
    multitaskelasticnet = "multitaskelasticnet "


app = typer.Typer()


@app.command()
def main(
    input_spectra: str,
    reference: Annotated[
        Reference,
        typer.Option(
            help="choose the reference scale: 'photometry' for 65 stars with Teff from Casagrande08 and [Fe/H] from Neves12, or 'interferometry' for 47 stars with Teff from Khata21 and Rabus19, and [Fe/H] from Neves12",
            case_sensitive=False,
        ),
    ] = Reference.interferometry,
    regression: Annotated[
        Regression,
        typer.Option(
            help="choose the ML model. Recommended: ridge", case_sensitive=False
        ),
    ] = Regression.ridge,
    rv_cor: Annotated[bool, typer.Option()] = True,
    verbose: Annotated[bool, typer.Option()] = False,
    skip_ew_measurements: Annotated[
        bool,
        typer.Option(
            help="If this step is already done, then it can be skipped in further analysis, as it is a but slow"
        ),
    ] = False,
):
    """Run ODUSSEAS with the arguments as listed below"""
    spectra = dict()
    with open(input_spectra, "r") as f:
        for spectrum in f:
            fname, resolution = spectrum.strip().split(" ")
            spectra[fname] = resolution
    if not skip_ew_measurements:
        spec_utils.EWmeasurements(spectra, rv_cor, verbose)
    machinelearning.ML(spectra, regression.value, reference.value)


if __name__ == "__main__":
    typer.run(main)
