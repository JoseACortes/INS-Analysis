import Spectrum

from calibration import evaluate

def apply(filename: str, calibration: dict, **kwargs) -> list[float]:
    # load the file
    bins, vals = Spectrum.read(filename, **kwargs)

    # apply the calibration
    return evaluate(bins, vals, calibration)