import Spectrum

from calibration import applyCalibrationAreas
from areas import calcPeakAreas

def apply(
        bins,
        vals,  
        calibration: dict, **kwargs) -> list[float]:
    
    # calculate the areas
    areas = calcPeakAreas(bins, vals, **kwargs)
    fits = None
    
    if isinstance(areas, tuple):
        areas, fits = areas

    # apply the calibration
    if fits:
        return applyCalibrationAreas(areas, calibration), fits
    else:
        return applyCalibrationAreas(areas, calibration)

def applyFromFile(filename: str, **kwargs) -> list[float]:
    # load the file
    bins, vals = Spectrum.read(filename, **kwargs)

    # apply the calibration
    return apply(bins, vals, **kwargs)