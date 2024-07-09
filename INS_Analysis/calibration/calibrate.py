from scipy.optimize import curve_fit
from ..tools import fitting_functions as ff

calibrations_indexing = {
    'original': {
        'function': ff.original_calibration,
        'params': ['k1', 'k2'],
        'p0': [1, 1]
    }
}

def calibrate(areas: list, true_concentrations: list, method: str = 'original', p0=None ,**kwargs) -> float:

    if p0 is None:
        p0 = calibrations_indexing[method]['p0']

    curve = calibrations_indexing[method]['function']

    r_params, c_covs = curve_fit(curve, areas, true_concentrations, p0=p0, maxfev=10000)
    
    weights = {}
    for i, param in enumerate(calibrations_indexing[method]['params']):
        weights[param] = r_params[i]
        
    calibration = {
        'method': method,
        'weights': weights
    }
    return calibration

def applyCalibrationAreas(areas: list, calibration: dict) -> list[float]:
    method = calibration['method']
    curve = calibrations_indexing[method]['function']
    weights = calibration['weights']
    return curve(areas, **weights)
