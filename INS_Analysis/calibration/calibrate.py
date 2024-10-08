from scipy.optimize import curve_fit
from ..tools import fitting_functions as ff
import numpy as np

two_variable_calibrations_indexing = {
    'original': {
        'function': ff.original_calibration,
        'params': ['k1', 'k2'],
        'p0': [1, 1]
    }
}

single_variable_calibrations_indexing = {
    'original': {
        'function': ff.k1,
        'params': ['k1'],
        'p0': [1]
    }
}

def calibrate(
        areas: list, 
        true_concentrations: list[list[float]], 
        Si1_method: str = 'original', 
        Si1_p0=None,
        Si2C1_method: str = 'original',
        Si2C1_p0=None,
        **kwargs):
    
    true_concentrations = np.array(true_concentrations)
    areas = np.array(areas)

    if Si1_p0 is None:
        Si1_p0 = single_variable_calibrations_indexing[Si1_method]['p0']
    Si1_curve = single_variable_calibrations_indexing[Si1_method]['function']
    
    r_params, c_covs = curve_fit(Si1_curve, areas[:, 0], true_concentrations[:,  0], p0=Si1_p0, maxfev=10000)
    
    Si1_weights = {}
    for i, param in enumerate(single_variable_calibrations_indexing[Si1_method]['params']):
        Si1_weights[param] = r_params[i]

    if Si2C1_p0 is None:
        Si2C1_p0 = two_variable_calibrations_indexing[Si2C1_method]['p0']
    Si2C1_curve = two_variable_calibrations_indexing[Si2C1_method]['function']
    r_params, c_covs = curve_fit(Si2C1_curve, areas, true_concentrations[:, 1], p0=Si2C1_p0, maxfev=10000)

    Si2C1_weights = {}
    for i, param in enumerate(two_variable_calibrations_indexing[Si2C1_method]['params']):
        Si2C1_weights[param] = r_params[i]
        

    calibration = {
        'methods': {
            'Si1':Si1_method,
            'Si2C1':Si2C1_method
            },
        'weights': {
            'Si1': Si1_weights,
            'Si2C1': Si2C1_weights
            }
    }
    return calibration

def applyCalibrationAreas(areas: list, calibration: dict):
    Si1_curve = single_variable_calibrations_indexing[calibration['methods']['Si1']]['function']
    Si2C1_curve = two_variable_calibrations_indexing[calibration['methods']['Si2C1']]['function']
    Si1_weights = calibration['weights']['Si1']
    Si2C1_weights = calibration['weights']['Si2C1']
    Si1_areas = [area[0] for area in areas]
    Si1_predictions = Si1_curve(Si1_areas, **Si1_weights).tolist()
    Si2C1_predictions = Si2C1_curve(areas, **Si2C1_weights).tolist()
    return [Si1_predictions, Si2C1_predictions]
