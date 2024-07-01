import numpy as np

from ..tools import fitting_functions
from scipy.optimize import curve_fit

defaultPeakWindows = {
    'si0_window': [1.6, 2.2],
    'c0si1_window': [4.2, 4.8],
}

methods = {
    'original': {
        'baseline_function': None,
        'baseline_bounds': {
            'lower': [-1e100, 0],
            'upper': [0, np.inf]
        },
        'baseline_starting_weights': [0, 0],
        'peak_function': None,
        'peak_bounds': None,
        'peak_starting_weights': None,
    }
}

def calcPeakAreas(bins, vals, peakWindows=defaultPeakWindows, peakFunctions, peakStartingWeights, returnFits=False):
    pass