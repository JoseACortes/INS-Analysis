import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad

from ..tools import fitting_functions as ff
from ..tools import window_maker as wm

defaultPeakWindows = {
    'Si1': [1.6, 2.2],
    'Si2C1': [4.2, 4.8],
}

baselineFunctions = {
    'point_slope': {
        'fn': ff.point_slope,
        'bounds': {
            'lower': [-np.inf, -np.inf],
            'upper': [np.inf, np.inf],
            'starting_weights': [1, 1],
            }
        },
    'point_slope_super': {
        'fn': ff.point_slope_super,
        'bounds': {
            'lower': [-1e100, 0, 'lower bound x'],
            'upper': [0, 'upper bound y', 'upper bound x'],
            'starting_weights': [0, 'lower bound y', 'center x']
            },
        },
    'exp_falloff': {
        'fn': ff.exp_falloff,
        'bounds': {
            'lower': ['lower bound x', 0, 0, 'lower bound y'],
            'upper': ['upper bound x', np.inf, np.inf, 'upper bound y'],
            'starting_weights': ['center x', 'window height', .1, 'lower bound y']
            }
        },
    'fat_tail': {
        'fn': ff.fat_tail,
        'bounds': {
            'lower': [0, 0, 0, 0],
            'upper': ['lower bound x', np.inf, 1, 'upper bound y'],
            'starting_weights': [0, 'lower bound y', .1, 'lower bound y']
            }
        },
}

peakFunctions = {
    'gaus': {
        'fn': ff.gaus,
        'bounds': {
            'lower': ['lower bound x', 0, 0],
            'upper': ['upper bound x', 'window height', 'upper bound sigma'],
            'starting_weights': ['center x', 'window height', 'sigma']
            }
        }
}

default_peak_area_config = {
    'Si1': {
        'targets': [1.78],
        'window': defaultPeakWindows['Si1'],
        'peaks': [peakFunctions['gaus']],
        'baseline': baselineFunctions['point_slope'],
    },
    'Si2C1': {
        'targets': [4.44, 4.5],
        'window': defaultPeakWindows['Si2C1'],
        'peaks': [peakFunctions['gaus'], peakFunctions['gaus']],
        'baseline': baselineFunctions['point_slope'],
    },
    'geb': None,
}

geb_fn = ff.geb

def lbx(bins, **kwargs):
    return min(bins)

def ubx(bins, **kwargs):
    return max(bins)

def lby(vals, **kwargs):
    return min(vals)

def uby(vals, **kwargs):
    return max(vals)

def cx(bins, **kwargs):
    return (min(bins) + max(bins)) / 2

def wh(vals, **kwargs):
    return (max(vals) - min(vals))*1.1

def sgm(target, geb, **kwargs):
    if geb is None:
        return 1
    return geb_fn(target, geb['a'], geb['b'], geb['c'])

def ubsgm(target, geb, **kwargs):
    if geb is None:
        return np.inf
    return geb_fn(target, geb['a'], geb['b'], geb['c']) * 2 
atbds = {
    'lower bound x': lbx,
    'upper bound x': ubx,
    'lower bound y': lby,
    'upper bound y': uby,
    'center x': cx,
    'window height': wh,
    'sigma': sgm,
    'upper bound sigma': ubsgm,
}


def autobound(bounds, bins, vals, target, geb):
    ret = []
    for ab in bounds:
        if isinstance(ab, str):
            ret.append(atbds[ab](bins=bins, vals=vals, target=target, geb=geb))
        else:
            ret.append(ab)
    return ret

###example config
# {
#   'element 1': {
#           targets:[], # Where you think the peak is centered
#           window:[], # minimum and maximum energy to regress calculate area
#           peaks:[{peak_config_1}, {peak_config_2}],
#           baseline:{baseline_config},
#                   },
#   'element 2': {},
#   'geb': {'a': -0.0073, 'b': 0.0078, 'c': 0},
#   'maxfev': 50000,
#   }

def theActualPeakAreaCalculation(bins, vals, config, returnFits=False, **kwargs):
    wave = np.array([bins, vals])
    geb = config['geb']
    # si1
    si1_config = config['Si1']
    si1_wave_bins, si1_wave_vals = wm.make_window(wave, si1_config['window'][0], si1_config['window'][1])
    si1_peak_config = si1_config['peaks']
    si1_baseline_config = si1_config['baseline']
    
    
    targets = si1_config['targets']

    si1_baseline_lower_bounds = autobound(bounds=si1_baseline_config['bounds']['lower'],bins=si1_wave_bins,vals=si1_wave_vals,target=targets[0],geb=geb)
    si1_baseline_upper_bounds = autobound(bounds=si1_baseline_config['bounds']['upper'],bins=si1_wave_bins,vals=si1_wave_vals,target=targets[0],geb=geb)
    si1_baseline_starting_weights = autobound(bounds=si1_baseline_config['bounds']['starting_weights'],bins=si1_wave_bins,vals=si1_wave_vals,target=targets[0],geb=geb)
    si1_baseline_function = si1_baseline_config['fn']
    si1_baseline_n_weights = len(si1_baseline_upper_bounds)


    si1_peak_upper_bounds = []
    si1_peak_lower_bounds = []
    si1_peak_starting_weights = []
    si1_peak_functions = []
    si1_peak_n_weights = []

    for _, peak in enumerate(si1_peak_config):
        peak_lower_bounds = autobound(bounds=peak['bounds']['lower'],bins=si1_wave_bins,vals=si1_wave_vals,target=targets[_],geb=geb)
        peak_upper_bounds = autobound(bounds=peak['bounds']['upper'],bins=si1_wave_bins,vals=si1_wave_vals,target=targets[_],geb=geb)
        peak_starting_weights = autobound(bounds=peak['bounds']['starting_weights'],bins=si1_wave_bins,vals=si1_wave_vals,target=targets[_],geb=geb)
        peak_function = peak['fn']
        peak_n_weights = len(peak_upper_bounds)

        si1_peak_upper_bounds += peak_upper_bounds
        si1_peak_lower_bounds += peak_lower_bounds
        si1_peak_starting_weights += peak_starting_weights
        si1_peak_functions.append(peak_function)
        si1_peak_n_weights.append(peak_n_weights)

    
    si1_peak_function = ff.generate_compound_sum(si1_peak_functions, si1_peak_n_weights)
    si1_peak_n_weights = sum(si1_peak_n_weights)

    si1_sum_function = ff.generate_compound_sum([si1_baseline_function, si1_peak_function], [si1_baseline_n_weights, si1_peak_n_weights])
    si1_sum_lower_bounds = si1_baseline_lower_bounds + si1_peak_lower_bounds
    si1_sum_upper_bounds = si1_baseline_upper_bounds + si1_peak_upper_bounds
    si_starting_weights = si1_baseline_starting_weights + si1_peak_starting_weights
    
    maxfev = 50000
    if 'maxfev' in config:
        maxfev = config['maxfev']
    
    si1_final_weights, si1_covs = curve_fit(
        si1_sum_function,
        si1_wave_bins,
        si1_wave_vals,
        p0=si_starting_weights,
        maxfev=maxfev,
        bounds=(
            si1_sum_lower_bounds,
            si1_sum_upper_bounds
            )
    )

    si1_baseline_final_weights = si1_final_weights[:si1_baseline_n_weights]
    si1_peak_final_weights = si1_final_weights[si1_baseline_n_weights:si1_baseline_n_weights+si1_peak_n_weights]

    si1_peak_area = quad(
        si1_peak_function,
        si1_wave_bins[0],
        si1_wave_bins[-1],
        args=tuple(si1_peak_final_weights.tolist())
    )[0]

    # si2C1

    si2C1_config = config['Si2C1']
    si2C1_wave_bins, si2C1_wave_vals = wm.make_window(wave, si2C1_config['window'][0], si2C1_config['window'][1])
    si2C1_peak_config = si2C1_config['peaks']
    si2C1_baseline_config = si2C1_config['baseline']
    targets = si2C1_config['targets']

    si2C1_baseline_lower_bounds = autobound(bounds=si2C1_baseline_config['bounds']['lower'],bins=si2C1_wave_bins,vals=si2C1_wave_vals,target=targets[0],geb=geb)
    si2C1_baseline_upper_bounds = autobound(bounds=si2C1_baseline_config['bounds']['upper'],bins=si2C1_wave_bins,vals=si2C1_wave_vals,target=targets[0],geb=geb)
    si2C1_baseline_starting_weights = autobound(bounds=si2C1_baseline_config['bounds']['starting_weights'],bins=si2C1_wave_bins,vals=si2C1_wave_vals,target=targets[0],geb=geb)
    si2C1_baseline_function = si2C1_baseline_config['fn']
    si2C1_baseline_n_weights = len(si2C1_baseline_upper_bounds)

    si2C1_peak_upper_bounds = []
    si2C1_peak_lower_bounds = []
    si2C1_peak_starting_weights = []
    si2C1_peak_functions = []
    si2C1_peak_n_weights = []

    for _, peak in enumerate(si2C1_peak_config):

        peak_lower_bounds = autobound(bounds=peak['bounds']['lower'],bins=si2C1_wave_bins,vals=si2C1_wave_vals,target=targets[_],geb=geb)
        peak_upper_bounds = autobound(bounds=peak['bounds']['upper'],bins=si2C1_wave_bins,vals=si2C1_wave_vals,target=targets[_],geb=geb)
        peak_starting_weights = autobound(bounds=peak['bounds']['starting_weights'],bins=si2C1_wave_bins,vals=si2C1_wave_vals,target=targets[_],geb=geb)
        peak_function = peak['fn']
        peak_n_weights = len(peak_upper_bounds)

        si2C1_peak_upper_bounds += peak_upper_bounds
        si2C1_peak_lower_bounds += peak_lower_bounds
        si2C1_peak_starting_weights += peak_starting_weights
        si2C1_peak_functions.append(peak_function)
        si2C1_peak_n_weights.append(peak_n_weights)

    si2C1_peak_function = ff.generate_compound_sum(si2C1_peak_functions, si2C1_peak_n_weights)
    si2C1_peak_n_weights = sum(si2C1_peak_n_weights)

    si2C1_sum_function = ff.generate_compound_sum([si2C1_baseline_function, si2C1_peak_function], [si2C1_baseline_n_weights, si2C1_peak_n_weights])
    si2C1_sum_lower_bounds = si2C1_baseline_lower_bounds + si2C1_peak_lower_bounds
    si2C1_sum_upper_bounds = si2C1_baseline_upper_bounds + si2C1_peak_upper_bounds
    si2C1_starting_weights = si2C1_baseline_starting_weights + si2C1_peak_starting_weights

    si2C1_final_weights, si2C1_covs = curve_fit(
        si2C1_sum_function,
        si2C1_wave_bins,
        si2C1_wave_vals,
        p0=si2C1_starting_weights,
        maxfev=maxfev,
        bounds=(
            si2C1_sum_lower_bounds,
            si2C1_sum_upper_bounds
            )
    )

    si2C1_baseline_final_weights = si2C1_final_weights[:si2C1_baseline_n_weights]
    si2C1_peak_final_weights = si2C1_final_weights[si2C1_baseline_n_weights:si2C1_baseline_n_weights+si2C1_peak_n_weights]

    si2C1_peak_areas = []
    for i in range(len(targets)):
        si2C1_peak_areas.append(quad(
            si2C1_peak_function,
            si2C1_wave_bins[0],
            si2C1_wave_bins[-1],
            args=tuple(si2C1_peak_final_weights.tolist())
        )[0])
    
    si2C1_peak_area = sum(si2C1_peak_areas)


    weights = {
        'Si1': {
            'baseline': si1_baseline_final_weights,
            'peak': si1_peak_final_weights,
        },
        'Si2C1': {
            'baseline': si2C1_baseline_final_weights,
            'peak': si2C1_peak_final_weights,
        }
    }

    fits = {
        'Si1': {
            'bins': bins,
            'baseline': si1_baseline_function(np.array(bins), *si1_baseline_final_weights).tolist(),
            'peak': si1_sum_function(np.array(bins), *si1_final_weights).tolist(),
        }, 
        'Si2C1': {
            'bins': bins,
            'baseline': si2C1_baseline_function(np.array(bins), *si2C1_baseline_final_weights).tolist(),
            'peak': si2C1_sum_function(np.array(bins), *si2C1_final_weights).tolist(),
        }
        }
    
    areas = {
        'Si1': si1_peak_area,
        'Si2C1': si2C1_peak_area,
    }

    if returnFits:
        return areas, fits
    return areas

windowlabellist = ['Si1', 'Si2C1']

common_fns = {
    'point_slope': ff.point_slope,
    'point_slope_super': ff.point_slope_super,
    'exp_falloff': ff.exp_falloff,
    'fat_tail': ff.fat_tail,
    'chi_2': ff.chi_2,
    'char': ff.char,
    'x': ff.x,
    'const': ff.const,
    'gaus': ff.gaus,
}


def calcPeakAreas(
        bins, 
        vals, 
        peakWindows=defaultPeakWindows, # list or dict
        peakFunctions='gaus', # str, list or dict
        peakStartingWeights=None, # list of lists or dict
        peakUpperBounds=None, # list of lists or dict
        peakLowerBounds=None, # list of lists or dict
        baselineFunction='point_slope', # str, list or dict
        baselineStartingWeights=None, # list of lists or dict
        baselineUpperBounds=None, # list of lists or dict
        baselineLowerBounds=None, # list of lists or dict
        geb = None, # list
        maxfev=None,
        returnFits:bool=False,
        **kwargs
        ): #a wrapper for theActualPeakAreaCalculation
    if isinstance(peakWindows, list):
        peakWindows = dict(zip(windowlabellist, peakWindows))

    if isinstance(peakFunctions, str):
        peakFunctions = [[peakFunctions], [peakFunctions, peakFunctions]]
    if isinstance(peakFunctions, list):
        peakFunctions = dict(zip(windowlabellist, peakFunctions))

    if isinstance(peakStartingWeights, list):
        peakStartingWeights = dict(zip(windowlabellist, peakStartingWeights))
    if isinstance(peakUpperBounds, list):
        peakUpperBounds = dict(zip(windowlabellist, peakUpperBounds))
    if isinstance(peakLowerBounds, list):
        peakLowerBounds = dict(zip(windowlabellist, peakLowerBounds))


    if isinstance(baselineFunction, str):
        baselineFunction = [baselineFunction]*len(windowlabellist)
    if isinstance(baselineFunction, list):
        baselineFunction = dict(zip(windowlabellist, baselineFunction))

    if isinstance(baselineStartingWeights, list):
        baselineStartingWeights = dict(zip(windowlabellist, baselineStartingWeights))
    if isinstance(baselineUpperBounds, list):
        baselineUpperBounds = dict(zip(windowlabellist, baselineUpperBounds))
    if isinstance(baselineLowerBounds, list):
        baselineLowerBounds = dict(zip(windowlabellist, baselineLowerBounds))
    if isinstance(geb, list):
        geb = dict(zip(windowlabellist, geb))

    config = default_peak_area_config

    if peakWindows is not None:
        for key in peakWindows:
            config[key]['window'] = peakWindows[key]

    if peakFunctions is not None:
        for key in peakFunctions:
            for i in range(len(config[key]['peaks'])):
                if isinstance(peakFunctions[key][i], str):
                    config[key]['peaks'][i]['fn'] = common_fns[peakFunctions[key][i]]
                else:
                    config[key]['peaks'][i]['fn'] = peakFunctions[key][i]

    if peakStartingWeights is not None:
        for key in peakStartingWeights:
            for i in range(len(config[key]['peaks'])):
                config[key]['peaks'][i]['bounds']['starting_weights'] = peakStartingWeights[key]
    if peakUpperBounds is not None:
        for key in peakUpperBounds:
            for i in range(len(config[key]['peaks'])):
                config[key]['peaks'][i]['bounds']['upper'] = peakUpperBounds[key]
    if peakLowerBounds is not None:
        for key in peakLowerBounds:
            for i in range(len(config[key]['peaks'])):
                config[key]['peaks'][i]['bounds']['lower'] = peakLowerBounds[key]

    if baselineFunction is not None:
        for key in baselineFunction:
            if isinstance(baselineFunction[key], str):
                config[key]['baseline'] = baselineFunctions[baselineFunction[key]]
            elif isinstance(baselineFunction[key], dict):
                config[key]['baseline'] =baselineFunction[key]

    if baselineStartingWeights is not None:
        for key in baselineStartingWeights:
            for i in range(len(config[key]['peaks'])):
                config[key]['baseline']['bounds']['starting_weights'] = baselineStartingWeights[key]
    if baselineUpperBounds is not None:
        for key in baselineUpperBounds:
            for i in range(len(config[key]['peaks'])):
                config[key]['baseline']['bounds']['upper'] = baselineUpperBounds[key]
    if baselineLowerBounds is not None:
        for key in baselineLowerBounds:
            for i in range(len(config[key]['peaks'])):
                config[key]['baseline']['bounds']['lower'] = baselineLowerBounds[key]
    if geb is not None:
        config['geb'] = geb
    
    if maxfev is not None:
        config['maxfev'] = maxfev

    return theActualPeakAreaCalculation(bins, vals, config, returnFits=returnFits, **kwargs)