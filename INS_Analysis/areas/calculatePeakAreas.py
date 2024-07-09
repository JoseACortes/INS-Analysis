import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad

from ..tools import fitting_functions as ff
from ..tools import window_maker as wm

defaultPeakWindows = {
    'Si1': [1.6, 2.2],
    'C1Si2': [4.2, 4.8],
}

baselineFunctions = {
    'point_slope': {
        'fn': ff.point_slope,
        'bounds': {
            'lower': [None, None],
            'upper': [None, None],
            'starting_weights': [None, None],
            }
        },
    'point_slope_super': {
        'fn': ff.point_slope_super,
        'bounds': {
            'lower': [-1e100, 0, 'lower bound x'],
            'upper': [0, 'upper bound y', 'upper bound x'],
            'starting_weights': [0, 'lower bound y', 'center x']
            },
        }
}

peakFunctions = {
    'gaus': {
        'fn': ff.gaus,
        'bounds': {
            'lower': ['lower bound x', 0, 0],
            'upper': ['upper bound x', 'window height', 'upper bound sigma'],
            'starting_weights': [None, None, 'sigma']
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
    'C1Si2': {
        'targets': [4.44, 4.5],
        'window': defaultPeakWindows['C1Si2'],
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
        return None
    return geb_fn(target, geb['a'], geb['b'], geb['c'])

def ubsgm(target, geb, **kwargs):
    if geb is None:
        return None
    return geb_fn(target, geb['a'], geb['b'], geb['c']) * 2 
atbds = {
    'lower bound x': lbx,
    'upper bound x': ubx,
    'lower bound y': lby,
    'upper bound y': uby,
    'center x': cx,
    'window height': wh,
    'sigma': sgm,
}


def autobound(bounds, bins, vals, target, geb):
    ret = []
    for ab in bounds:
        if isinstance(ab, str):
            ret.append(atbds[ab](bins=bins, vals=vals, target=target, geb=geb))
        else:
            ret.append(ab)
    return ret


def theActualPeakAreaCalculation(bins, vals, config, returnFits=False, **kwargs):
    wave = np.array([bins, vals])
    geb = config['geb']
    # si1
    si1_config = config['Si1']
    si1_wave_bins, si1_wave_vals = wm.make_window(wave, si1_config['window'][0], si1_config['window'][1])
    si1_peak_config = si1_config['peaks']
    si1_baseline_config = si1_config['baseline']
    target = si1_config['targets'][0]

    si1_baseline_lower_bounds = autobound(bounds=si1_baseline_config['bounds']['lower'],bins=si1_wave_bins,vals=si1_wave_vals,target=target,geb=geb)
    si1_baseline_upper_bounds = autobound(bounds=si1_baseline_config['bounds']['upper'],bins=si1_wave_bins,vals=si1_wave_vals,target=target,geb=geb)
    si1_baseline_starting_weights = autobound(bounds=si1_baseline_config['starting_weights'],bins=si1_wave_bins,vals=si1_wave_vals,target=target,geb=geb)
    si1_baseline_function = si1_baseline_config['fn']
    si1_baseline_n_weights = len(si1_baseline_upper_bounds)

    si1_peak_lower_bounds = autobound(bounds=si1_peak_config['bounds']['lower'],bins=si1_wave_bins,vals=si1_wave_vals,target=target,geb=geb)
    si1_peak_upper_bounds = autobound(bounds=si1_peak_config['bounds']['upper'],bins=si1_wave_bins,vals=si1_wave_vals,target=target,geb=geb)
    si1_peak_starting_weights = autobound(bounds=si1_peak_config['starting_weights'],bins=si1_wave_bins,vals=si1_wave_vals,target=target,geb=geb)
    si1_peak_function = si1_peak_config['fn']
    si1_peak_n_weights = len(si1_peak_upper_bounds)

    si1_sum_function = ff.generate_compound_sum([si1_baseline_function, si1_peak_function], [si1_baseline_n_weights, si1_peak_n_weights])
    si1_sum_lower_bounds = si1_baseline_lower_bounds + si1_peak_lower_bounds
    si1_sum_upper_bounds = si1_baseline_upper_bounds + si1_peak_upper_bounds
    si_starting_weights = si1_baseline_starting_weights + si1_peak_starting_weights

    si1_final_weights, si1_covs = curve_fit(
        si1_sum_function,
        si1_wave_bins,
        si1_wave_vals,
        p0=si_starting_weights,
        maxfev=500000,
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
        args=si1_peak_final_weights
    )[0]

    # c1si2
    c1si2_config = config['C1Si2']
    c1si2_wave_bins, c1si2_wave_vals = wm.make_window(wave, c1si2_config['window'][0], c1si2_config['window'][1])
    c1si2_peak_config = c1si2_config['peaks']
    c1si2_baseline_config = c1si2_config['baseline']
    targets = c1si2_config['targets']

    c1si2_baseline_lower_bounds = autobound(bounds=c1si2_baseline_config['bounds']['lower'],bins=c1si2_wave_bins,vals=c1si2_wave_vals,target=targets,geb=geb)
    c1si2_baseline_upper_bounds = autobound(bounds=c1si2_baseline_config['bounds']['upper'],bins=c1si2_wave_bins,vals=c1si2_wave_vals,target=targets,geb=geb)
    c1si2_baseline_starting_weights = autobound(bounds=c1si2_baseline_config['starting_weights'],bins=c1si2_wave_bins,vals=c1si2_wave_vals,target=targets,geb=geb)
    c1si2_baseline_function = c1si2_baseline_config['fn']
    c1si2_baseline_n_weights = len(c1si2_baseline_upper_bounds)

    c1si2_peak_lower_bounds = autobound(bounds=c1si2_peak_config['bounds']['lower'],bins=c1si2_wave_bins,vals=c1si2_wave_vals,target=targets,geb=geb)
    c1si2_peak_upper_bounds = autobound(bounds=c1si2_peak_config['bounds']['upper'],bins=c1si2_wave_bins,vals=c1si2_wave_vals,target=targets,geb=geb)
    c1si2_peak_starting_weights = autobound(bounds=c1si2_peak_config['starting_weights'],bins=c1si2_wave_bins,vals=c1si2_wave_vals,target=targets,geb=geb)
    c1si2_peak_function = c1si2_peak_config['fn']
    c1si2_peak_n_weights = len(c1si2_peak_upper_bounds)

    c1si2_sum_function = ff.generate_compound_sum([c1si2_baseline_function, c1si2_peak_function], [c1si2_baseline_n_weights, c1si2_peak_n_weights])
    c1si2_sum_lower_bounds = c1si2_baseline_lower_bounds + c1si2_peak_lower_bounds
    c1si2_sum_upper_bounds = c1si2_baseline_upper_bounds + c1si2_peak_upper_bounds
    c1si2_starting_weights = c1si2_baseline_starting_weights + c1si2_peak_starting_weights

    c1si2_final_weights, c1si2_covs = curve_fit(
        c1si2_sum_function,
        c1si2_wave_bins,
        c1si2_wave_vals,
        p0=c1si2_starting_weights,
        maxfev=500000,
        bounds=(
            c1si2_sum_lower_bounds,
            c1si2_sum_upper_bounds
            )
    )

    c1si2_baseline_final_weights = c1si2_final_weights[:c1si2_baseline_n_weights]
    c1si2_peak_final_weights = c1si2_final_weights[c1si2_baseline_n_weights:c1si2_baseline_n_weights+c1si2_peak_n_weights]

    c1si2_peak_areas = []
    for i in range(len(targets)):
        c1si2_peak_areas.append(quad(
            c1si2_peak_function,
            c1si2_wave_bins[0],
            c1si2_wave_bins[-1],
            args=c1si2_peak_final_weights
        )[0])
    
    c1si2_peak_area = sum(c1si2_peak_areas)


    weights = {
        'Si1': {
            'baseline': si1_baseline_final_weights,
            'peak': si1_peak_final_weights,
        },
        'C1Si2': {
            'baseline': c1si2_baseline_final_weights,
            'peak': c1si2_peak_final_weights,
        }
    }


    fits = {
        'Si1': {
            'bins': bins,
            'baseline': si1_baseline_function(bins, *si1_baseline_final_weights),
            'peak': si1_sum_function(bins, *si1_final_weights),
        }, 
        'C1Si2': {
            'bins': bins,
            'baseline': c1si2_baseline_function(bins, *c1si2_baseline_final_weights),
            'peak': c1si2_sum_function(bins, *c1si2_final_weights),
        }
        }
    
    areas = {
        'Si1': si1_peak_area,
        'C1Si2': c1si2_peak_area,
    }

    if returnFits:
        return areas, fits
    return areas

windowlabellist = ['Si1', 'C1Si2']

def calcPeakAreas(
        bins, 
        vals, 
        peakWindows=defaultPeakWindows, # list or dict
        peakFunctions='gauss', # str, list or dict
        peakStartingWeights=None, # list of lists or dict
        peakUpperBounds=None, # list of lists or dict
        peakLowerBounds=None, # list of lists or dict
        baselineFunction='point_slope', # str, list or dict
        baselineStartingWeights=None, # list of lists or dict
        baselineUpperBounds=None, # list of lists or dict
        baselineLowerBounds=None, # list of lists or dict
        geb = None, # list
        returnFits:bool=False,
        **kwargs
        ): #a wrapper for theActualPeakAreaCalculation
    
    if isinstance(peakWindows, list):
        peakWindows = dict(zip(windowlabellist, peakWindows))
    if isinstance(peakFunctions, str):
        peakFunctions = dict(zip(windowlabellist, peakFunctions))
    if isinstance(peakStartingWeights, list):
        peakStartingWeights = dict(zip(windowlabellist, peakStartingWeights))
    if isinstance(peakUpperBounds, list):
        peakUpperBounds = dict(zip(windowlabellist, peakUpperBounds))
    if isinstance(peakLowerBounds, list):
        peakLowerBounds = dict(zip(windowlabellist, peakLowerBounds))
    if isinstance(baselineFunction, str):
        baselineFunction = dict(zip(windowlabellist, [baselineFunction]*len(windowlabellist)))
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
                config[key]['peaks'][i]['fn'] = peakFunctions[key]
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
            for i in range(len(config[key]['peaks'])):
                config[key]['baseline']['fn'] = baselineFunction[key]
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
        for key in geb:
            config[key]['geb'] = geb[key]

    return theActualPeakAreaCalculation(bins, vals, config, returnFits=returnFits, **kwargs)