import numpy as np

def gaus(x,x0,a,sigma):
    return a*np.exp(-(x-x0)**2/(sigma**2))

def double_gaus(x,x0,a,sigma,x1,a1,sigma1):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+a1*np.exp((-((x-(x1))))**2/(2*sigma1**2))

def lorentz(x,x0,a,gamma):
    return a/(1+((x-x0)/gamma)**2)

def pvoigt(x,x0,a,sigma,gamma):
    return (1-a)*gaus(x,x0,a,sigma)+a*lorentz(x,x0,a,gamma)

def point_slope(x,a,b):
    return a*x+b

def point_slope_super(x,a,b,x0):
    return a*(x-x0)+b

def exp_falloff(x,x0,a,p,b):
    return (a*np.exp(-p*(x-x0)))+b

def fat_tail(x,x0,a,p,b):
    return (a/(x-x0)**p)+b

from scipy.stats import chi2

def chi_2(x,x0,a,df):
    return chi2.pdf(x-x0, df)*a

def char(x,x0,x1):
    return np.ma.masked_inside(x,x0,x1).astype(int).astype(float)
    
def x(x, x0):
    return x0*x

def const(x, x0):
    return x0

def original_calibration(x, k1, k2):
    c = x[0]
    si = x[1]
    return (c-k1*si)/k2

def generate_compound_sum(fitting_functions, weight_lens):
    """
    fitting_functions: list of fitting functions
    weight_lens: list of lengths of weights
    """
    
    # generate slice masks for each fitting function from weight_lens
    indexs = []
    for _ in range(len(weight_lens)):
        indexs.append((sum(weight_lens[:_]),sum(weight_lens[:_+1])))

    def compound_sum(x, *args):
        """
        x: array of x values
        args: list of weights
        """
        compound_sum = 0
        for _ in range(len(fitting_functions)):
            weights = args[indexs[_][0]:indexs[_][1]]
            compound_sum += fitting_functions[_](x, *weights)
        return compound_sum
    return compound_sum

def fixed_weight_function(fit_function, weights, fixed_weights_mask):
    new_weights_index = np.argwhere(fixed_weights_mask == False)
    new_weights_index = new_weights_index.flatten()
    new_weights = weights[new_weights_index]
    
    def fixed_weight_fn(x, *newargs):

        weights[new_weights_index] = newargs    
        return fit_function(x, *weights)
    
    return fixed_weight_fn, new_weights