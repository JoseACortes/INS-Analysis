from .Spectrum import read
from .areas import calcPeakAreas
from .calibration import calibrate
from .calibration import applyCalibrationAreas
from .process import apply
from .process import applyFromFile


class Analyzer():
    def __init__(self):
        self.calibration = None
        self.spectrums = {}

    def addSpectrum(self, spec, label, **kwargs):
        
        if isinstance(spec, str):
            spec = read(spec, **kwargs)

        self.spectrums[label] = {
            'bins': spec[0], 
            'vals': spec[1],
            'areas': {
                'Si1': None,
                'Si2C1': None,
            },
            'fits': {
                'Si1': {
                    'bins': None,
                    'baseline': None,
                    'peak': None,
                },
                'Si2C1': {
                    'bins': None,
                    'baseline': None,
                    'peak': None,
                },
            },
            'true_comp': {
                'Si1': None,
                'Si2C1': None,
            },
            'pred_comp': {
                'Si1': None,
                'Si2C1': None,
            },
            'for_calib': False,
            'area_calc_failed': False,
        }

    def addSpectrums(self, specs, labels, **kwargs):
        for spec, label in zip(specs, labels):
            self.addSpectrum(spec, label, **kwargs)

    def toDataFrame(self):
        import pandas as pd
        df = pd.DataFrame(self.spectrums).T
        return df.reset_index().rename(columns={'index':'label'})
    
    def calcPeakArea(self, label, returnFits=False, **kwargs):
        self.spectrums[label]['area_calc_failed'] = False
        bins = self.spectrums[label]['bins']
        vals = self.spectrums[label]['vals']
        try:
            res = calcPeakAreas(bins, vals, returnFits=returnFits, **kwargs)
        except:
            self.spectrums[label]['area_calc_failed'] = True
            return None
        areas = res
        fits = None
        if returnFits:
            areas, fits = areas
        self.spectrums[label]['areas']['Si1'] = areas['Si1']
        self.spectrums[label]['areas']['Si2C1'] = areas['Si2C1']
        if returnFits:
            self.spectrums[label]['fits']['Si1']['bins'] = fits['Si1']['bins']
            self.spectrums[label]['fits']['Si1']['baseline'] = fits['Si1']['baseline']
            self.spectrums[label]['fits']['Si1']['peak'] = fits['Si1']['peak']
            self.spectrums[label]['fits']['Si2C1']['bins'] = fits['Si2C1']['bins']
            self.spectrums[label]['fits']['Si2C1']['baseline'] = fits['Si2C1']['baseline']
            self.spectrums[label]['fits']['Si2C1']['peak'] = fits['Si2C1']['peak']
        return res

    def calcPeakAreas(self, labels, returnFits=False, **kwargs):
        for label in labels:
            self.calcPeakArea(label, returnFits=returnFits, **kwargs)
        ret = {}
        for label in labels:
            ret[label] = self.spectrums[label]['areas']
        if returnFits:
            ret_fits = {}
            for label in labels:
                ret_fits[label] = self.spectrums[label]['fits']
            return ret, ret_fits
        return ret
    
    def calibrate(
            self, 
            labels: list,
            concentrations=None,
            Si1_method: str = 'original', 
            Si1_p0=None,
            Si2C1_method: str = 'original',
            Si2C1_p0=None,
            **kwargs):
        
        for spectrum in self.spectrums.values():
            spectrum['for_calib'] = False
        
        if concentrations is None:
            concentrations = []
            for label in labels:
                concentrations.append(self.spectrums[label]['true_comp'])
        if isinstance(concentrations[0], dict):
            concentrations = [[c['Si1'], c['Si2C1']] for c in concentrations]
        
        areas = [self.spectrums[label]['areas'] for label in labels]
        areas = [[a['Si1'], a['Si2C1']] for a in areas]

        res = calibrate(
            areas=areas, 
            true_concentrations=concentrations,
            Si1_method=Si1_method,
            Si1_p0=Si1_p0,
            Si2C1_method=Si2C1_method,
            Si2C1_p0=Si2C1_p0,
            **kwargs
            )
        self.calibration = res

        for label, concentration in zip(labels, concentrations):
            self.spectrums[label]['true_comp']['Si1'] = concentration[0]
            self.spectrums[label]['true_comp']['Si2C1'] = concentration[1]
        for label in labels:
            self.spectrums[label]['for_calib'] = True
        return res
    
    def applyCalibrationAreas(self, labels=None, **kwargs):
        if labels is None:
            labels = list(self.spectrums.keys())
            labels = [label for label in labels if not self.spectrums[label]['area_calc_failed']] # removes failed area calculations
        areas = []
        for label in labels:
            a = self.spectrums[label]['areas']
            areas.append([a['Si1'], a['Si2C1']])
        res = applyCalibrationAreas(areas=areas, calibration = self.calibration, **kwargs)
        for _, label in enumerate(labels):
            self.spectrums[label]['pred_comp']['Si1'] = res[0][_]
            self.spectrums[label]['pred_comp']['Si2C1'] = res[1][_]
        ret = {}
        for label in labels:
            ret[label] = self.spectrums[label]['pred_comp']
        return ret
    
    def apply(self, **kwargs):
        return apply(calib_params = self.calibration, **kwargs)
    
    def applyFromFile(self, **kwargs):
        return applyFromFile(calib_params = self.calibration, **kwargs)