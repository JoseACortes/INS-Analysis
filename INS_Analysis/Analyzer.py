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
                    'vals': None,
                },
                'Si2C1': {
                    'bins': None,
                    'vals': None,
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
        }

    def addSpectrums(self, specs, labels, **kwargs):
        for spec, label in zip(specs, labels):
            self.addSpectrum(spec, label, **kwargs)

    def calcPeakAreas(self, labels, **kwargs):
        res = {}
        for label in labels:
            bins = self.spectrums[label]['bins']
            vals = self.spectrums[label]['vals']
            res[label] = calcPeakAreas(bins, vals, **kwargs)
        return res
    
    def calibrate(
        self, 
        labels: list,
        Si1_method: str = 'original', 
        Si1_p0=None,
        Si2C1_method: str = 'original',
        Si2C1_p0=None,
        **kwargs):
        true_concentrations = {}
        for label in labels:
            true_concentrations[label] = self.spectrums[label]['true_comp']
        true_concentrations = list(true_concentrations.values())
        areas = self.calcPeakAreas(labels, **kwargs)
        return self.calibrate(
            areas=areas, 
            true_concentrations=true_concentrations,
            Si1_method=Si1_method,
            Si1_p0=Si1_p0,
            Si2C1_method=Si2C1_method,
            Si2C1_p0=Si2C1_p0,
            **kwargs
            )
    
    def applyCalibrationAreas(self, areas, **kwargs):
        return applyCalibrationAreas(areas=areas, calib_params = self.calibration, **kwargs)
    
    def apply(self, **kwargs):
        return apply(calib_params = self.calibration, **kwargs)
    
    def applyFromFile(self, **kwargs):
        return applyFromFile(calib_params = self.calibration, **kwargs)