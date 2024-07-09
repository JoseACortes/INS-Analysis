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
    
    def calibrate(self, areas, **kwargs):
        self.calibration = calibrate(areas, **kwargs)
        return self.calibration
    
    def applyCalibrationAreas(self, areas, **kwargs):
        return applyCalibrationAreas(areas=areas, calib_params = self.calibration, **kwargs)
    
    def apply(self, **kwargs):
        return apply(calib_params = self.calibration, **kwargs)
    
    def applyFromFile(self, **kwargs):
        return applyFromFile(calib_params = self.calibration, **kwargs)