import Spectrum as Spectrum
import calibration as calibration
import areas
import process

calcPeakAreas = areas.calcPeakAreas
calibrate = calibration.calibrate
applyCalibrationAreas = calibration.applyCalibrationAreas
apply = process.apply
applyFromFile = process.applyFromFile

class Analyzer():
    def __init__(self):
        self.calibration = None

    def calcPeakAreas(self, bins, vals, **kwargs):
        return calcPeakAreas(bins, vals, **kwargs)
    
    def calibrate(self, areas, **kwargs):
        self.calibration = calibrate(areas, **kwargs)
        return self.calibration
    
    def applyCalibrationAreas(self, areas, **kwargs):
        return applyCalibrationAreas(areas=areas, calib_params = self.calibration, **kwargs)
    
    def apply(self, **kwargs):
        return apply(calib_params = self.calibration, **kwargs)
    
    def applyFromFile(self, **kwargs):
        return applyFromFile(calib_params = self.calibration, **kwargs)