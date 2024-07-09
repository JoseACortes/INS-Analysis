# INS-Analysis
Python package for analyzing MINS results

# Use

```python
import INS-Analysis as insa
```

## Data Preperation
```
bins, vals = [], []
bins, vals = isna.read(filepath)
```


```
MINS = insa.Analyzer()
```

```
# the following three give the same result
# areas is a list of numbers (int or float)

areas = insa.calcPeakAreas(bins, vals)

areas, fits = insa.calcPeakAreas(bins, vals, return_fits=True)
```

This can also be called with the MINS object

```
areas = MINS.calcPeakAreas(bins, vals)
```


```
MINS.calibrate(areas, true_concentrations)
```

```
calib_params = MINS.calibrate(areas, true_concentrations)
calib_params = insa.calibrate(areas, true_concentrations)
```

```
pred_concentrations = MINS.applyCalibrationAreas(areas)
pred_concentrations = insa.applyCalibrationAreas(areas, calib_params)

pred_concentrations = MINS.apply(bins, vals)
pred_concentrations = insa.apply(bins, vals, calib_params)

pred_concentrations = MINS.applyFromFile(filenames)
pred_concentrations = insa.applyFromFile(filenames, calib_params)
```


# package structure

```
INS-Analysis/
├── README.md
├── setup.py
├── requirements.txt
├── LICENSE
├── .gitignore
├── Documenation/
├── Examples/
└── INS-Analysis/
    ├── __init__.py
    ├── apply.py
    ├── Analyzer.py
    ├── utils.py
    ├── areas/
    ├── Spectrum/
    └── calibration/
```