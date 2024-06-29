# INS-Analysis
Python package for analyzing MINS results

# Use

```python
import INS-Analysis as insa
```

## Data Preperation
```
bins, vals = [], []
bins, vals = isna.Spectrum.readMCTAL(mctal_filepath)
bins, vals = isna.Spectrum.readMCA(mca_filepath)
```


```
MINS = insa.Analyzer()
```

```
# the following three give the same result
# areas is a list of numbers (int or float)

areas = insa.Calc_Peak_areas(bins, vals)

areas, fits = insa.Calc_Peak_areas(bins, vals, return_fits=True)
```

This can also be called with the MINS object

```
areas = MINS.Calc_Peak_areas(bins, vals)
```


```
MINS.calibrate(areas, true_concentrations)
```

```
calib_params = MINS.calibrate(areas, true_concentrations)
calib_params = insa.calibrate(areas, true_concentrations)
```

```
pred_concentrations = MINS.evaluate_areas(areas)
pred_concentrations = insa.evaluate_areas(areas, calib_params)

pred_concentrations = MINS.evaluate(bins, vals)
pred_concentrations = insa.evaluate(bins, vals, calib_params)

pred_concentrations = MINS.apply(filenames)
pred_concentrations = insa.apply(filenames, calib_params)
```