# GRACE input formats

`grace2plmt_new` reads monthly GRCOF2 spherical-harmonic files. Set
`ORIGINALGRACEDATA` to the raw-data root, or use `$GRACEDATA/raw`
(falling back to `$IFILES/GRACE/raw`). Place files under `<release>/<center>`.
Processed files and logs go under `GRACEDATA` or `$IFILES/GRACE`.

| Center | Release | Filename pattern | Input bandwidth |
| --- | --- | --- | --- |
| CSR | RL05 | `GSM-2_2*-2*_*_UTCSR_0060_0005*` | 60 |
| GFZ | RL05 | `GSM-2_2*-2*_*_EIGEN_G---_005a` | 90 |
| JPL | RL05 | `GSM-2_2*-2*_*_JPLEM_0000_0005*` | 60 or 90, stored at 90 |
| CSR, GFZ, JPL | RL06 | `GSM-2_*_BA01_06*` / `GSM-2_*_BB01_06*` | 60 / 96 |
| CSR mascon | RL06 | `GSU-2_*_B---_06*` | 720 |

RL05 and mascon bandwidths override `Ldata` with a warning when needed.
Use `Loutput` to truncate or pad the returned arrays. Missing coefficients
and uncertainties remain **NaN**, with valid degree/order labels in both
arrays. In particular, JPL RL05 months of degree 60 do not supply measurements
at degrees 61–90. Use `Loutput=60` when downstream work needs complete months
at a common bandwidth. RL04 loading is not implemented.

```matlab
% Default RL05 behavior: no degree-1, C20 or C30 replacements.
[p, sigma, dates] = grace2plmt_new('JPL', 'RL05', 90, 'SD', 'Loutput', 60);

% Raw mascon updates relative to GGM05C, with replacements disabled.
[p, sigma, dates] = grace2plmt_new('CSR mascon', 'RL06', 720, 'GRAV', ...
    'Deg1Correction', false, 'C20Correction', false, 'C30Correction', false);
```

## Reference fields and corrections

Legacy headers supply GM and radius on the `EARTH` record; RL06 supplies them
as YAML attributes. Missing or invalid reference metadata raises an error.
This preserves GFZ RL05's radius of 6378136.46 m instead of substituting
CSR's 6378136.30 m.

GSM processing retains the existing WGS84 zonal subtraction. CSR mascon GSU
files are **updates relative to GGM05C**, so no WGS84 zonals are subtracted.
With replacements enabled, TN-14 full-field C20/C30 values are first expressed
relative to GGM05C: subtract C20 = −4.841694573200e−4 and C30 =
9.571647583412e−7. These are the zero-tide model coefficients, also listed as
means in TN-14. Standard deviations are unchanged by subtracting this fixed
reference. Degree-1 replacement uses CSR TN-13 (GGM05C degree 1 is zero).
See the [CSR GSU signal definition](https://www2.csr.utexas.edu/grace/RL06_mascons.html),
[GGM05C coefficients](https://download.csr.utexas.edu/pub/grace/GGM05/GGM05C.ICGEM)
and [TN-14](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt).

Degree-1/C20/C30 replacement defaults to **true for RL06**, including mascons.
For RL05, omitted or empty flags default to **false**; explicitly requesting
`true` raises `ULMO:grace2plmt:UnsupportedCorrection`. Release-compatible
RL05 replacement series are not implemented. Setting flags to false avoids
loading correction datasets. Missing monthly replacements retain the input
coefficient and are recorded in the processing log.

GSU output is not the fully corrected CSR NetCDF mascon product: this reader
does not restore GAD, remove GIA, apply the gridded ellipsoidal correction,
or remove a temporal mean. Choose a common anomaly baseline before comparing
GSM and GSU time series.

## Caches and validation

RL05 and CSR mascon processed files use an `_inputV2.mat` suffix so caches
created by the initial reader cannot silently retain incorrect radii, NaN
indices or mascon reference offsets. Older caches are left intact; the first
call with the revised reader requires raw inputs even with `ForceNew=false`.
Standard RL06 GSM cache names are unchanged.

With ULMO and its Slepian dependencies on the MATLAB path, run:

```matlab
results = runtests('unittests/graceInputFormatsTest.m');
assertSuccess(results);
```

The tests create temporary GRACE, correction and Love-number fixtures, restore
the environment afterward, and do not require downloads or personal datasets.
