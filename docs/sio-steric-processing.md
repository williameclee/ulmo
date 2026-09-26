# SIO Roemmich–Gilson processing

`processStericDataSio` uses the same density, climatology, component-density,
vertical-integration, and aggregation functions as the EN4/CMEMS pipelines.
It reads the 2019-release `RG_ArgoClim_Temperature_2019.nc` and
`RG_ArgoClim_Salinity_2019.nc`, followed by `RG_ArgoClim_YYYYMM_2019.nc`
extensions. Do not mix these with the older Atlas `Temp.nc` / `Psal.nc` files.

## Source interpretation

- Temperature: in-situ ITS-90 degrees Celsius. NetCDF says
  `degree celcius (ITS-90)` but does not explicitly distinguish in-situ from
  potential temperature. The interpretation is supported by NASA PO.DAAC's
  [HOMaGE SIO reader](https://github.com/podaac/HOMaGE/blob/main/mod_dl_cv_SIO.jl),
  which uses `gsw_ct_from_t` after adding the mean and anomaly.
  For the downloaded files, the 2019 temperature mean is numerically identical
  to the Atlas temperature mean used by that reader (maximum difference zero).
- Salinity: NetCDF explicitly specifies `Practical Salinity Scale 78`.
  This is dimensionless practical salinity, commonly called PSU, not g/kg.
- Pressure: dbar, with 58 sample levels from 2.5 to 1975 dbar. Integration
  extends to the 2000-dbar boundary, converted to positive depth by latitude.
- Time: fractional calendar months since 2004-01-01; `.5` is the actual
  calendar-month midpoint, including half-days and leap years.
- Mean plus monthly anomaly reconstructs each absolute T/S field.
  Both temperature and salinity bathymetry masks and mapping pressure limits
  are applied; missing columns remain NaN rather than zero sea level.

The [SIO product page](https://sio-argo.ucsd.edu/RG_Climatology.html) describes
2004–2018 base fields and later monthly extensions. Its potential-temperature
plot is a derived diagnostic and does not identify the raw T variable as PT.

The TEOS-10 sequence is `gsw_SA_from_SP` (Absolute Salinity, g/kg),
`gsw_CT_from_t` (Conservative Temperature, Celsius), then `gsw_rho`.
The default reference period is January 2004–December 2018. The shared
pipeline averages monthly density, SA and CT over that period and integrates
`(rhoClim/rho - 1) * layerThickness`, using positive-depth midpoint layers.
Thermosteric density holds SA at its reference value; halosteric density
holds CT at its reference value. These are the existing shared-pipeline
component definitions; their sum need not exactly equal total steric height.

## Running

```matlab
processStericDataSio(inputFolder, outputFolder, ...
    fullfile(outputFolder, 'SIO-StericSeaLevel.mat'), ...
    [datetime(2004,1,1), datetime(2018,12,31)], UseParallel=false);
```

Outputs are `SIO-MyyyyMM.mat`, the reference climatology, and the aggregate.
The aggregate contains coordinates, datetime timestamps, calendar bounds,
source metadata, and total/thermosteric/halosteric sea level in metres.
The full column is shallower than 2000 m, so shallow fields equal full fields
and deep fields are NaN.

After validation, install the aggregate at
`$IFILES/STERIC/SIO/SIO-StericSeaLevel.mat` and archive old interpolated
`SIO-StericSeaLevel-*` caches. `steric2lonlatt('SIO')` then reads the rebuild.
Monthly stages can resume; use `ForceNew=true` if raw files or scientific
settings change. New extensions are discovered on every run.
