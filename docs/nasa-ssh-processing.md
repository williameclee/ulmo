# NASA-SSH reference-mission simple grid

Use product name `NASASSH` in ULMO. This identifies
`NASA_SSH_REF_SIMPLE_GRID_V1`, DOI https://doi.org/10.5067/NSREF-SG0V1.
It is sea-surface-height anomaly in metres relative to DTU21, using the
TOPEX/Poseidon, Jason and Sentinel-6 reference-mission series.

| Product | Grid | Time sampling | Mapping / observations |
| --- | --- | --- | --- |
| NASA-SSH simple grid v1 | 0.5° | Every 7 days, overlapping 10-day windows | Gaussian spatial weighting; reference missions only |
| MEaSUREs v2205 | 1/6° | Every 5 days | Kriging; reference and complementary missions |
| Local DT2024 source | 0.125° | Calendar-month means | DUACS all-satellite mapped SLA |

The NASA download inspected on 2026-09-29 contains 1,377 weekly epochs,
1999-12-27 through 2026-05-11. The collection itself begins in 1992;
this download is a subset. The local DT2024 NetCDF identifies its SLA
baseline as 1993–2012. NASA uses DTU21, so absolute anomaly offsets should
not be interpreted as product errors without a common reference period.

These products can support basin-scale variability and trend comparisons,
but share satellite observations and are not independent measurements.
Match time periods and temporal aggregation, document the spatial support,
and check sensitivity to coverage and smoothing. A finer interpolated grid
does not add resolution. Preserve the product's missing cells and epochs.
Standard altimetric corrections are already upstream; the converter adds
no GIA, ocean-bottom-deformation, drift, or other trend correction.

The NASA catalogue describes a 100 km Gaussian width; the actual downloaded
file metadata specifies `roi=600000.0, sigma=175000.0, neighbours=500`.
The converter records the source gridding metadata without regridding.
The catalogue calls this L4, while the examined NetCDF global attribute
says Level 3. Product identity is established from its short name and DOI.

## Conversion

```matlab
processSSHDataNasa(getenv('NASA_SSH_INPUT_FOLDER'));  % directory of source NetCDF files
[ssh, sigma, dates, lon, lat] = ssh2lonlatt('NASASSH');
```

The default destination is `fullfile(getenv('IFILES'),'SSH','NASASSH.mat')`.
An explicit output path can be passed as the second converter argument.
Existing files are not overwritten. Conversion writes a temporary v7.3
MAT file and publishes it only after all input maps pass validation.

The MAT file contains `sshs` (single precision, longitude × latitude × time,
metres), increasing `lon` and `lat`, column `datetime` values in `dates`,
empty `sshErrors`, and provenance `metadata`. Time comes from NetCDF
`time` (seconds since 1990-01-01), cross-checked against filenames.
Fill and invalid values become NaN. Entirely missing maps retain their
epochs. No uncertainty field is available; counts are not uncertainties.
Source filenames, sizes, modification timestamps, coverage intervals and
finite-cell counts are stored in metadata.

```matlab
% ULMO's existing temporal/spatial interpolation, output in millimetres:
[ssh, sigma, dates, lon, lat] = ssh2lonlatt( ...
    'NASASSH', 'midmonth', 1, [], 'Unit', 'mm');
```

For weekly NASA data, the current `interptemporal` linear/midmonth path
samples at month midpoints; it does not compute calendar-month means.
For comparison with DT2024 monthly means, explicitly aggregate over
calendar months and account for overlapping source windows. Native data
are retained so that a scientific aggregation choice can be made later.

## Sources

- NASA collection: https://podaac.jpl.nasa.gov/dataset/NASA_SSH_REF_SIMPLE_GRID_V1
- MEaSUREs: https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205
- DUACS DT2024: https://duacs.cls.fr/duacs-system-description/operational-news/duacs-my-system-impacting-version-changes/nov-2024-duacs-dt-2024/

## Installed download validation (2026-09-29)

All 1,377 output maps and NaN masks were verified against the NetCDF sources.
Maximum single-precision rounding error was 2.98014e-8 m. Nine all-missing
maps were preserved: 2001-12-10, 2003-11-24, 2005-09-26, 2006-11-06,
2006-11-13, 2013-04-01, 2013-09-09, 2020-02-03 and 2020-02-10.
Reader checks passed for native loading, coordinate order, units, dates,
spatial/temporal interpolation, cache reuse and the installed data file.

## Calendar-month product

`NASASSHMonthly` is a separate product made with `processSSHMonthlyNasa`:

```matlab
processSSHMonthlyNasa;  % reads IFILES/SSH/NASASSH.mat
[ssh, sigma, dates, lon, lat] = ssh2lonlatt('NASASSHMonthly');
% Already monthly: 'midmonth' preserves these values and dates.
[ssh, sigma, dates, lon, lat] = ssh2lonlatt('NASASSHMonthly','midmonth',1);
```

Each cell is the equal-weight mean of finite weekly values whose central
NetCDF timestamps lie in the calendar month. Accumulation uses double
precision; results are single-precision metres on the native 0.5° grid.
Zero valid values give NaN. Dates are exact halfway points between month
starts, including fractional days for odd-length months. The weekly
`NASASSH` product and its interpolation behavior are unchanged.

The monthly MAT file additionally contains:

- `validMapCounts`: number of finite maps per cell and month (uint8).
- `mapCounts`: total available source epochs per month, including empty maps.
- `expectedMapCounts`: number of weekly epochs expected in the full month,
  using the source's seven-day sampling phase.
- `incompleteCoverage`: true where valid counts are below expected counts;
  this also flags permanently missing cells. This flag does not mask means.
- `partialMonths`: true where expected weekly epochs extend beyond the source
  date span. December 1999 and May 2026 are retained and flagged.
- `monthStarts`, `monthEnds`: inclusive start and exclusive end boundaries.
- `metadata`: aggregation policy, source dates/indices and native provenance.

Coverage variables are accessible with `load` or `matfile`; `ssh2lonlatt`
returns its usual five outputs. Counts are sampling diagnostics, not an
uncertainty estimate or a count of independent observations. Means remain
approximate calendar-month means because the 10-day source windows overlap.
No minimum-count threshold is imposed: downstream analyses can choose one
using the saved counts and flags. For complete-month comparisons, omit the
two flagged edge months.
