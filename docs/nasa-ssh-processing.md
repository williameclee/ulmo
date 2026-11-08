# NASA-SSH reference-mission sea level anomalies

ULMO reads weekly NASA-SSH maps as `NASASSH` and calendar-month means as
`NASASSHMonthly`. The source variable `ssha` gives sea level anomalies
in metres relative to the DTU21 mean sea surface. Each 0.5° map combines
ten days of observations from the TOPEX/Poseidon, Jason and Sentinel-6
reference missions. Consecutive maps are seven days apart and share
observations. [NASA-SSH reference-mission simple grid, version 1](https://doi.org/10.5067/NSREF-SG0V1)

## Prepare the weekly maps

Set `NASA_SSH_INPUT_FOLDER` to the directory of source NetCDF files and
`IFILES` to the ULMO data directory. Run the following in MATLAB.

```matlab
weeklyPath = processSSHDataNasa(getenv('NASA_SSH_INPUT_FOLDER'));
[ssh, sigma, dates, lon, lat] = ssh2lonlatt('NASASSH');
```

The converter reads files named `NASA-SSH_alt_ref_simple_grid_v1_YYYYMMDD.nc`.
It checks the product, reference surface, units and grid, and matches the
central NetCDF epoch to the filename date. Fill values, nonfinite values and
values outside the source validity range become NaN. An entirely missing map
retains its epoch. The converter keeps the native grid and applies no
additional GIA, ocean-bottom-deformation or drift corrections.

The default output is `IFILES/SSH/NASASSH.mat`. To choose another destination,
pass a filename or path as the second argument. Both converters finish a
temporary v7.3 MAT file in the destination directory, then publish it with a
hard link. Publication uses MATLAB's JVM and requires a filesystem with
hard-link support. It preserves any existing destination, including one
created by a concurrent writer.
Pass `BeQuiet=true` to suppress progress and completion messages.

The weekly file contains single-precision `sshs` in longitude × latitude ×
epoch order, increasing `lon` and `lat` axes, a `dates` column of central
epochs and an empty `sshErrors`. The `metadata` structure records the source
product, DOI, reference surface, units, grid order, gridding method, map
coverage intervals and processing description. `creationDate` records when
the file was written, in UTC. By default, `ssh2lonlatt` returns latitude ×
longitude × epoch arrays in metres. Use `OutputFormat='ndgrid'` for the
on-disk order or `Unit='mm'` for millimetres.

## Prepare calendar-month means

```matlab
monthlyPath = processSSHMonthlyNasa(weeklyPath, ...
    fullfile(getenv('IFILES'), 'SSH', 'NASASSHMonthly.mat'), BeQuiet=true);
[ssh, sigma, dates, lon, lat] = ssh2lonlatt('NASASSHMonthly');
coverage = load(monthlyPath, 'mapCounts', 'expectedMapCounts', 'partialMonths');
```

Each finite weekly value has equal weight in the calendar month containing
its central epoch. The converter accumulates in double precision and writes
single-precision means on the native grid. A cell with no finite values has
a NaN mean. All calendar months between the first and last source epochs are
present, including empty months and partial edge months. Each output date
is halfway between consecutive month starts.

| Monthly variable | Meaning |
| --- | --- |
| `mapCounts` | Available source maps per month, including entirely missing maps. |
| `expectedMapCounts` | Scheduled weekly epochs per full calendar month, using the phase of the source dates. |
| `validMapCounts` | Finite weekly values per cell and month. |
| `incompleteCoverage` | True where a cell has fewer finite samples than the full month’s expected count. |
| `partialMonths` | True where scheduled epochs fall before the first or after the last input date. A false value can still accompany missing maps or cells. |
| `monthStarts`, `monthEnds` | Inclusive start and exclusive end of each calendar month. |
| `metadata` | Source dates, weekly indices per month, native metadata and averaging policy. |

The coverage flags leave the means intact. Use `partialMonths` to identify
truncated edge months, and `validMapCounts` or `incompleteCoverage` for
cell-level screening. The overlapping ten-day windows limit the interpretation
of these means as calendar-month observations. Neither NASA-SSH product
supplies an uncertainty estimate.

A request such as `ssh2lonlatt('NASASSH', 'midmonth')` also averages weekly
samples by calendar month. The reader returns only midpoints inside the
native date span, so it can omit edge months that the monthly converter
retains. If requested, the reader interpolates the spatial grid before
averaging. Use the prepared monthly file for the native-grid means and
sampling diagnostics.

## Comparisons with other altimetry products

| Product used for comparison | Grid spacing | Time sampling | Mapping and observations |
| --- | --- | --- | --- |
| NASA-SSH simple grid v1 | 0.5° | Seven-day interval, ten-day windows | Gaussian weighting of reference-mission observations |
| MEaSUREs v2205 | 1/6° | Five-day interval | Kriging of reference and complementary mission observations |
| DUACS DT2024 monthly SLA | 0.125° | Calendar-month means | All-satellite mapped anomalies |

NASA-SSH uses DTU21. The DT2024 SLA baseline is the 1993–2012 mean.
Align the anomaly baseline, time span, spatial coverage and averaging
procedure before comparing the products. For NASA-SSH, consult the saved
`metadata.griddingMethod` for the parameters recorded in the source file.
Shared satellite observations limit the independence of the products.
[MEaSUREs v2205](https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205) · [DUACS DT2024](https://duacs.cls.fr/duacs-system-description/operational-news/duacs-my-system-impacting-version-changes/nov-2024-duacs-dt-2024/)
