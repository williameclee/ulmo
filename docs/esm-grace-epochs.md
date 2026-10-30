# Cached ESAESM2 fields and monthly means

`esm2plmt` returns six-hour ESAESM2 fields or calendar-month means, with time
selection and the coefficient layout and units used by `grace2plmt_new`. Both
the native inputs and monthly products use NetCDF files on disk. It supports land-only, ocean-only
and combined synthetic recovery experiments.

The [ESAESM2 archive](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/)
contains six-hour Stokes fields through degree/order 180 for 1995–2006. Component
identifiers are A (IB-corrected atmosphere), H (terrestrial hydrology), I (land
ice), O (ocean), S (solid Earth) and AOHIS (their sum). See the
[archive README](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/README.txt)
and [Dobslaw et al. (2015)](https://doi.org/10.1007/s00190-014-0787-8).

## Download the full record

Export `IFILES` with the staging data root, which may be on an external drive.
The downloader stores GFC files under `$IFILES/ESM/ESAESM2/raw/COMPONENT/`.
This root can differ from the MATLAB `IFILES` used for prepared NetCDF files.
Run from the repository root with Python 3.9 or newer; no third-party Python
packages are needed:

```sh
# All six components, every six-hour sample from 1995 through 2006.
python3 data/download_esm.py
# Full time record for the components used in ocean/land recovery experiments.
python3 data/download_esm.py --components H I O
# Optional whole-year subset.
python3 data/download_esm.py --years 2002 2003 --components O
```

The default covers all 17,532 samples per component in the six-hour,
degree-180 `mtmshc` product. The separate three-hour, degree-360 product and
other ESAESM2 variants are not included. No GRACE file or epoch is required.
The programmatic entry point is
`download_dataset(components, data_root, years=YEARS)` with explicit inputs;
CLI/environment parsing is separate.

GFZ packages the samples in yearly `.tar.gz` archives. The downloader streams
one component-year at a time, retaining every expected sample without saving
the archive itself. Members need not be chronologically ordered. Each saved
file has a local SHA-256 hash and source URL/member record. These are local
integrity records, not publisher-supplied digests; the full gzip archive is
not retained or verified.

`$IFILES/ESM/ESAESM2/manifests/ESAESM2_YEAR.json` is saved after each component
completes. Restarting verifies all recorded sizes and hashes before skipping
completed component-years. Missing samples trigger another archive read;
modified files raise an error and are preserved. A failed component-year is
retried from its archive start, reusing matching files already saved. Completed
years/components survive later failures. Adding components preserves existing
manifest entries. Proof-of-concept samples share the same raw paths and can
be reused when downloading the complete record.

### Optional single-epoch proof-of-concept demo

The January 2003 proof-of-concept example demonstrates ingestion, conversion
and spatial recovery for one GRACE interval. To reproduce it, set `GRACE_FILE`
to a local GRACE source file:

```sh
python3 data/download_esm.py --grace-file "$GRACE_FILE" --components H I O
```

The downloader reads the exact start and exclusive end from GRCOF2 records,
rather than the inclusive day-of-year range in the filename. Boundaries must
align with the six-hour grid. `--grace-file` and `--years` are mutually
exclusive. The demo defaults to H/I/O when no components are supplied.
Its programmatic entry point remains
`download_epoch(grace_file, components, data_root)`.

The stream stops after every requested timestamp has been collected; later
intervals may require reading most of an annual archive.
`manifests/GRACE_FILENAME.json` records the GRACE input checksum, interval
and saved samples. Demo runs reread the relevant archives, require existing
files to match, and preserve completed components. Overlapping intervals
share raw files but retain separate source manifests. A changed GRACE input
or interval cannot reuse its previous manifest.

CSR RL06 January 2003 uses `[2003-01-01, 2003-02-01)`, with 124 samples per
component and midpoint `2003-01-16 12:00`. H, I and O together occupy
approximately 339 MiB after extraction. These figures describe the demo,
not the full production download.

## MATLAB processing

Configure ULMO and its dependencies on the MATLAB path, including `data/`
and the usual Earth constants and Love-number data. Raw GFC ingestion is a
separate operation:

```matlab
% ESM_SOURCE identifies the extracted files, possibly on an external drive.
sourceFolder = getenv("ESM_SOURCE");
for component = ["H", "I", "O"]
    prepareesm(sourceFolder, component);
end
% Return every available monthly epoch from the prepared ocean record.
[monthly, monthlyStd, monthlyDates] = esm2plmt("O", [], "midmonth", "Loutput", 60);
```

`prepareesm` searches the supplied folder recursively for the requested
component across all years, sorts timestamps, validates sampling and harmonic rows, and writes
`$IFILES/ESM/ESAESM2/ESAESM2_O.nc`. It reads at most 32 files per batch and
writes a temporary NetCDF-4 file, replacing the destination only on success.
An existing destination requires `ForceNew=true`. Reimport the entire desired
collection when adding data; imports replace rather than append. The original
GFC files are never modified and can be kept offline after preparation. Neither
the importer nor the production loader restricts data to one GRACE epoch.

### Prepared input and reader interface

The file has three dimensions: `lm` (the triangular list of degree/order
pairs), `cs` (cosine then sine), and `time`. MATLAB arrays use this order:

| Variable | Type and representation |
| --- | --- |
| `degree`, `order` | `uint16` vectors on `lm`, stored once |
| `cs` | `uint8` coordinate: 0 cosine, 1 sine |
| `coefficients` | M-by-2-by-N double Stokes values |
| `time` | N double seconds since 2000-01-01 00:00:00 UTC; proleptic Gregorian calendar |
| `uncertainty` (optional) | M-by-2-by-N double formal Stokes uncertainties |

`equatorRadius` is a scalar global attribute in metres. The importer verifies
that every source file uses the same reference radius. `coefficientUnit` is
`GRAV` for native inputs and identifies the converted units in derived files.
ESAESM2 has no formal errors, so prepared files omit `uncertainty`; the reader
returns NaNs for the selected samples. Coefficients retain double precision.

NetCDF-4 uses lossless deflate compression (level 4) and shuffle, with MATLAB
chunk sizes `[min(M,1024), 2, min(N,8)]`. Only coordinate labels use uint16;
labels are converted to doubles when reconstructing public lmcosi outputs.
The `lm` axis lists `(0,0), (1,0), (1,1), ...` through the native degree.
Explicit dimension declarations preserve a three-dimensional variable even
for one epoch. MATLAB's NetCDF interface presents dimensions in MATLAB array
order; C/Python tools may display the reversed dimension order.

The shared helpers expose lmcosi arrays by default while retaining compact
storage. `readshnc` returns M-by-4-by-N arrays with degree/order in the first
two columns. `OutputFormat="compact"` returns only the two coefficient columns;
the fourth output supplies the shared labels. Both layouts read time first
and load only the selected time/degree slab:

```matlab
[plmt, sigma, dates] = readshnc(inputFile, timelim, Loutput);
[coefficients, errors, dates, labels] = readshnc(inputFile, timelim, Loutput, ...
    "OutputFormat", "compact");
writeshnc(outputFile, plmt, sigma, dates, "GRAV");
% Equivalent compact input; degree/order are supplied separately.
writeshnc(outputFile, coefficients, errors, dates, "GRAV", "Labels", labels);
```

`writeshnc` infers labels from lmcosi input and verifies that they agree across
epochs and uncertainty arrays. Pass `[]` for uncertainty to omit that variable.
Both `createshnc` and `writeshnc` accept additional global attributes through
an alternating name/value cell, for example
`Attributes={"equatorRadius",6378137,'source','ESAESM2'}`. Attribute names must
be unique; `coefficientUnit` is reserved for the unit argument. `esm2plmt`
explicitly requests compact reads for conversion and averaging.

`esm2plmt` expects the fixed prepared NetCDF file and never searches for or
reads GFC files, even with `ForceNew=true`.

The interface is `esm2plmt(component, timelim, timestep, ...)`. Component
defaults to `'O'`, limits to `[]`, and time step to `''`. Text inputs accept
both char and string scalars.

```matlab
timelim = [datetime(2003, 1, 1), datetime(2003, 1, 31, 23, 59, 59)];
[raw, rawStd, rawDates] = esm2plmt("O", timelim);
[esmPlmt, esmStdPlmt, dates] = esm2plmt("O", timelim, "midmonth", ...
    "Unit", "SD", "Loutput", 60);
```

`timelim` is an increasing two-element datetime or datenum vector, selecting
output timestamps inclusively; empty returns all dates. Raw mode allows gaps
and reads the time coordinate first, identifies the inclusive selection,
then reads only the corresponding coefficient/uncertainty slab and requested
harmonic rows. Empty selections do not read either coefficient variable.
`'midmonth'` averages all
samples in `[month start, next month start)` and returns the exact midpoint.
A missing sample invalidates that month with NaN coefficients/uncertainties
and an `IncompleteEpoch` warning, including partially downloaded boundary
months and missing interior months. Download complete calendar months for
this option. Calendar months need not match GRACE reporting windows.

### Averaged caches

Like `ssh2lonlatt`, the reader returns an existing averaged cache immediately.
Cache hits also read time first and load only the selected dates/degrees.
On a miss, it reads the prepared NetCDF file one epoch at a time, converts units,
averages every epoch at native degree, saves, and then selects/formats output.
**Averaged caches always cover all epochs at native degree**, regardless of
`timelim` or `Loutput`. Unit conversions of raw requests are in memory; the
reader does not save additional full raw copies for each unit.

Files live in `$IFILES/ESM/ESAESM2/`. Examples:

| Product | Filename |
| --- | --- |
| Prepared native Stokes input | `ESAESM2_O.nc` |
| Monthly Stokes | `ESAESM2_O_GRAV-Tmidmonth.nc` |
| Monthly surface density | `ESAESM2_O_SD_ISSM-Tmidmonth.nc` |
| Monthly potential | `ESAESM2_O_POT-Tmidmonth.nc` |

Averaged caches use the same compact `degree`, `order`, `cs`, `time`,
`coefficients` and `uncertainty` schema, at native degree. They do not repeat
labels for each time or for the uncertainty array. Names use `-TSTEP` without
a method suffix;
SD includes the Love-number source. A cache hit works even when the prepared
input is temporarily unavailable.

`prepareesm(..., "ForceNew", true)` replaces the prepared input and invalidates
all its component's `-T*.nc` caches, across units and Love-number sources.
`esm2plmt(..., "ForceNew", true)` only recomputes averages from the prepared
NetCDF file. If replacing a NetCDF file manually, force a rebuild of affected
averaged caches. `SaveData=false` prevents cache writes, while allowing reads;
combine with `ForceNew=true` to recompute without changing disk products.
Neither option changes the prepared input. `BeQuiet=true` suppresses messages.

Existing MAT files are not read or deleted. Run `prepareesm` to generate
NetCDF inputs; averaged `.nc` caches are built on demand. The importer only
invalidates derived NetCDF caches for the component it replaces.

### Epoch generation and future GRACE support

`calendarmonthepochs([firstDate, lastDate])` supplies N-by-2 `[start, end)`
windows spanning the inclusive date range. Equal endpoints select one month.
`averageepochs(nativeDates, data, epochs, expectedStep, ...)` accepts those windows explicitly,
returns their exact midpoints, and has no calendar-month logic. Windows can
overlap and their boundaries need not align with the native sampling grid.
The helper averages generic arrays with time in their third dimension; the
caller handles harmonic labels and layout conversion. The first two outputs
are means and midpoint dates, retaining the interface on the GRACE epoch
averaging branch. Its third output is the uncertainty estimate of each mean.

`expectedStep` accepts a duration or a numeric number of days.
`CheckCompleteness=false` is the default and averages available nonmissing
values. Set it to true to reject incomplete timestamp grids or missing values.
`SamplingOrigin` specifies the grid phase (default: first native timestamp).
ESM callers supply `hours(6)`, an explicit UTC grid origin and
`CheckCompleteness=true`. The calendar helper takes only the first/last dates and makes no
sampling-rate or product assumptions. The demo uses these same averaging and
completeness rules for its explicit GRACE source window.

This uses the epoch-array convention on `feature/grace-epoch-averaging`.
A future `'GRACE'` provider can supply actual GRACE windows where available
and fallback windows elsewhere without changing the averaging or formatting
code. The public time-step options are `''` (or `""`) and `'midmonth'`; `'GRACE'` is
not enabled in this PR. Any future reference-product choices will also need
to distinguish their caches.

### Output representation

| Output | Default representation |
| --- | --- |
| `esmPlmt` | N-by-M-by-4, columns degree, order, C and S |
| `esmStdPlmt` | Same labels; supplied formal errors or empirical standard errors of means |
| `dates` | N-by-1 MATLAB datenums at sample times or month midpoints |

M is `(Loutput+1)*(Loutput+2)/2`: 16,471 rows at degree 180 or 1,891 at degree
60. Omitting `Loutput` retains the native degree. A larger requested degree pads
coefficients with NaN, without changing the valid degree/order labels. Options
also support `OutputFormat="traditional"` (M-by-4-by-N),
`TimeFormat="datetime"`, and `Unit="GRAV"` or `"POT"`.

`GRAV` returns dimensionless Stokes coefficients without requiring gravity
conversion data. `SD` uses the existing `convertgravity` implementation, source
radius, project density of 5517 kg/m³, and the selected load Love numbers
(default ISSM, with CF degree 1). Its units are kg/m², numerically equal to mm
of fresh water. `POT` follows the existing radius-times-Stokes convention.
The importer reads a positive finite source radius for unit conversion.
Descriptive header metadata are not checked or returned; inputs are assumed
to follow ESAESM2's fully normalized coefficient convention. The importer
rejects missing or duplicate harmonic rows.

### Uncertainty definition

Raw `esmStdPlmt` contains formal coefficient errors like `graceStdPlmt`, or
NaN when absent. ESAESM2 supplies no formal errors. Averaged errors use one
of two estimates, separately for each coefficient:

- If every included sample has a supplied formal error,
  `sigmaMean = sqrt(sum(sigma.^2))/n`.
- Otherwise, the empirical standard error is `std(samples,0)/sqrt(n)`, with
  sample variance normalized by n−1. At least two valid samples are required;
  fewer yield NaN. Known and unknown formal errors are not mixed or set to zero.

Mean and uncertainty use the same valid-data mask. Rejected/empty means have
NaN uncertainty. Both formulas assume independent samples. Six-hour signals
can be temporally correlated, so the empirical estimate is descriptive and
may underestimate uncertainty; it does not measure model bias or provide an
autocorrelation-adjusted trend uncertainty.

HI/HIO demo means and empirical errors are calculated from summed raw samples,
retaining their observed cross-component covariance. When all formal errors
are supplied, component errors are combined in quadrature under the separate
assumption of independent component errors.

## Save and plot components

The demo reads inputs from the same fixed location. Its output directory is
specified separately. It selects the exact GRACE reporting window from the
prepared NetCDF files and averages those samples in the demo, so it also supports windows
that differ from calendar months:

```matlab
outputFolder = fullfile(getenv("IFILES"), "ESM", "ESAESM2", "processed");
graceFile = string(getenv("GRACE_FILE"));
products = esm2plmt_demo(graceFile, outputFolder);
```

This saves native GRAV fields for H/I/O, native and degree-60 SD fields for
H/I/O/HI/HIO, and a two-column PNG comparing native with truncated H/I/O fields.
Demo NetCDF files use the same compact schema as averaged caches; the figure
reader reconstructs the required lmcosi field from coefficients and labels.
For CSR January 2003, the GRACE interval is the complete calendar month, so
the demo and `esm2plmt(component, [], "midmonth", ...)` produce the same coefficients.
`HI` and `HIO` add coefficient columns while preserving one set of harmonic
labels. Use `Loutput=45` to choose
a different truncated bandwidth or `MakeFigure=false` to skip plotting.

The figure synthesizes the complete fields on a 0.5-degree grid; it does not
mask out spectral coastal tails. Each row shares a symmetric color scale.
Display saturation is stated in the caption and can be changed without
reprocessing the coefficients:

```matlab
plot_esm_epoch(products, fullfile(outputFolder, "comparison.png"), ...
    "ColorLimits", [500, 500, 100]);
```

The figure requires `plm2xyz` dependencies and MATLAB's `coastlines` example
data. Plotting does not change the coefficients or the source data.

## Scientific conventions

- ESM fields are anomalies. No WGS84 subtraction, temporal demeaning, filtering,
  real-data low-degree replacement, de-aliasing or GIA removal is applied.
- C00 and degree 1 are retained for every component. Individual components
  need not have zero total mass.
- ESM supplies no formal uncertainties. Raw errors are NaN; means can use
  empirical standard errors from sample spread as described above.
- ISSM conversion matches the ULMO processing convention; it does not prove
  exact inversion of the Love-number spectrum used to generate ESM.
- The supplied ocean field retains its uniform mass-compensation layer. No
  SLE-based redistribution is introduced.
- HIO omits atmosphere and solid Earth and is not a closed whole-Earth mass
  budget. It is also not a GSM-like residual: no DEAL has been subtracted.
- SD conversion of S describes equivalent mass, not a physical surface load.
- The single-epoch proof of concept demonstrates ingestion and spatial recovery. Trend experiments need
  a multi-epoch record and a fixed truth across the candidate buffers.

## Validation

```sh
python3 -B -m pytest unittests/test_esm_download.py unittests/test_esm_verify.py
python3 data/verify_esm_epoch.py \
  "$IFILES/ESM/ESAESM2/manifests/${GRACE_FILE##*/}.json" \
  --output "$IFILES/ESM/ESAESM2/verification"
```

The downloader and verifier tests require pytest and run offline using
temporary directories and mocked archive responses. They check full-record
scheduling, leap years, resumable annual downloads, local corruption detection,
and optional GRACE-interval selection. The demo verifier checks file hashes,
complete six-hour sampling and all degree-180
harmonics, and independently averages each component using Python's standard
library. It writes `reference_COMPONENT.bin` as little-endian float64 values,
with C/S interleaved in ascending degree/order, and `raw_validation.json`.
Compare a saved native GRAV product in MATLAB as follows:

```matlab
verificationFolder = fullfile(getenv("IFILES"), "ESM", "ESAESM2", "verification");
fid = fopen(fullfile(verificationFolder, "reference_O.bin"), "r", "ieee-le");
reference = fread(fid, [2, 16471], "double")';
fclose(fid);
coefficients = readshnc(products.O.nativeGRAV);
max(abs(coefficients(:, 3:4, 1) - reference), [], "all")
```

Run the MATLAB regression tests with:

```matlab
assertSuccess(runtests(["unittests/esm2plmtTest.m", "unittests/epochStatisticsTest.m", ...
    "unittests/shncTest.m"]));
```

They exercise inclusive time selection, calendar boundaries and leap years,
prepared NetCDF and monthly cache reuse, full cache coverage before truncation, forced
refresh and invalidation, disabled saving, missing months, unit conversion,
low degrees, layouts, padding, malformed inputs, explicit irregular/overlapping
epochs, arbitrary cadences and formal/empirical uncertainty estimates through saved
monthly and demo products. NetCDF tests check stored integer labels, singleton
dimensions, exact round trips for lmcosi and compact layouts, inferred labels,
additional attributes, preservation of existing files on rejected writes, and
the actual time/degree slab arguments used by native and cached reads.
