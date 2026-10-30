# Cached ESAESM2 fields and monthly means

`esm2plmt` returns six-hour ESAESM2 fields or calendar-month means, with time
selection and the coefficient layout and units used by `grace2plmt_new`. Both
the native inputs and monthly products use MAT files on disk. It supports land-only, ocean-only
and combined synthetic recovery experiments.

The [ESAESM2 archive](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/)
contains six-hour Stokes fields through degree/order 180 for 1995–2006. Component
identifiers are A (IB-corrected atmosphere), H (terrestrial hydrology), I (land
ice), O (ocean), S (solid Earth) and AOHIS (their sum). See the
[archive README](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/README.txt)
and [Dobslaw et al. (2015)](https://doi.org/10.1007/s00190-014-0787-8).

## Download one interval

Set `GRACE_FILE` to a local GRACE source file and export `IFILES` with the same
staging data root for the download, which may be on an external drive.
The downloader stores GFC files under `$IFILES/ESM/ESAESM2/raw/COMPONENT/`.
This download root can differ from the MATLAB `IFILES` used for prepared MAT files. Run from the repository root with Python
3.9 or newer:

```sh
python3 data/download_esm_epoch.py \
  --grace-file "$GRACE_FILE" --components H I O
```

No third-party Python dependencies are needed. The downloader reads the exact
start and exclusive end from the GRCOF2 records, rather than interpreting the
inclusive day-of-year range in the filename. Boundaries must align with the
six-hour grid. The six-hour series covers 1995–2006; the separate three-hour,
degree-360 product is not supported.

GFZ packages the samples in yearly `.tar.gz` archives. The downloader streams
each required archive and saves only the requested samples under
`$IFILES/ESM/ESAESM2/raw/H`, `raw/I` and `raw/O`. Members need not be
chronologically ordered. The stream
stops only after every expected timestamp has been collected; later intervals
may require reading most of the annual archive. It does not retain or verify
the full yearly gzip archive.

`$IFILES/ESM/ESAESM2/manifests/GRACE_FILENAME.json` records the GRACE input checksum, interval, component file
names, timestamps, source URLs, archive member names and local SHA-256 hashes.
Completed components are recorded individually, and subsequent component
downloads preserve earlier entries. Rerunning an existing component reads its
archive again; existing files must match the downloaded content. Overlapping
GRACE intervals share the same raw sample files and retain separate manifests
for their source filenames. A manifest cannot be reused after changing its
GRACE source contents or interval. The local checksums
are integrity records, not publisher-supplied digests.

For example, CSR RL06 January 2003 uses `[2003-01-01, 2003-02-01)`, with 124
samples per component and midpoint `2003-01-16 12:00`. H, I and O together occupy
approximately 339 MiB after extraction.

## MATLAB processing

Configure ULMO and its dependencies on the MATLAB path, including `data/`
and the usual Earth constants and Love-number data. Raw GFC ingestion is a
separate operation:

```matlab
% ESM_SOURCE identifies the extracted files, possibly on an external drive.
sourceFolder = getenv("ESM_SOURCE");
prepareesm(sourceFolder, "O");
```

`prepareesm` searches the supplied folder recursively for the requested
component, sorts timestamps, validates sampling and harmonic rows, and writes
`$IFILES/ESM/ESAESM2/ESAESM2_O.mat`. It reads at most 32 files per batch and
writes a temporary v7.3 MAT file, replacing the destination only on success.
An existing destination requires `ForceNew=true`. Reimport the entire desired
collection when adding data; imports replace rather than append. The original
GFC files are never modified and can be kept offline after preparation.

### Prepared input and reader interface

The prepared file contains:

| Variable | Representation |
| --- | --- |
| `esmPlmt` | N-by-M-by-4 native-degree Stokes lmcosi, columns degree/order/C/S |
| `dates` | N-by-1 sorted, unique, unzoned six-hour datetime timestamps |
| `equatorRadius` | N-by-1 source radii in metres, used for unit conversion |

No all-NaN uncertainty array is stored for ESAESM2. The reader also accepts an
optional `esmStdPlmt` with formal Stokes uncertainties in the same layout.
Use MAT v7.3 for partial array reads. `esm2plmt` expects this fixed prepared
file and never searches for or reads GFC files, even with `ForceNew=true`.

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
and reads only the selected dates and degrees. `'midmonth'` averages all
samples in `[month start, next month start)` and returns the exact midpoint.
A missing sample invalidates that month with NaN coefficients/uncertainties
and an `IncompleteEpoch` warning, including partially downloaded boundary
months and missing interior months. Download complete calendar months for
this option. Calendar months need not match GRACE reporting windows.

### Averaged caches

Like `ssh2lonlatt`, the reader returns an existing averaged cache immediately.
On a miss, it reads the prepared MAT file one epoch at a time, converts units,
averages every epoch at native degree, saves, and then selects/formats output.
**Averaged caches always cover all epochs at native degree**, regardless of
`timelim` or `Loutput`. Unit conversions of raw requests are in memory; the
reader does not save additional full raw copies for each unit.

Files live in `$IFILES/ESM/ESAESM2/`. Examples:

| Product | Filename |
| --- | --- |
| Prepared native Stokes input | `ESAESM2_O.mat` |
| Monthly Stokes | `ESAESM2_O_GRAV-Tmidmonth.mat` |
| Monthly surface density | `ESAESM2_O_SD_ISSM-Tmidmonth.mat` |
| Monthly potential | `ESAESM2_O_POT-Tmidmonth.mat` |

Averaged caches contain `esmPlmt`, `esmStdPlmt` and datetime `dates`, in
native-degree, time-first form. Names use `-TSTEP` without a method suffix;
SD includes the Love-number source. A cache hit works even when the prepared
input is temporarily unavailable.

`prepareesm(..., "ForceNew", true)` replaces the prepared input and invalidates
all its component's `-T*.mat` caches, across units and Love-number sources.
`esm2plmt(..., "ForceNew", true)` only recomputes averages from the prepared
MAT file. If replacing a MAT file manually, force a rebuild of affected
averaged caches. `SaveData=false` prevents cache writes, while allowing reads;
combine with `ForceNew=true` to recompute without changing disk products.
Neither option changes the prepared input. `BeQuiet=true` suppresses messages.

Earlier draft versions wrote unit-specific raw caches such as
`ESAESM2_O_GRAV.mat`. These are no longer used: run `prepareesm` once to create
the canonical input and invalidate the old averaged products. The import
also ensures old NaN-only uncertainty caches are rebuilt with the new
empirical estimate; obsolete raw copies may be removed separately.

### Epoch generation and future GRACE support

`calendarmonthepochs(nativeDates)` supplies N-by-2 `[start, end)` windows.
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
`CheckCompleteness=true`. The calendar helper takes only dates and makes no
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
prepared MAT files and averages those samples in the demo, so it also supports windows
that differ from calendar months:

```matlab
outputFolder = fullfile(getenv("IFILES"), "ESM", "ESAESM2", "processed");
graceFile = string(getenv("GRACE_FILE"));
products = esm2plmt_demo(graceFile, outputFolder);
```

This saves native GRAV fields for H/I/O, native and degree-60 SD fields for
H/I/O/HI/HIO, and a two-column PNG comparing native with truncated H/I/O fields.
All demo MAT files contain the three outputs above, without a metadata variable.
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
- A single epoch tests ingestion and spatial recovery. Trend experiments need
  a multi-epoch record and a fixed truth across the candidate buffers.

## Validation

```sh
python3 -B -m pytest unittests/test_esm_download.py
python3 data/verify_esm_epoch.py \
  "$IFILES/ESM/ESAESM2/manifests/${GRACE_FILE##*/}.json" \
  --output "$IFILES/ESM/ESAESM2/verification"
```

The downloader tests require pytest and use its temporary-directory and
monkeypatch fixtures. The verifier checks file hashes, complete six-hour sampling and all degree-180
harmonics, and independently averages each component using Python's standard
library. It writes `reference_COMPONENT.bin` as little-endian float64 values,
with C/S interleaved in ascending degree/order, and `raw_validation.json`.
Compare a saved native GRAV product in MATLAB as follows:

```matlab
verificationFolder = fullfile(getenv("IFILES"), "ESM", "ESAESM2", "verification");
fid = fopen(fullfile(verificationFolder, "reference_O.bin"), "r", "ieee-le");
reference = fread(fid, [2, 16471], "double")';
fclose(fid);
saved = load(products.O.nativeGRAV);
max(abs(squeeze(saved.esmPlmt(1,:,3:4)) - reference), [], "all")
```

Run the MATLAB regression tests with:

```matlab
assertSuccess(runtests(["unittests/esm2plmtTest.m", "unittests/epochStatisticsTest.m"]));
```

They exercise inclusive time selection, calendar boundaries and leap years,
prepared MAT and monthly cache reuse, full cache coverage before truncation, forced
refresh and invalidation, disabled saving, missing months, unit conversion,
low degrees, layouts, padding, malformed inputs, explicit irregular/overlapping
epochs, arbitrary cadences and formal/empirical uncertainty estimates through saved
monthly and demo products.
