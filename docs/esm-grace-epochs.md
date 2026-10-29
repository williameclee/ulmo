# ESAESM2 fields for GRACE intervals

`esm2plmt` averages separate ESAESM2 components over specified GRACE intervals
and returns the array layout, units and midpoint dates used by `grace2plmt_new`.
It supports land-only, ocean-only and combined synthetic recovery experiments.

The [ESAESM2 archive](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/)
contains six-hour Stokes fields through degree/order 180 for 1995–2006. Component
identifiers are A (IB-corrected atmosphere), H (terrestrial hydrology), I (land
ice), O (ocean), S (solid Earth) and AOHIS (their sum). See the
[archive README](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/README.txt)
and [Dobslaw et al. (2015)](https://doi.org/10.1007/s00190-014-0787-8).

## Download one interval

Set `GRACE_FILE` to a local GRACE source file and export `IFILES` with the same
data root used by MATLAB. The fixed ESM location is
`$IFILES/ESM/ESAESM2/raw/COMPONENT/`; no input-folder argument is needed in
either the downloader or `esm2plmt`. Run from the repository root with Python
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

Configure ULMO and its dependencies on the MATLAB path, including the usual
`IFILES` Earth constants and Love-number data for SD conversion. Read the
interval from the same GRACE file used by the downloader:

```matlab
graceFile = string(getenv("GRACE_FILE"));
[~, ~, ~, intervals] = parsegracesourcefile(graceFile);
[esmPlmt, esmStdPlmt, dates] = esm2plmt(intervals, ...
    "Component", "O", "Unit", "SD", "Loutput", 60);
```

The processor searches only the requested component's fixed raw folder,
including subdirectories. It reports an error if `IFILES` is unset or that
folder is missing; it does not search the current working directory.

Intervals are an N-by-2 unzoned UTC datetime array, with exclusive ends. Every
six-hour sample in each interval receives equal weight. Missing, duplicate or
off-grid timestamps fail. The interval matches GRACE's reporting window; this
does not reproduce orbit-dependent observation weights or missing observation
days within that window.

| Output | Default representation |
| --- | --- |
| `esmPlmt` | N-by-M-by-4, columns degree, order, C and S |
| `esmStdPlmt` | Same labels; coefficient uncertainties are NaN |
| `dates` | N-by-1 MATLAB datenums at interval midpoints |

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
The source radius is read for unit conversion. Header metadata are not checked
or returned; the input is assumed to follow ESAESM2's fully normalized
coefficient convention. Missing or duplicate harmonic rows still fail.

## Save and plot components

The demo reads inputs from the same fixed location. Its output directory is
specified separately:

```matlab
outputFolder = fullfile(getenv("IFILES"), "ESM", "ESAESM2", "processed");
products = esm2plmt_demo(graceFile, outputFolder);
```

This saves native GRAV fields for H/I/O, native and degree-60 SD fields for
H/I/O/HI/HIO, and a two-column PNG comparing native with truncated H/I/O fields.
All MAT files contain the three outputs above, without a metadata variable.
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
- ESM supplies no formal uncertainties. NaN means unavailable, not zero error.
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

Run the MATLAB regression tests with `assertSuccess(runtests("unittests/esm2plmtTest.m"))`.
They exercise interval boundaries, averaging, coefficient sorting, unit
conversion, low degrees, layouts, padding, and malformed/incomplete inputs.
