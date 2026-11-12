# Independent GRACE earthquake corrections

`grace2seismic` estimates a permanent step for each earthquake separately in
its own circular cap, including both land and ocean. Each fit uses the full
available GRACE record with a quadratic background and annual/semiannual
harmonics. The static event fields are converted to spherical harmonics
before addition; coefficients in different Slepian bases are never added.

```matlab
product = {'CSR', 'RL06', 60};
events = {'Sumatra', [142.373, 38.297]};
eventDates = [datetime(2004, 12, 26), datetime(2011, 3, 11)];
[correction, sigma, dates, eventFields] = grace2seismic( ...
    product, events, 20, eventDates, 60, TimeFormat="datetime");
```

`earthquakes` accepts a scalar name, a string or cell vector of names, an
`nEvents`-by-2 coordinate matrix, or a cell vector mixing scalar names with
numeric `[longitude, latitude]` rows. Every entry is one event. Dates and
vector radii pair in input order; a scalar radius is shared. The fourth
output preserves event order. Unknown names and malformed pairs fail.

Event times round to the nearest UTC midnight **before both fitting and
interval weighting**. Noon rounds forward. Zoned datetimes are converted
to UTC; unzoned times and datenums are interpreted as UTC. Supplied times
are authoritative, including for named events.

## Catalogue

Names are case-insensitive and ignore spaces, hyphens and punctuation.
Tōhoku is accepted with its macron. These are fixed documented epicentres;
the function does not query a changing catalogue during a fit.

| Normalised name | Aliases | Longitude | Latitude | Expected rounded date |
| --- | --- | ---: | ---: | --- |
| Sumatra | Sumatra2004, Sumatra-Andaman, Indian Ocean | 95.982 | 3.295 | 2004-12-26 |
| Tohoku | Tohoku2011, Japan | 142.373 | 38.297 | 2011-03-11 |
| Chile | Chile2010, Maule | -72.733 | -35.909 | 2010-02-27 |

Sources: [USGS Sumatra event](https://earthquake.usgs.gov/earthquakes/eventpage/official20041226005853450_30),
[USGS Tohoku event](https://earthquake.usgs.gov/earthquakes/eventpage/official20110311054624120_30),
and the [USGS 2010 Chile report](https://pubs.usgs.gov/of/2010/1152/presentations/of2010-1152_20100720i_kanamori.pdf).

Numeric coordinates matching a catalogue entry within `1e-10` degrees use
its name and coordinates. Longitudes outside `[-180, 180)` are wrapped;
longitude at a pole becomes zero. Other coordinates remain arbitrary
events and require explicit dates.

## Fitting and uncertainty

The polar-cap basis is sorted by decreasing concentration and rotated to
the epicentre using the Slepian package's harmonic convention. The default
retained count is `min((L+1)^2, max(1, round(N)))`, as in
`prepareslepianbasis`. This is explicitly `round(N)`, not `ceil(N)`.
Projection, fitting and synthesis use exactly those retained columns.

Each call to `slept2seismic` receives only one event. Ordinary least squares
matches that wrapper's default; GRACE formal coefficient uncertainties are
not used as temporal weights. The jump's standard errors come from the
joint temporal fit covariance and its residual variance. With retained
basis `G` and jump standard errors `s`, harmonic variance is
`(G.^2)*(s.^2)`. Slepian modes and events are assumed independent even when
they share observations; GLDAS is a fixed correction with no model-error
term. No additional physical-unit or basis-normalisation factor is applied.

`graceepochs` provides the solution windows, synthetic-window flags and
optional source reference radius. `grace2seismic` excludes the synthetic
windows before fitting or assembling the correction.
The loader's midpoint dates must match every window. Missing months and
the mission gap are not filled. `TimeRange` selects only returned midpoints
and does not restrict a fit or change an event cache. For window `[a,b]`
and rounded event day `t`, the output weight is
`max(0,b-max(a,t))/(b-a)`. Weighted event amplitudes add, and their variances
add with squared weights. The result is zero before an event, fractional
across it, and constant afterwards. The temporal fit itself uses steps at
the observation midpoints.

Both `timefirst` and `traditional` layouts retain degree/order labels in
the field and uncertainty arrays. Units follow `grace2plmt_new`; the default
is `SD` in kg/m², and the default returned time format is MATLAB datenum.
RL05 GSM, RL06 GSM and CSR mascon source discovery follows the loader.
`L` cannot exceed the product's bandwidth, and every selected harmonic must
be available at every observation.

## NetCDF event caches

The directory is `GRACEDATA/Seismic`, falling back to
`IFILES/GRACE/Seismic`. The stable filename is:

```text
grace2seismic-<centre><releaseDigits><productDegree>-<eventId>_<radius>[-TYYYYMMDD]-L<degree>-<unit>[-GLDAS].nc
```

Examples:

```text
grace2seismic-CSR0660-Sumatra_20-L60-SD.nc
grace2seismic-CSR0660-Sumatra_20-T20041227-L60-SD.nc
grace2seismic-CSR0660-121p48_n30p72_20-T20091028-L60-SD-GLDAS.nc
```

`-TYYYYMMDD` is omitted only for a resolved catalogue event whose rounded
input date is its expected day. `-L` is always present; `-GLDAS` appears only
when removal is enabled. `num2filename` chooses sufficient significant
digits to round-trip a double, replaces decimal points with `p` and minus
signs with `n`, and keeps ordinary coordinate/radius magnitudes in decimal
notation. There are no hashes or cache versions.

Files use `writeshnc`/`readshnc`: double `coefficients` and `uncertainty`
arrays in MATLAB `(lm, cs, time)` order, separate `degree`, `order`, `cs`
and `time` coordinates, and a single time entry containing the event day.
`uncertainty` contains standard deviations. Global attributes record the
epicentre, cap radius, actual retained mode count, Shannon number, fit
background and date extent, Love-number source and GLDAS treatment.

Completeness checks require those six variables; attributes are provenance,
not cache validation. Missing variables trigger regeneration to the same
path when saving is enabled. Users clear caches after changing the source
record or scientific defaults. No combined-event cache is written.

`ForceNew=true` bypasses event and processed GRACE caches. `SaveData=false`
prevents event writes and is forwarded to upstream routines with save
controls. The legacy `plm2rot` dependency `dlmb` has no save control and may
write its reusable rotation matrix under `IFILES/DLMB`; that directory must
exist when the matrix is first computed. Geometric caches are not forced
by an event refresh.

Other fit parameters and corrections have no filename identifiers. Changes
to `truncation`, `polynomialDegree`, `periods`, `FitRange`, Love numbers or
low-degree corrections therefore require both `ForceNew=true` and
`SaveData=false`. The function rejects conflicting cache settings:

```matlab
correction = grace2seismic(product, "Sumatra", 20, ...
    datetime(2004, 12, 26), 60, polynomialDegree=1, periods=[], ...
    ForceNew=true, SaveData=false);
```

## GLDAS subtraction

`RemoveGLDAS=true` uses only **GLDAS_NOAH10_M, release 2.1**, the final
monthly 1-degree Noah product. Place its original files
`GLDAS_NOAH10_M.AYYYYMM.021.nc4` under the environment-configured
`GLDASDATA` directory, optionally in year subdirectories. The fallback is
`IFILES/EARTHMODELS/GLDAS2.1`. The function reads local files and does not
download or generate a second GLDAS cache.

Storage is the sum of four soil layers (0–10, 10–40, 40–100 and 100–200 cm),
canopy interception and snow water equivalent, all in kg/m². There is no
groundwater term. NASA documents the [Noah layers](https://ldas.gsfc.nasa.gov/gldas/specifications)
and [storage variables](https://ldas.gsfc.nasa.gov/gldas/model-output).
The [product archive](https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH10_M.2.1/)
provides the monthly files. Other GLDAS models/releases and early products
are not substituted.

The complete native grid covers all longitudes and latitudes north of
60°S. The six components must share a time-invariant native missing-data
mask; those masked cells and latitudes south of 60°S contribute zero.
No additional glacier mask is applied. Missing variables, incompatible
units or grids, changing component masks, duplicate months and incomplete
temporal coverage raise errors.

Monthly values are treated as constant within their calendar month and
averaged by overlap duration over each actual GRACE window. No temporal
extrapolation or midpoint-only approximation is used. Cell-constant storage
is integrated into the standard real harmonic basis with exact longitude
integrals and eight-point latitude quadrature. This supports `L < 180` for
the native 1-degree grid. Degree zero is removed to match GRACE's zero C00;
other low degrees are retained. Unit conversion uses the GRACE reference
radius and requested Love numbers. These harmonics are subtracted before
projection and fitting, which is equivalent to subtracting their projection
in exactly the same cap basis.

The older `gldas2TWSt` helper is not used: it assumes a different product
layout, fixes degree 60 and unconditionally saves generated data.
