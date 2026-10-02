# Jason 3 correction implementation plan

DT2024 omits the supplementary Jason-3 radiometer drift adjustment. Add an
explicit, opt-in correction to `ssh2lonlatt`, implemented by a separate temporal
helper. Original data and all existing interpolated SSH caches remain
uncorrected. This PR contains the plan only; it does not change numerical
results, add the helper, or enable any correction.

## Agreed interface and storage behavior

Proposed public interface:

```matlab
[correction, dates] = sshcorrection2t('Jason3', timestep, timelim, ...
    Interpolation='linear', Unit='m');

[ssh, sshSigma, dates, lon, lat] = ssh2lonlatt('DT2024', ...
    ApplyJason3Correction=true);
```

`ApplyJason3Correction` defaults to `false`. The helper follows existing ULMO
conventions for timestep, time range, interpolation methods, units, date format,
verbosity and call-chain reporting. Its output is a column vector of range
corrections to subtract from SSH, not a replicated spatial raster. Provide an
explicit-epoch route, provisionally `QueryDates`, so the loader can request its
actual epochs; verify exact date equality instead of independently generating
and assuming matching grids. Conflicting epoch/timestep requests must be rejected.

Use the existing correction text file under `IFILES/SSH` and document its release
date. No helper disk cache is necessary. Add no cache version numbers, hashes,
new cache naming scheme, or additional invalidation framework. Existing native
and interpolated cache paths and contents retain their current meaning.
`SaveData=true` saves uncorrected SSH exactly as before; the correction modifies
only the returned or plotted array. Repeated corrected calls must never accumulate
the adjustment. Existing saved downstream analysis results require explicit
recomputation when the caller changes the option; this PR does not redesign
their caching.

## Focused implementation commits

### 1. Define the DT2024 correction recipe

- Document the JPL file's four columns, centimeter units, release date and valid
  time coverage. Define positive output as a correction subtracted from SSH;
  convert centimeters to meters by multiplying by 0.01.
- Specify the global, spatially uniform approximation for merged monthly L4
  data separately from the exact pass-wise L2 correction. The ESA approach
  interpolates onto tracks and uses spatial averaging before global averaging;
  a simple mean of all pass values must not be presented as the same algorithm.
- Establish whether the local text file is sufficient for the selected method
  or whether a provider-derived global series/track coordinates are needed.
  Record that choice and its limitations before wiring it into the loader.
- Resolve product-specific start/handover dates, offset continuity through the
  transition to Sentinel-6, partial months, and behavior outside source coverage.
  Do not silently turn the correction off at handover, extrapolate indefinitely,
  or use the September audit's sensitivity cases as a validated recipe.
- Define calendar-month means first for monthly DT2024 data, then any requested
  temporal interpolation. Existing nearest-midpoint downsampling in
  `interptemporal` must not accidentally substitute for calendar-month averaging.

Review checkpoint: the scientific recipe is explicit and reproducible, including
any remaining approximation. Production integration depends on resolving it.

### 2. Implement the temporal helper with focused tests

- Add `sshcorrection2t.m` and a small parser/averaging helper if needed.
- Support the loader's interpolation choices, `Unit`, `TimeFormat`, `BeQuiet`,
  `CallChain`, timestep/time-range conventions and exact requested epochs.
- Reject malformed files, invalid or duplicate timestamps, unsupported correction
  names, conflicting requests and unsupported dates with clear ULMO errors.
- Keep monthly integration distinct from point interpolation; ensure requesting
  a shorter interval does not change endpoint monthly averages by discarding
  source support too early.
- Add small synthetic fixtures testing sign and unit conversion, nonconstant
  calendar-month averages, leap February, interpolation, actual-epoch alignment,
  and the transition/coverage rules from commit 1. Tests require no downloads or
  user-specific `IFILES` archive.

Review checkpoint: helper behavior can be assessed independently of SSH rasters.

### 3. Unify loader completion paths without changing results

- Refactor the cache-hit and recomputation branches to converge after loading
  or computing the uncorrected interpolated arrays, before formatting/plotting.
- Preserve the existing save location, filenames, variable names, default output,
  `ForceNew`, `SaveData`, units, output orientation and date filtering semantics.
- Preserve source dates and coordinates for uncertainty interpolation; the
  current uncached path reuses transformed axes and includes a stray positional
  argument. Include the minimal fix if exercised by this refactor, with a
  separate regression case; do not expand into general interpolation changes.
- Use temporary synthetic native files to prove cached and fresh default calls
  agree, including uncertainty arrays, missing cells and the plotting path.

Review checkpoint: correction-disabled behavior remains compatible with the
existing loader and caches; any necessary bug fix is individually identified.

### 4. Integrate the opt-in correction and cache-safety tests

- Add `ApplyJason3Correction=false` and initially allow it only for the identified
  DT2024 product. Reject other products, including already-corrected NASA-SSH
  when available; do not depend on unmerged NASA product-enum additions.
- After the existing uncorrected cache save/load step, obtain the correction
  aligned to the actual output epochs and subtract it once from finite SSH cells.
  Apply it in meters before output-unit conversion. Preserve NaNs and the meaning
  of `sshSigma`; do not add a systematic drift uncertainty as independent cell noise.
- Ensure plots use the corrected output when requested, just as returned arrays do.
- Test warm/cold caches, `ForceNew`, `SaveData`, on/off/on calls, both output
  orientations, units and date formats. Assert source/cache bytes are unchanged
  on cache hits and a cold saved cache contains uncorrected values.

Review checkpoint: the option changes only the intended in-memory SSH values,
with no cache contamination or double application.

### 5. Document usage and validate against the real data

- Update function help and the relevant ULMO documentation with the correction's
  meaning, supported product, source file, averaging and transition rules,
  out-of-range behavior, and uncorrected-cache semantics.
- Document global-approximation limits and shared systematic uncertainty; distinguish
  formal grid errors from uncertainty in the correction model.
- Run the focused helper/loader tests and relevant existing temporal-interpolation
  checks. Perform a read-only real-data comparison over 2003–2022, checking both
  return paths and verifying that corrected-minus-original equals the negative
  correction at finite cells. Compute the actual trend effect from the chosen
  recipe; do not force agreement with the prior sensitivity estimates.
- Record the validation results in the PR. Do not commit downloaded datasets,
  production MAT files, credentials, or changes to analysis defaults.

Review checkpoint: the documented method, numerical effect and cache guarantees
are demonstrated before the draft is marked ready for review.

## Sources

- [JPL supplementary correction](https://podaac.jpl.nasa.gov/dataset/JASON_3_PD_CORRECTION)
- [DT2024 paper and author response on the omitted correction](https://essd.copernicus.org/preprints/essd-2025-604/)
- [ESA sea-level budget algorithm, section 2.1.3](https://climate.esa.int/media/documents/2024-10-31_SLBC_CCI-DT-041-MAG_ATBD_D2-4_V2_0_signed.pdf)
- [AVISO mean-sea-level post-processing](https://www.aviso.altimetry.fr/en/data/products/ocean-indicators-products/mean-sea-level/processing-and-corrections.html)
