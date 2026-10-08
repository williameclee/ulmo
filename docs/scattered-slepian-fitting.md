# Fitting scattered observations

`xyzs2slep` fits a scalar field observed at paired longitude/latitude points directly to a truncated regional Slepian basis. No gridding is performed, and each observation receives equal weight.

```matlab
% Inputs are vectors of equal length. Coordinates use degrees.
[falpha, V, N] = xyzs2slep(values, lon, lat, domain, L, truncation=J);
```

Use a `GeoDomain`, a region name, a buffered-region cell specification, or a closed `[lon, lat]` polygon for `domain`. Complex regions retain the geometry supported by `glmalpha_new`; its integration limitations still apply. A scalar domain is a north-polar cap radius in degrees, evaluated through Slepian Alpha's `glmalpha` because ULMO's extracted scalar-cap helper is incomplete.

`domain` defines where the basis concentrates. It does not filter observations: explicitly supplied points outside the region participate in the fit. Longitudes are periodic, and latitudes must lie between -90 and 90 degrees. NaNs raise an error unless `missingPolicy="omit"` is specified; infinite inputs always raise an error.

## Reproducible example

With ULMO and Slepian Alpha on the MATLAB path and `IFILES` configured for the existing basis caches:

```matlab
% An asymmetric polygon straddling longitude zero.
domain = [340 -25; 380 -20; 375 30; 345 20; 340 -25];
L = 3;
J = 5;
truth = (1:J)'/10;

% Deterministic scattered points, including points outside the region.
k = (1:80)';
lon = mod(round(k*137.507764), 360);
lat = round(asind(-0.95 + 1.9*(k-0.5)/80));

% Generate observations using the existing ULMO synthesis path.
lmcosi = slep2plm_new(truth, domain, L, 0, 0, 0, J);
grid = plm2xyz(lmcosi, 1, 'BeQuiet', true);
values = grid(sub2ind(size(grid), 91-lat, lon+1));

[falpha, V, N, usedJ, nObs, usedIndices, numericalRank, singularValues, ...
    conditionNumber, fittedValues, residuals, residualDof] = ...
    xyzs2slep(values, lon, lat, domain, L, truncation=J);
maxCoefficientError = max(abs(falpha-truth))
```

The example selects irregular paired samples from a synthesised grid only to create independently verifiable test data. ULMO's current `plm2xyz` input parser rejects vector coordinates despite its documented scattered mode. The new fitter does not call `plm2xyz` or grid observations: real inputs may have arbitrary noninteger coordinates.

## Outputs and numerical behaviour

Outputs are separate positional arguments, beginning with coefficients, concentration eigenvalues, and Shannon number, following the Slepian transform convention. `help xyzs2slep` lists all twelve. The coefficients match ULMO's $4\pi$-normalised harmonic transforms; the implementation explicitly corrects `ylm` phase and normalisation. They should not be compared directly to Bravo's unit-normalised cap coefficients without converting conventions and ensuring the same basis.

The default truncation is `max(1, round(N))`, limited to the basis dimension. Specify `truncation` to choose it explicitly. The sampled design matrix must have full numerical column rank, with at least as many observations as retained coefficients. Invalid truncations and rank-deficient fits raise errors instead of silently discarding modes. Exactly determined full-rank fits are supported.

An economy SVD solves ordinary least squares without normal equations. `rankTolerance` is a relative singular-value cutoff (default `max(nObs,J)*eps`). Full-rank fits with condition number above `1/sqrt(eps)` warn. A low concentration eigenvalue and a small sampled-design singular value describe different limitations; spatial sampling can constrain a concentrated basis poorly.

Basis construction requests the full basis before sorting by concentration. Harmonic evaluation uses point blocks controlled by `blockSize` (default 1024), but the dense observation-by-truncation design matrix and its SVD still consume memory. The full basis itself scales as `(L+1)^4` elements. Existing basis/kernel routines may read and write their normal caches under `IFILES`.

This version has no uncertainty propagation, weighted fitting, area weighting, or regularisation. Those belong to the next implementation stage. Equal observation weights give densely sampled regions more influence. `xyz2slep` is not changed or shadowed: with Bravo installed it still resolves to Bravo's scattered-cap routine; a future gridded ULMO interface is separate work.
