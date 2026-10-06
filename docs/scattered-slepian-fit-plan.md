# Scattered observations to regional Slepian coefficients

Status: proposed implementation plan, awaiting review. This PR adds documentation only; implementation starts after approval.

## Purpose and naming

Add `xyzs2slep` to fit scalar observations at paired, irregular longitude–latitude locations directly to a truncated regional Slepian basis, without gridding. Support arbitrary geographic concentration regions through ULMO's existing domain machinery and propagate observation uncertainty to coefficient covariance.

Recommend interpreting the `s` as **scattered**, not sparse: it describes the sampling and avoids confusion with MATLAB sparse matrix storage. `points2slep` is a clearer alternative but is less consistent with the existing `xyz` transform family. Use `xyzs2slep` provisionally for this plan, subject to review.

Reserve `xyz2slep` for a future ULMO gridded-data interface. ULMO currently has no such wrapper: Slepian Bravo's existing `xyz2slep` handles scattered cap data. This change will not rename or shadow that dependency, so a bare `xyz2slep` call can still resolve to Bravo. Document that distinction; implementing or changing the gridded interface is a separate task.

## Proposed public contract

Stage 1 delivers ordinary least-squares fitting without uncertainty inputs or propagation:

```matlab
[falpha, V, N] = xyzs2slep(values, lon, lat, domain, L, truncation=J);

[falpha, V, N, J, nObs, usedIndices, numericalRank, singularValues, ...
    conditionNumber, fittedValues, residuals, residualDof] = ...
    xyzs2slep(values, lon, lat, domain, L, truncation=J);
```

Stage 2 adds uncertainty options and appends uncertainty outputs without changing the Stage 1 output positions:

```matlab
[falpha, V, N, J, nObs, usedIndices, numericalRank, singularValues, ...
    conditionNumber, fittedValues, residuals, residualDof, ...
    coefficientCovariance, coefficientStd, uncertaintySource] = ...
    xyzs2slep(values, lon, lat, domain, L, truncation=J, dataStd=sigma);
% Alternatively use dataCovariance=Cdata instead of dataStd=sigma.
```

Each result is a separate output argument, consistent with the Slepian packages; there is no `info` structure. Callers request only the leading outputs they need or use `~` to skip outputs. The common coefficient/eigenvalue/Shannon-number outputs follow `plm2slep_new` ordering.

- `values`, `lon`, and `lat`: equal-length real vectors, normalized to columns; one scalar field per call initially. Coordinates are paired points, never a Cartesian-product grid. Public coordinates use longitude and latitude in degrees, with explicit conversion to colatitude in radians internally.
- `domain`: a `GeoDomain`, named region, buffered-domain cell specification, or closed polygon in `[lon, lat]` degrees, following `glmalpha_new`. Numeric scalar caps may use the same existing machinery. Complex ocean geometry must retain the topology supported by `GeoDomain`; do not replace it with a bounding cap or silently invent polygon semantics.
- `L`: initially a nonnegative integer maximum degree. Bandpass and multiple epochs are deferred to keep the first implementation focused.
- `truncation`: positive integer number of retained functions; default `max(1, round(N))`, bounded by the available basis dimension. Report the actual choice. Require at least as many usable observations as coefficients and full numerical column rank; do not silently reduce the requested truncation to fit the data.
- `dataStd` (Stage 2): optional vector of finite, strictly positive observation standard deviations, in the same units as `values`.
- `dataCovariance` (Stage 2): optional finite, symmetric positive-definite observation covariance, in squared data units. Mutually exclusive with `dataStd`. Validate dimensions and symmetry with a documented numerical tolerance, then use Cholesky factorization; do not silently add jitter or discard correlations. Singular covariance and exact constraints are deferred.
- `rankTolerance`: optional relative singular-value threshold, with a documented dimension-aware numerical default. This diagnoses identifiability; it is not an implicit regularization parameter.
- `estimateNoiseVariance` (Stage 2): default false. Without supplied uncertainty, optionally estimate a common iid noise variance from residuals when there are positive residual degrees of freedom. Reject this option together with supplied absolute uncertainty.
- `missingPolicy`: default `"error"`; optional `"omit"` removes rows with missing values/coordinates and, in Stage 2, consistently subsets standard deviations or both covariance axes. Invalid covariance entries and nonpositive standard deviations remain errors. Report retained indices.

`falpha` is a J-by-1 vector in descending concentration-eigenvalue order, with normalization compatible with `plm2slep_new` and `slep2plm_new`. `coefficientCovariance` is J-by-J, in squared coefficient units; return `[]` when uncertainty is unspecified and residual-based estimation was not requested. Never interpret missing uncertainty as zero uncertainty.

The output contract is:

| Output | Meaning |
| --- | --- |
| `falpha` | J-by-1 fitted coefficients. |
| `V` | J-by-1 concentration eigenvalues for the retained basis, in matching order. |
| `N` | Shannon number of the concentration problem. |
| `J` | Actual truncation. |
| `nObs` | Number of usable observations. |
| `usedIndices` | Column vector of retained original observation indices. |
| `numericalRank` | Numerical column rank of the fitted design matrix. |
| `singularValues` | Singular values of A in Stage 1, or the whitened matrix B when weighted fitting is requested in Stage 2. |
| `conditionNumber` | Ratio of largest to smallest retained singular value of that same matrix. |
| `fittedValues` | Model values at retained observations, in data units. |
| `residuals` | Observed minus fitted values at retained observations, in data units, even for a weighted fit. |
| `residualDof` | Residual degrees of freedom, `nObs - J` for the required full-rank fit. |
| `coefficientCovariance` | Stage 2: J-by-J covariance, or `[]` when uncertainty is unavailable. |
| `coefficientStd` | Stage 2: J-by-1 square roots of covariance diagonal entries, or `[]`. |
| `uncertaintySource` | Stage 2: `"dataStd"`, `"dataCovariance"`, `"estimatedIid"`, or `"none"`. |

Avoid returning the full design matrix by default. Compute optional output-only quantities only when requested, while always performing validation and rank checks needed for the fit.

## Basis and spatial evaluation

1. Obtain the geographic basis from `glmalpha_new`, reusing its domain handling and cache controls. Keep basis computation out of production data paths in tests.
2. Sort basis columns and concentration eigenvalues together, then retain the requested J functions. Verify that any partial-basis request returns the globally intended functions before using that optimization.
3. Construct the M-by-J design matrix `A(i, alpha) = g_alpha(lon(i), lat(i))` at the paired observation points. With a consistent spherical-harmonic convention, this is `A = Y * GJ`, where Y is the M-by-K harmonic evaluation matrix and GJ is K-by-J.
4. Audit coefficient ordering, real-harmonic phase, and normalization against ULMO's synthesis functions before choosing the evaluation helper. In particular, `ylm` and `plm2xyz` document different phase conventions, and `plm2xyz`'s scattered-coordinate documentation and implementation need verification. Do not copy a `Y * GJ` expression without a synthesis-equivalence test.
5. Evaluate points in blocks to avoid materializing all M-by-K harmonics at once. The initial dense solve still stores M-by-J values and, if supplied, an M-by-M covariance. Document these memory limits.

The concentration domain defines the basis, independently of the observation coordinates. Use all explicitly supplied valid observations, including those outside the concentration region; do not silently clip them. Document this so callers can subset their observations intentionally.

## Stage 1: standalone fitting without uncertainty propagation

Construct A as above and solve `min ||A*c - d||^2` using an economy SVD. For `A = U*S*Q'`, compute `c_hat = Q * diag(1 ./ s) * U' * d`. Use Q here to distinguish the right singular vectors from the public concentration-eigenvalue output V. Return the first twelve outputs specified above. This stage includes arbitrary-domain support, paired-point evaluation, truncation, missing-data handling, rank/conditioning diagnostics, documentation, and its own tests. It does not accept uncertainty options or calculate coefficient covariance, standard deviations, or residual-based noise variance.

Fail clearly for insufficient observations or numerical rank deficiency, and document a warning criterion for poorly conditioned full-rank fits. An exactly determined full-rank fit is supported. Stage 1 must be complete, tested, and independently reviewable before starting Stage 2.

## Stage 2: uncertainty propagation and weighted fitting

Extend the tested Stage 1 function with the uncertainty options and three appended outputs. Preserve the Stage 1 call behavior, coefficient results, and output order when uncertainty options are absent. Reuse the basis evaluation and solver rather than introduce a second inconsistent fitting path.

The model is `d = A*c + epsilon`. With no uncertainty, solve ordinary least squares. With standard deviations, whiten by dividing each row of A and d by its standard deviation. With full covariance `Cdata = R*R'` (lower Cholesky factor), whiten using `B = R\A` and `b = R\d`.

Use an economy SVD of B (or A for unweighted fitting) to solve stably without explicitly inverting covariance or forming normal equations. For `B = U*S*Q'` and full column rank:

```text
c_hat = Q * diag(1 ./ s) * U' * b
Ccoef = Q * diag(1 ./ s.^2) * Q'
```

The covariance expression applies when the supplied observation covariance is absolute and correctly specified. It is equivalent to `(A' * Cdata^-1 * A)^-1` but must be evaluated through the factorization. More generally, for the linear estimator H, propagation is `Ccoef = H * Cdata * H'`.

When no uncertainty is supplied and `estimateNoiseVariance=true`, use `sigma2 = sum(residual.^2)/(M-J)` and multiply the unweighted covariance factor by sigma2. Require `M > J`. Label this as an iid residual-based estimate, which can absorb model mismatch, rather than propagated measurement uncertainty. Do not rescale covariance supplied by the caller using residual variance.

For rank-deficient geometry, fail with an informative identifier and suggest fewer retained functions or improved spatial coverage. A minimum-norm pseudoinverse can otherwise conceal unconstrained coefficient directions behind misleading finite variances. Report conditioning for full-rank but poorly constrained fits and define a documented warning criterion during implementation.

This uncertainty is conditional on the selected basis, truncation, and observation error model. It excludes truncation bias, leakage, coordinate uncertainty, and uncertainty in the region/basis. Irregular sampling alone does not supply area weights: equal-error observations receive equal weight, and densely sampled areas have more influence. Area weighting and robust fitting are deferred.

## Implementation sequence after review

1. Finalize the function name and public contract in this plan.
2. **Implement standalone fitting (Stage 1).** Add `xyzs2slep.m` with an `arguments (Input)` block, camelCase name-value options, and ULMO-style documentation. Implement basis evaluation, ordinary least squares, validation, and the twelve individual fitting outputs. Add a narrowly scoped `aux/` helper only if reusable paired-point evaluation needs it. No uncertainty calculation is part of this step.
3. **Validate and review Stage 1 independently.** Add native MATLAB tests in `unittests/xyzs2slepTest.m`, a reproducible no-uncertainty example, and a `docs/functions.md` entry. Run fitting tests and relevant basis/transform tests. Deliver this as a self-contained implementation PR before Stage 2; it must be usable without uncertainty support.
4. **Add uncertainty support (Stage 2, separate follow-up PR).** Extend the established fit with diagonal/full-covariance whitening, coefficient covariance propagation, optional residual-variance estimation, and the three appended uncertainty outputs. Do not reorder the existing outputs.
5. **Validate Stage 2 separately.** Add uncertainty-specific analytic and Monte Carlo tests, uncertainty documentation/examples, and regression checks proving that calls without uncertainty still match Stage 1. Record validation and dependency limitations in the follow-up PR.

Use temporary cache fixtures and synthetic data in both stages. Keep production data untouched and machine-specific paths out of the repository. This planning PR remains documentation-only and awaiting approval.

## Acceptance tests

### Stage 1: fitting only

- Recover known coefficients from noiseless irregular points for an asymmetric geographic polygon, and exercise a `GeoDomain` with supported ocean topology plus a cap regression case. Include points on both sides of the longitude seam.
- Check basis evaluation against an independent existing synthesis route at low degree, including nonzonal cosine and sine terms to expose phase/order mistakes. Test compatibility with ULMO transforms rather than only recovering data generated by the new helper.
- For a well-conditioned cap, compare fitted spatial values with Bravo `xyz2slep` after reconciling coordinate and normalization conventions; compare coefficients only with an identical basis, since signs and degenerate eigenspaces can differ.
- Match ordinary least squares on a small known system; exercise duplicate/clustered points, exact rank deficiency, insufficient observations, and an exactly determined full-rank system.
- Check individual output ordering, shapes and meanings, vector orientation, invalid coordinates, missing-data omission, and out-of-domain observations retained as documented.
- Verify concentration ordering and block-size-independent evaluation without requiring exact eigenvector signs across independent eigensolves.

### Stage 2: uncertainty and regression coverage

- Match independently whitened generalized least squares on small known systems; show that constant standard deviation leaves coefficients unchanged and diagonal covariance matches `dataStd`.
- Match coefficient covariance to an independent analytic small-system result. A fixed-seed Monte Carlo test checks empirical covariance against propagated covariance, including correlated observation errors, with statistically justified tolerances.
- Verify scaling observation standard deviations by k preserves the weighted fit and scales covariance by k^2; scaling data and uncertainty together scales coefficients by k and covariance by k^2.
- Confirm an exactly determined full-rank system accepts known measurement covariance but rejects residual-variance estimation. Test residual-based noise estimation separately from supplied absolute uncertainty.
- Check missing-data subsetting of both covariance axes, incompatible options, nonsymmetric/non-positive-definite covariance, appended output shapes, and uncertainty provenance.
- Re-run all Stage 1 tests; verify calls without uncertainty preserve the established output positions and coefficient results, with empty covariance/standard-deviation outputs and `uncertaintySource="none"` when requested.

## Review decisions

The proposed defaults are `xyzs2slep` (s = scattered), degree-based lon/lat inputs, full covariance support, strict rank checking, and no inferred uncertainty unless explicitly requested. The plan now uses individual outputs and independently deliverable fitting and uncertainty stages. Please review the revised API, truncation default, and missing/out-of-domain behavior before implementation.

No executable MATLAB changes, package installation, or production data processing are included in this planning PR.
