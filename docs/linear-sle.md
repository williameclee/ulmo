# Linear sea-level and low-degree solvers

`solvesle` and `solvedegree1` compute new results with the linear solver. Both use the same
elastic Earth response, fixed ocean geometry and mass constraint. Internal
coefficients have kernel order, fields are columns, and loads use kg/m².

## Shared matrices

```matlab
L = 96;
ocean = GeoDomain('alloceans', Buffer=0.5);
oceanKernel = kernelcp_new(L, ocean, BeQuiet=true);
[earthResponseKernel, forcingKernel, waterDensity] = ...
    prepareslekernels(L, oceanKernel, Frame="CF", ...
    LoveNumSource="ISSM", RotationFeedback=true);
```

The first two outputs are ordinary sparse numeric matrices: `earthResponseKernel`
is `G+R`, and `forcingKernel` is the mass-constrained forcing matrix `KF`.
The optional third output, `waterDensity`, is the scalar conversion factor
from metres of sea level to kg/m². The ocean matrix stays with the caller;
its first column is the ocean function and its `(1,1)` entry is the ocean
area fraction. No kernel struct or model-metadata output is returned.

Writing `O = oceanKernel` and `e` for the constant-harmonic unit vector, the system is

```text
P  = I - e*O(1,:)/O(1,1)
KF = P*(G+R) - e*e.'/O(1,1)
A  = I - P*(G+R)*O
A*S = KF*F
```

`solveslesystem` uses products of these matrices in GMRES, without allocating
another dense matrix for `A`. This is a concrete memory saving even when the
ocean matrix is stored explicitly. Preparation returns no function handles.
The finite-bandwidth ocean matrix is not treated as an idempotent projection.

```matlab
% F is an M-by-N matrix of forcing coefficients, M=(L+1)^2.
[solutionCoefficients, relativeResiduals, iterationCounts] = ...
    solveslesystem(oceanKernel, earthResponseKernel, forcingKernel*F);
```

`relativeResiduals` and `iterationCounts` are separate numeric row vectors,
with one entry per right-hand side. Both entries are zero for a zero
right-hand side. The same output convention applies to `Transpose=true`.

The public `solvesle` function handles lmcosi/cosi packing, units and formatting.
It does not localize the supplied forcing to land.

## Selected low-degree estimates

```matlab
[~, unknownIds] = lowdegreecoefficients(false, false, true);
[knownToUnknownMap, knownIds] = preparelowdegreemap( ...
    oceanKernel, earthResponseKernel, forcingKernel, unknownIds);
unknownLoads = knownToUnknownMap*knownLoads(knownIds,:);
unknownStds = sqrt((knownToUnknownMap.^2)*knownVariances(knownIds,:));
```

Here `knownLoads` and `knownVariances` are full-bandwidth packed arrays; the
unknown rows are ignored. GIA removal and restoration remain the caller's
responsibility when using this low-level interface.

The map builder solves `H*A = O(unknownIds,:)` using only one transposed system
per selected coefficient, then forms `oceanResponseMap = H*KF*(I-O)` and

```matlab
unknownConstraint = oceanResponseMap(:,unknownIds) - oceanKernel(unknownIds,unknownIds);
knownConstraint = oceanKernel(unknownIds,knownIds) - oceanResponseMap(:,knownIds);
knownToUnknownMap = unknownConstraint \ knownConstraint;
```

The matrices, `knownIds`, reciprocal condition number, relative residuals and
iteration counts are returned as separate outputs. The first two outputs are
sufficient for applying the map; additional outputs expose the ocean-response
and constraint matrices for inspection.

`solvedegree1` performs these steps once for all epochs and handles GRACE
loading, GIA, coefficient selection, units and dates. It never
calls `solvesle`. `Method="uniform"` builds the analytical uniform-ocean map
without an SLE solve; its low-level call is
`preparelowdegreemap(oceanKernel, [], [], unknownIds, Method="uniform")`.
Both methods support joint C20, C30 and C21/S21 estimates.

`preparelowdegreemap` builds a map on every call and performs no file I/O.
The returned map can be reused in memory while the kernels, coefficient
selection and numerical settings remain fixed.

`solvedegree1` retains its coefficient-result cache. `ForceNew=true` bypasses
that cache and prepares a fresh map; `SaveData=false` prevents new coefficient
result writes. Coefficient results use the original filenames, with no solver,
version or accuracy suffix. Existing cached results are loaded as stored,
including results from the previous solver. The caller decides when to use
`ForceNew`, including after solver, numerical-setting, observation or model changes.

## Uncertainty and numerical settings

The linear solver propagates independent input variances through the **complete**
response. It retains correlations induced by localization, ocean feedback and
mass conservation. Its standard deviations therefore differ from the legacy
intermediate diagonal approximation. GIA/model uncertainty is not included.

For full SLE fields, `slevariance` computes response columns in blocks and omits
exactly zero input variances. This bounds temporary memory but may still require
one large solve for every uncertain input harmonic. Low-degree uncertainty uses
only the small prepared `knownToUnknownMap` and is inexpensive.

`RelativeTolerance` controls the relative Euclidean residual of each linear
solve; its default is `1e-12`. `MaxIter` in `solvesle`, and `SleMaxIter` in
`solvedegree1`, cap the GMRES iterations per right-hand side. Unconverged solves
raise an error. `Tolerance` additionally checks ocean mass balance in metres for
`solvesle`, and the low-degree constraint in kg/m² for `solvedegree1`.

The public functions no longer accept `Solver` or zero `Tolerance`. Legacy-only
options (`InitialRSL` and `OceanVarianceKernel` in `solvesle`; `MaxIter` and
`SleTolerance` in `solvedegree1`) have also been removed.

For regression testing, the previous implementations live in
`legacy/solvesle_iteration.m` and `legacy/solvedegree1_iteration.m`.
Their diagonal-variance helper is `legacy/prepareoceanvariancekernel.m`.
The suffix identifies the solvers’ fixed-point iteration method. These references retain
fixed-step execution and the previous variance recurrences. Tests add `legacy/`
to the MATLAB path through a path fixture; the production functions do not call
these references.
