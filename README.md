# Unified Localisation of Mass over Oceans (ULMO)

A pipeline for using GRACE satellite gravimetry and altimetry-steric (in development) measurements to study ocean mass change over ocean basins, with extended functionalities for Slepian functions and sea level equation modelling.

![Example of Slepian functions over the Pacific Ocean](images/cover.svg)

## Dependencies

Functions in this repository may call or overwrite functions from the following packages. Please ensure they are installed and the paths are configured properly before running the functions in this repository.

- [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha.git): Computation of Slepian functions on the sphere
- [slepian_bravo](https://github.com/csdms-contrib/slepian_bravo.git): Conversion between spherical harmonics and Slepian functions
- [slepian_delta](https://github.com/csdms-contrib/slepian_delta.git): Processing GRACE data and GIA models
- [Gibbs Seawater Toolbox (GSW)](https://www.teos-10.org/software.htm): Computation of steric sea level from temperature and salinity profiles
<!-- - [MatlabColourmapGenerator](https://github.com/williameclee/MatlabColourmapGenerator) -->

<!-- [MatlabColourmapGenerator](https://github.com/williameclee/MatlabColourmapGenerator) is not necessary, but it generates better-looking colours and colourmaps (like the one above). -->

### Geographic domains

The following ocean basins are supported:

- `alloceans`: All oceans, including everything
- `oceans`: All oceans, excluding the Arctic Ocean at the moment
- `pacific`: The Pacific Ocean
- `npacific`: The North Pacific Ocean
- `spacific`: The South Pacific Ocean
- `atlantic`: The Atlantic Ocean
- `natlantic`: The North Atlantic Ocean
- `satlantic`: The South Atlantic Ocean
- `indian`: The Indian Ocean
- `arctic`: The Arctic Ocean, which is by default rotated to the equator

The boundaries of these ocean basins are given by the International Hydrographic Organisation (IHO)'s *Limits of Oceans and Seas*.
Differing from the slepian_alpha package, the coastline data in this package is from GSHHG (A Global Self-consistent, Hierarchical, High-resolution Geography Database).

Along with the ocean basins, the following geographic domains are also supported:

- `earthquakes`: A mask of coastal megathrust earthquakes.

A new class is introduced as an interface to geographic domains:

- `GeoDomain`: A class that supports a few new methods that make fetching vertices and defining file names easier.

### Computing sea-level fingerprints

- `solvesle`: Solves the elastic sea-level equationfor a given forcing and ocean function
- `grace2fingerprint`: Uses GRACE data and GIA models to compute the sea-level fingerprints

### Projection of other fields to Slepian functions

The CSR GRACE mascon data (all components), various GIA models, and the sea-level fingerprints (all components) are natively supported to be projected to Slepian functions. The functions are:

- `gia2plmt`/`giaz2plmt`: Fetch the the geoid deformation/vertical displacement due to GIA projected to SH coefficients.
- `gia2slept`/`giaz2slept`: Fetch the geoid deformation/vertical displacement due to GIA projected to Slepian coefficients. This and the above function replace `correct4gia` from slepian_delta.
- `grace2trend`: Compute the GRACE mass trend in a given domain.
- `mascon2plmt`: Fetch the CSR GRACE mascon data projected to Slepian coefficients.
- `mascon2slept`: Fetch the CSR GRACE mascon data projected to Slepian coefficients.
- `mass2weq`: Converts mass data to water equivalent.
- `slep2xyz`: Converts Slepian functions to a mesh on a sphere.
- `slfa2plmt`: Fetch the sea-level fingerprints calculated by [Adhikari et al. (2019)](https://doi.org/10.5194/essd-11-629-2019).
- `slfa2slept`: Project the sea-level fingerprints calculated by [Adhikari et al. (2019)](https://doi.org/10.5194/essd-11-629-2019) to Slepian coefficients.

### Modifications to functions from other packages

The following functions have also been modified (mainly to support `GeoDomain`), and can probably safely replace the original functions:

- `correct4gia` (from slepian_delta) → `correct4gia_new`, although this function should be archived
- `glmalpha` (from slepian_alpha) → `glmalpha_new`
- `grace2plmt` (from slepian_delta)
- `grace2slept` (from slepian_delta) → `grace2slept_new`
- `integratebasis` (from slepian_delta) → `integratebasis_new` (not fully tested yet)
- `kernelcp` (from slepian_alpha) → `kernelcp_new`
- `plm2slep` (from slepian_bravo) → `plm2slep_new`
- `plm2xyz` (from slepian_alpha)
- `slep2plm` (from slepian_bravo) → `slep2plm_new`
- `slep2resid` (from slepian_bravo) → `slep2resid_new`

The ultimate goal is to drop the `_new` suffixes and replace the original functions with the modified ones. Note that these functions are only tested for geographic domains, and may not work as expected for, e.g. circular caps.

Additionally, `correct4gia` (from slepian_delta) has been split into two functions: `gia2plmt` and `gia2slept`.

### Visualisations

- `eigenwmesh`: Returns a mesh of the eigenvalue-weighted power map of the given Slepian functions.
- `equalearth` and `equalearthd`: Project a longitude-latitude vertex to the Equal Earth projection.
- `loadbasemap`: Loads an axesm-based map of the given geographic domain.
- `plotqdm`: Plots a quick-and-dirty map (i.e. on a normal axes object) of the given coordinates.

Other supporting functions for visualisation include:
`formatlonticks`, `fotmatlatticks`, `loadcbar`, `loadcmap`, and more.

---
Last modified by:

- [En-Chi Lee (@williameclee)](https://github.com/williameclee), 2025/05/19
