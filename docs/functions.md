# Functions

## Overview

The public functions in ULMO are organised into several categories based on their purpose. Each function is documented with its name, a brief description, any package dependencies, and data dependencies.

Helper functions are not listed here but can be found in the source code (under `aux/`).

Many funcions have been modified from their original versions in the dependent packages (notably `slepian_alpha` and `slepian_delta`).
The ultimate goal is to drop the `_new` suffixes and replace the original functions with the modified ones. Note that these functions are only tested for geographic domains, and may not work as expected for, e.g. circular caps.

## Functions by Category

### Domain Processing

| Function name | Description | Package dependencies | Data dependencies |
|:---|:---|:---|:---|
| `GeoDomain` (class) | Stores and provides methods for geographic domain information. | [slepian_alpha]((https://github.com/csdms-contrib/slepian_alpha.git)), MATLAB Mapping Toolbox | none |
| `gshhscoastline` | Returns the global coastline data. | | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/) |
| `mass2weq` | Converts total mass change to equivalent sea level change. | | |
| `oceans` | Returns the boundary of the 'global' ocean domain *for GMOM* computation as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `alloceans` | Returns the boundary of the true global ocean domain *for solving the SLE* as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `npacific` | Returns the boundary of the North Pacific Ocean as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `spacific` | Returns the boundary of the South Pacific Ocean as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `natlantic` | Returns the boundary of the North Atlantic Ocean as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `satlantic` | Returns the boundary of the South Atlantic Ocean as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `indian` | Returns the boundary of the Indian Ocean as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `earthquakes` | Returns the boundary of the megathrust earthquake mask as used in Lee & Harig. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `arctic` | Returns the boundary of the Arctic Ocean. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `pacific` | Returns the boundary of the Pacific Ocean. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |
| `atlantic` | Returns the boundary of the Atlantic Ocean. | MATLAB Mapping Toolbox | [GSHHG](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/), [Limits of oceans and seas](https://doi.pangaea.de/10.1594/PANGAEA.777975) |

### Spherical Harmonic/Slepian Function Processing

| Function name | Description | Package dependencies | Data dependencies | Notes |
|:---|:---|:---|:---|:---|
| `integratebasis_new` | Integrate the Slepian functions within the domain. | slepian_delta | none | Replacement for `integratebasis` in slepian_alpha, with support for `GeoDomain` objects. |
| `kernelcp_new` | Computes the localisation matrix of Slepian functions in parallel. | slepian_alpha | none | Replacement for `kernelcp` in slepian_alpha, with support for `GeoDomain` objects. |
| `glmalpha_new` | Computes the projection matrix of Slepian functions. | slepian_alpha | none | Replacement for `glmalpha` in slepian_alpha, with support for `GeoDomain` objects. |
| `plm2slep_new` | Projects spherical harmonics onto Slepian basis. | slepian_alpha | none | Replacement for `plm2slep` in slepian_alpha, with support for `GeoDomain` objects. |
| `plm2xyz` | Projects spherical harmonics onto latitude-longitude grids/points. | slepian_alpha | none | Replacement for `plm2xyz` in slepian_alpha. |
| `slep2plm_new` | Projects Slepian functions back to spherical harmonics. | slepian_alpha | none | Replacement for `slep2plm` in slepian_alpha, with support for `GeoDomain` objects. |
| `slep2xyz` | Projects Slepian functions onto latitude-longitude grids/points. | slepian_alpha | none | |
| `slept2resid_new` | Fits a time series | slepian_alpha | none | Replacement for `slept2resid` in slepian_delta, with support for `GeoDomain` objects. |
| `sphgaussfilt` | Apply Gaussian filter to spherical harmonics | none | none | |

### GRACE Data Processing

| Function name | Description | Package dependencies | Data dependencies | Notes |
|:---|:---|:---|:---|:---|
| `aod1b2plmt` | Reads and converts GRACE GAC/GAD data to `lmcosi`-format spherical harmonic fields. | | GRACE | |
| `gracedeg1` | Reads and converts GRACE TN-13 degree-1 data to appropriate field and format. | | GRACE | Replacement for `gracedeg1` in slepian_delta. |
| `gracedeg2` | Reads and converts GRACE TN-14 $\Delta C_{2,0}$/$\Delta C_{3,0}$ data to appropriate field and format. | | GRACE | Replacement for `gracedeg2` in slepian_delta. |

### Sea Level Fingerprint/Geocentric Motion Computation

| Function name | Description | Package dependencies | Data dependencies | Notes |
|:---|:---|:---|:---|:---|
| `lovenumber` | Reads the elastic Love numbers for geoid/VLM/HLM. | | [ISSM](https://github.com/ISSMteam/ISSM.git) | Replacement for `lovenums` in slepian_delta. |
| `solvesle` | Solves the elastic/transient SLE for a given mass load and land water mask. | | none | |
| `solvedegree1` | Solves the Stokes coefficients for the geocentric motion from GRACE data. | | GRACE | |
| `grace2fingerprint` | Solves the SLF from GRACE data on land. | | GRACE | |
| `slfa2plmt` | Reads and converts the *RSL* SLF calculated by [Adhikari et al. (2019)](https://doi.org/10.5194/essd-11-629-2019) to spherical harmonic fields. | | [Adhikari et al. (2019)](https://doi.org/10.5194/essd-11-629-2019) | |
| `slfa2slept` | Projects the *RSL* SLF calculated by [Adhikari et al. (2019)](https://doi.org/10.5194/essd-11-629-2019) to Slepian coefficients. | slepian_alpha | See `slfa2plmt` | |

### GIA Model Processing

| Function name | Description | Package dependencies | Data dependencies | Notes |
|:---|:---|:---|:---|:---|
| `gia2plmt` | Reads and converts *geoid* data from different GIA models to spherical harmonic fields. | | | Replacement and extension for `correct4gia` in [slepian_delta](https://github.com/csdms-contrib/slepian_delta.git). |
| `gia2slept` | Projects *geoid* data from different GIA models to Slepian functions of a specific domain. | | See `gia2plmt` | Replacement and extension for `correct4gia` in [slepian_delta](https://github.com/csdms-contrib/slepian_delta.git). |
| `giaz2plmt` | Reads and converts *VLM* data from different GIA models to spherical harmonic fields. | | | |
| `gia2slept` | Projects *VLM* data from different GIA models to Slepian functions of a specific domain. | | See `giaz2plmt` | |

### Altimetric/Steric Data Processing

| Function name | Description | Package dependencies | Data dependencies |
|:---|:---|:---|:---|
| `ssh2lonlatt` | Reads and converts altimetric SSH data from different missions to gridded data. | | |
| `steric2lonlatt` | Reads and pre-processed steric sea level data from different datasets to gridded data. | | |

### Visualisation

| Function name | Description | Package dependencies | Data dependencies |
|:---|:---|:---|:---|
| `eigenwmesh` | Returns a mesh of the eigenvalue-weighted power map of the given Slepian functions. | | |
| `loadbasemap` | Loads an axesm-based map of the given geographic domain. | | |
| `plotqdm` | Plots a quick-and-dirty map (i.e. on a normal axes object) of the given coordinates. | | |

### Miscellaneous

| Function name | Description | Package dependencies | Data dependencies |
|:---|:---|:---|:---|
| `periodictimeseries` | Fits polynomials and sinusoidal functions to a time series. | | |

## Appendices

### Acronyms and Abbreviations

| Abbreviation | Definition |
|:---|:---|
| GIA | glacial isostatic adjustment |
| GMOM | global mean ocean mass |
| SLE | sea level equation |
| SLF | sea level fingerprint |
| SSH | sea surface height |
| TN | GRACE Technical Note |

---

Last modified: 2025/10/16, [@williameclee](mailto:williameclee@arizona.edu)
