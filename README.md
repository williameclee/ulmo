# Unified Localisation of Mass over Oceans (ULMO)

A pipeline for using GRACE satellite gravimetry, altimetry-steric measurements, and forward modelling to study ocean mass change over ocean basins, with extended functionalities for Slepian functions and sea level equation modelling.

Research findings using this code were presented at the [AGU24 Fall Meeting](https://agu.confex.com/agu/agu24/meetingapp.cgi/Paper/1579620). A manuscript is in preparation.

![Example of Slepian functions over the Pacific Ocean](images/cover.svg)

## Functionalities

See the [Functions](docs/functions.md) page for a list of public functions in this repository, along with their descriptions and dependencies.

## Installation

1. Download or clone this repository to your local machine.
2. Add the repository folder and subfolders to your MATLAB path.
3. Install the dependencies listed below and ensure their paths are configured properly, and the search order is set correctly (ULMO should be at the top of the order in `pathdef.m`).

## Dependencies

Functions in this repository may call or overwrite functions from the following packages. Please ensure they are installed and the paths are configured properly before running the functions in this repository.

- [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha.git): Computation of Slepian functions on the sphere
- [slepian_bravo](https://github.com/csdms-contrib/slepian_bravo.git): Conversion between spherical harmonics and Slepian functions
- [slepian_delta](https://github.com/csdms-contrib/slepian_delta.git): Processing GRACE data and GIA models
- [Gibbs Seawater Toolbox (GSW)](https://www.teos-10.org/software.htm): Computation of steric sea level from temperature and salinity

---

Last modified: 2025/11/03, [En-Chi Lee (@williameclee)](mailto:williameclee@arizona.edu)
