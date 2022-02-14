# GeoSimulation.jl

[![][build-img]][build-url] [![][codecov-img]][codecov-url]

Geostatistical simulation solvers for the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

### LUGS

This solver provides an implementation of direct Gaussian simulation (a.k.a. LU simulation)
as described in [Alabert 1987](https://link.springer.com/article/10.1007/BF00897191). In this
method, the full covariance matrix is built to include all locations of the simulation domain,
and samples from the multivariate Gaussian are drawn via LU factorization.

The method, which is widely implemented in many packages for Gaussian processes (e.g.
[GaussianProcesses.jl](https://github.com/STOR-i/GaussianProcesses.jl)),
is appropriate for relatively small simulation domains (e.g. 100x100 grids) where it is feasible
to factorize the full covariance. For larger domains (e.g. 3D grids), other methods are available
such as sequential Gaussian simulation, spectral methods, and FFT moving averages.

### FFTGS

This solver provides an implementation of spectral Gaussian simulation (a.k.a. FFT simulation)
as described in [Gutjahr 1997](https://link.springer.com/article/10.1007/BF02769641).
In this method, the covariance function is perturbed in the frequency
domain after a fast Fourier transform. White noise is added to the phase
of the spectrum, and a realization is produced by an inverse Fourier transform.

The method is limited to simulations on regular grids, and care must be taken
to make sure that the correlation length is small enough compared to the grid
size. As a general rule of thumb, avoid correlation lengths greater than 1/3
of the grid. The method is extremely fast, and can be used to generate large
3D realizations.

### SGS

This solver provides an implementation of sequential Gaussian simulation as described in
[Gomez-Hernandez & Journel 1993](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8).
The method traverses all locations of the spatial domain according to a path, approximates the
conditional distribution function within a neighborhood using Kriging, and assigns a value to
the center of the neighborhood by sampling from this distribution.

### SPDEGS

This solver provides an implementation of the stochastic partial differential equations
(SPDE) approach to Gaussian simulation [Lindgren et al. 2011.](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x).
The method relies on a discretization of the Laplace-Beltrami operator on meshes and is
adequate for highly curved domains (e.g. surfaces).

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add GeoSimulation
```

## Usage

This package is part of the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

For a simple example of usage, please check the main documentation.

## Asking for help

If you have any questions, please contact our community on the [gitter channel](https://gitter.im/JuliaEarth/GeoStats.jl).

[build-img]: https://img.shields.io/github/workflow/status/JuliaEarth/GeoSimulation.jl/CI?style=flat-square
[build-url]: https://github.com/JuliaEarth/GeoSimulation.jl/actions

[codecov-img]: https://img.shields.io/codecov/c/github/JuliaEarth/GeoSimulation.jl?style=flat-square
[codecov-url]: https://codecov.io/gh/JuliaEarth/GeoSimulation.jl

