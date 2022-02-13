# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoSimulation

using Meshes
using GeoStatsBase
using Variography
using KrigingEstimators

using SpecialFunctions
using Distributions
using LinearAlgebra
using Statistics
using Random
using Tables
using FFTW
using CpuId

import GeoStatsBase: preprocess, solve, solvesingle

include("lu.jl")
include("fft.jl")
include("seq.jl")
include("sgs.jl")
include("spde.jl")
include("cookie.jl")

export
  # generic solvers
  SeqSim,

  # concrete solvers
  LUGS,
  FFTGS,
  SGS,
  SPDEGS,

  # meta solvers
  CookieCutter

end
