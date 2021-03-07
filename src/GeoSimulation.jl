# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoSimulation

using Meshes
using GeoStatsBase
using Variography
using KrigingEstimators

using Distributions
using LinearAlgebra
using Statistics
using FFTW
using CpuId

import GeoStatsBase: preprocess, solve, solvesingle

include("lu.jl")
include("fft.jl")
include("seq.jl")
include("sgs.jl")
include("cookie.jl")

export
  # generic solvers
  SeqSim,

  # concrete solvers
  LUGS, FFTGS, SGS,

  # meta solvers
  CookieCutter

end
