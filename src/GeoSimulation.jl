# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoSimulation

using GeoStatsBase
using Variography
using KrigingEstimators
using Distributions
using LinearAlgebra
using Statistics
using FFTW
using CpuId

import GeoStatsBase: preprocess, solvesingle

include("lu.jl")
include("fft.jl")
include("seq.jl")

export LUGS, FFTGS, SGS

end
