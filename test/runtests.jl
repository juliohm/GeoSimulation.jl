using GeoSimulation
using Meshes
using GeoStatsBase
using Variography
using Distributions
using LinearAlgebra
using Plots; gr(size=(600,400))
using ReferenceTests, ImageIO
using Test, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

# dummy definitions
include("dummy.jl")

# list of tests
testfiles = [
  "lu.jl",
  "fft.jl",
  "seq.jl",
  "sgs.jl",
  "cookie.jl"
]

@testset "GeoSimulation.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
