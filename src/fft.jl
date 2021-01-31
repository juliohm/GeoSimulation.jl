# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    FFTGS(var₁=>param₁, var₂=>param₂, ...)

FFT Gaussian simulation.

## Parameters

* `variogram` - theoretical variogram (default to `GaussianVariogram()`)
* `mean`      - mean of Gaussian field (default to `0`)

## Global parameters

* `threads` - number of threads in FFT (default to all physical cores)

### References

Gutjahr 1997. *General joint conditional simulations using a fast
Fourier transform method.*
"""
@simsolver FFTGS begin
  @param variogram = GaussianVariogram()
  @param mean = 0.0
  @global threads = cpucores()
end

function preprocess(problem::SimulationProblem, solver::FFTGS)
  hasdata(problem) && @warn "Conditional spectral Gaussian simulation is not currently supported"
  
  # retrieve problem info
  pdomain = domain(problem)
  npts = nelms(pdomain)
  dims = size(pdomain)
  center = CartesianIndex(dims .÷ 2)
  c = LinearIndices(dims)[center]

  # number of threads in FFTW
  FFTW.set_num_threads(solver.threads)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  # result of preprocessing
  preproc = Dict()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine value type
      V = mactypeof[var]

      # determine variogram model and mean
      γ = varparams.variogram
      μ = varparams.mean

      # check stationarity
      @assert isstationary(γ) "variogram model must be stationary"

      # compute covariances between centroid and all locations
      covs = sill(γ) .- pairwise(γ, pdomain, [c], 1:npts)
      C = reshape(covs, dims)

      # move to frequency domain
      F = sqrt.(abs.(fft(fftshift(C))))
      F[1] = zero(V) # set reference level

      # save preprocessed inputs for variable
      preproc[var] = (γ=γ, μ=μ, F=F)
    end
  end

  preproc
end

function solvesingle(problem::SimulationProblem, covars::NamedTuple, ::FFTGS, preproc)
  # retrieve problem info
  pdomain = domain(problem)
  dims = size(pdomain)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  varreal = map(covars.names) do var
    # unpack preprocessed parameters
    γ, μ, F = preproc[var]

    # determine value type
    V = mactypeof[var]

    # perturbation in frequency domain
    P = F .* exp.(im .* angle.(fft(rand(V, dims))))

    # move back to time domain
    Z = real(ifft(P))

    # adjust mean and variance
    σ² = Statistics.var(Z, mean=zero(V))
    Z .= √(sill(γ) / σ²) .* Z .+ μ

    # flatten result
    var => vec(Z)
  end

  Dict(varreal)
end
