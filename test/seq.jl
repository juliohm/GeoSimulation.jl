@testset "SGS" begin
  𝒮 = georef((z=[1.,0.,1.],), [25. 50. 75.; 25. 75. 50.])
  𝒟 = RegularGrid(100,100)
  N = 3

  𝒫₁ = SimulationProblem(𝒮, 𝒟, :z, N)
  𝒫₂ = SimulationProblem(𝒟, :z=>Float64, N)

  solver = SGS(
    :z => (variogram=GaussianVariogram(range=35.),
           neighborhood=BallNeighborhood(10.))
  )

  Random.seed!(2017)
  sol₁ = solve(𝒫₁, solver)
  sol₂ = solve(𝒫₂, solver)

  # basic checks
  reals = sol₁[:z]
  inds = LinearIndices(size(𝒟))
  @test all(reals[i][inds[26,26]] == 1. for i in 1:N)
  @test all(reals[i][inds[51,76]] == 0. for i in 1:N)
  @test all(reals[i][inds[76,51]] == 1. for i in 1:N)
end
