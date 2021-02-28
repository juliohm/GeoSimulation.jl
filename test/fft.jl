@testset "FFTGS" begin
  𝒟 = CartesianGrid(100,100)
  problem = SimulationProblem(𝒟, :z=>Float64, 3)

  Random.seed!(2019)
  solver = FFTGS(:z => (variogram=GaussianVariogram(range=10.),))
  solution = solve(problem, solver)

  if visualtests
    @test_reference "data/FFT-iso.png" plot(solution,size=(900,300))
  end

  Random.seed!(2019)
  solver = FFTGS(:z => (variogram=GaussianVariogram(distance=aniso2distance([20.,5.],[0.])),))
  solution = solve(problem, solver)

  if visualtests
    @test_reference "data/FFT-aniso.png" plot(solution,size=(900,300))
  end
end
