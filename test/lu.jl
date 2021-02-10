@testset "LUGS" begin
  𝒮 = georef((z=[0.,1.,0.,1.,0.],), [0. 25. 50. 75. 100.])
  𝒟 = RegularGrid(100)

  # ----------------------
  # conditional simulation
  # ----------------------
  problem = SimulationProblem(𝒮, 𝒟, :z, 2)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),))

  Random.seed!(2018)
  solution = solve(problem, solver)

  if visualtests
    @test_ref_plot "data/LU-condsim.png" plot(solution,size=(600,400),layout=(2,1))
  end

  # ------------------------
  # unconditional simulation
  # ------------------------
  problem = SimulationProblem(𝒟, :z=>Float64, 2)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),))

  Random.seed!(2018)
  solution = solve(problem, solver)

  if visualtests
    @test_ref_plot "data/LU-uncondsim.png" plot(solution,size=(600,400),layout=(2,1))
  end

  # -------------
  # co-simulation
  # -------------
  𝒟 = RegularGrid(500)
  problem = SimulationProblem(𝒟, (:z=>Float64,:y=>Float64), 1)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),),
                :y => (variogram=GaussianVariogram(range=10.),),
                (:z,:y) => (correlation=0.95,))

  Random.seed!(2020)
  solution = solve(problem, solver)

  if visualtests
    @test_ref_plot "data/LU-cosim.png" plot(solution,size=(600,400),layout=(2,1))
  end
end
