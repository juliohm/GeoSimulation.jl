@testset "LUGS" begin
  𝒮 = georef((z=[0.,1.,0.,1.,0.],), [0. 25. 50. 75. 100.])
  𝒟 = CartesianGrid(100)

  # ----------------------
  # conditional simulation
  # ----------------------
  problem = SimulationProblem(𝒮, 𝒟, :z, 2)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),))

  Random.seed!(2018)
  sol = solve(problem, solver)

  if visualtests
    @test_reference "data/LU-condsim.png" plot(sol,layout=(2,1))
  end

  # ------------------------
  # unconditional simulation
  # ------------------------
  problem = SimulationProblem(𝒟, :z=>Float64, 2)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),))

  Random.seed!(2018)
  sol = solve(problem, solver)

  if visualtests
    @test_reference "data/LU-uncondsim.png" plot(sol,layout=(2,1))
  end

  # -------------
  # co-simulation
  # -------------
  𝒟 = CartesianGrid(500)
  problem = SimulationProblem(𝒟, (:z=>Float64,:y=>Float64), 1)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),),
                :y => (variogram=GaussianVariogram(range=10.),),
                (:z,:y) => (correlation=0.95,))

  Random.seed!(2020)
  sol = solve(problem, solver)

  if visualtests
    @test_reference "data/LU-cosim.png" plot(sol,layout=(2,1))
  end

  # ---------------------
  # custom factorization
  # ---------------------
  𝒮 = georef((z=[0.,1.,0.,1.,0.],), [0. 25. 50. 75. 100.])
  𝒟 = CartesianGrid(100)
  problem = SimulationProblem(𝒮, 𝒟, :z, 1)
  solver1 = LUGS(:z => (variogram=SphericalVariogram(range=10.),factorization=lu))
  solver2 = LUGS(:z => (variogram=SphericalVariogram(range=10.),factorization=cholesky))

  Random.seed!(2021)
  solution1 = solve(problem, solver1)

  Random.seed!(2021)
  solution2 = solve(problem, solver2)

  if visualtests
    p1 = plot(solution1)
    p2 = plot(solution2)
    @test_reference "data/LU-factorization.png" plot(p1, p2, layout=(2,1))
  end
end
