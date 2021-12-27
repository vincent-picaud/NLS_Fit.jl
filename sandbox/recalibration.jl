using Revise
using NLS_Fit
using DelimitedFiles
using LinearAlgebra

rootDir = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
dataFile = joinpath(dataDir,"recalibration.txt")

# Getting started example ****************
#

# Data ================
#
XY=readdlm(dataFile);
X=XY[:,1];
Y=XY[:,2];

# Model without recalibration ================
#
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak() 


θ1 = Float64[1.0,5,1]
θ2 = Float64[1.5,10,2]
θ3 = Float64[0.75,20,3]
θ_uncalibrated_model = vcat(θ1,θ2,θ3)

# Add an affine calibration ================
#
recalibration_map = Map_Affine_Monotonic(X[1],X[end])
recalibration_model = Model2Fit_Recalibration(model,recalibration_map)

θ_map = Float64[1,1]
θ_init_recalibration_model = vcat(θ_uncalibrated_model, θ_map)

# Solve a bound constrained nonlinear least squares ================
#
nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y);
conf = NLS_Solver.LevenbergMarquardt_Conf()
result = NLS_Solver.solve(nls,θ_init,conf)

# Fit result ================
#
@assert NLS_Solver.converged(result)

θ_true = Float64[2,25,3] # cheating...
θ_solver = NLS_Solver.solution(result)

println("Found $θ_solver, error is $(norm(θ_solver-θ_true,Inf))")
println("(caveat: error is mainly due to noisy data, not due to solver accuracy)")

# Check solver accuracy by using a noise free input ****************
#

# Generate noise free input ================
#
n = 50
model = Gaussian_Peak()
θ = Float64[2,n/2,3]

X=Float64[1:n;];
Y=eval_y(model,X,θ);

# Solve and print new noise free results ================
#
nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y);
result = NLS_Solver.solve(nls,θ_init,conf)
θ_solver = NLS_Solver.solution(result)

println("Found $θ_solver, error is $(norm(θ_solver-θ_true,Inf))")
