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
recalibration_map = Map_Affine(X[1]=>X[1],X[end]=>X[end])
recalibration_model = Model2Fit_Recalibration(model,recalibration_map)

θ_map = Float64[1,1]
θ_init_recalibration_model = vcat(θ_uncalibrated_model, θ_map)

# Solve a bound constrained nonlinear least squares ================
#

# Objective function ----------------
#
nls = NLS_ForwardDiff_From_Model2Fit(recalibration_model,X,Y);

# Bound constraints ----------------
#
ε = eps(Float64)
#                       h1 μ1   σ1,  h2  μ2   σ2   h3   μ3   σ3  θA   θB
lower_bound = Float64[   0, 5,   ε,   0, 10,   ε,   0,  20,   ε, 0.5, 0.5]
upper_bound = Float64[ Inf, 5, Inf, Inf, 10, Inf, Inf, 20, Inf, 1.5, 1.5]

bc = NLS_Solver.BoundConstraints(lower_bound,upper_bound)

# Solver ----------------
#
conf = NLS_Solver.LevenbergMarquardt_Conf()

# Solve the problem ----------------
#
result = NLS_Solver.solve(nls,θ_init_recalibration_model,conf)

# Fit result ================
#
@assert NLS_Solver.converged(result)

θ_fit_recalibration_model = NLS_Solver.solution(result)
Y_fit_recalibration_model = eval_y(recalibration_model,X,θ_fit_recalibration_model)

# Retrieve calibrated X ----------------
#
X_recalibrated = eval_calibrated_x(recalibration_model,X,θ_fit_recalibration_model);
X_true = (X .- 0.2)/1.1;

norm(X_recalibrated-X,Inf)
norm(X_recalibrated-X_true,Inf)


get_calibrated_model(recalibration_model)
get_calibrated_model_θ(recalibration_model,θ_fit_recalibration_model)
