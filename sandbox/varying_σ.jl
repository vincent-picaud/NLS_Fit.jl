using BenchmarkTools
using Revise
using NLS_Fit
using NLS_Solver
using DelimitedFiles
using StatsBase: sample


# ****************************************************************

XY=readdlm("../data/varying_σ.txt")
X=XY[:,1]
Y=XY[:,2]

model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]

θ_model = vcat(θ1,θ2,θ3)

# Create an evolution law for shape factor σ, depending on peak position

map_pos2sigma = Map_Affine{Float64}(1.0 => 1,30 => 2)
θ_map = ones(NLS_Fit.parameter_size(map_pos2sigma))
             
indices = [3,6,9]
pos = [5.0, 10.0, 20.0]

model_with_σ_law = Model2Fit_Mapped_Parameters(model,map_pos2sigma,indices,pos)
θ_model_with_σ_law = vcat(deleteat!(θ_model,indices),θ_map)

eval_y(model_with_σ_law,X,θ_model_with_σ_law)

# θ = vcat(θ1,θ2,θ3,θc)

# nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
# conf = Levenberg_Marquardt_Conf()
# result = solve(nls,θ,conf)

# eval_x(model,X,solution(result))
# eval_x(model,X,θ)




# @benchmark eval_x($model,$X,$θ)

