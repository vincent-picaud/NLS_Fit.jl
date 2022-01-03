using Revise
using NLS_Fit
using DelimitedFiles
using LinearAlgebra

rootDir = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
dataFile = joinpath(dataDir,"mapped_parameters.txt")

# Getting started example ****************
#

# Data ================
#
XY=readdlm(dataFile);
X=XY[:,1];
Y=XY[:,2];


# Model
#


# σ(μ)
map_μ_to_σ = Map_Affine(1.0=>1.0,30.0=>5.0)

# map (μ,σ(μ))
struct Varying_σ_Map <: NLS_Fit.Abstract_Map
    _n_peak::Int
    _map_μ_to_σ::NLS_Fit.Abstract_Map
end

NLS_Fit.parameter_size(varying_σ_map::Varying_σ_Map) =
    2*varying_σ_map._n_peak + NLS_Fit.parameter_size(varying_σ_map._map_μ_to_σ)

function NLS_Fit.eval_map(varying_σ_map::Varying_σ_Map,hat_θ::AbstractVector)
    @assert parameter_size(varying_σ_map) == length(hat_θ)

    n_peak = varying_σ_map._n_peak
    map_μ_to_σ = varying_σ_map._map_μ_to_σ

    input_θ =  @view hat_θ[1:(2*varying_σ_map._n_peak)]
    map_μ_to_σ_θ = @view hat_θ[(2*varying_σ_map._n_peak + 1):end]

    θ = Vector{eltype(hat_θ)}(undef,3*varying_σ_map._n_peak)

    # create input views
    view_h = @view input_θ[1:2:2*n_peak]
    view_μ = @view input_θ[2:2:2*n_peak]
    
    # copy height, μ and σ(μ)
    θ[1:3:3*n_peak] .= view_h
    θ[2:3:3*n_peak] .= view_μ
    θ[3:3:3*n_peak] .= eval_map(map_μ_to_σ,view_μ,map_μ_to_σ_θ)

    θ
end


t = Varying_σ_Map(4,map_μ_to_σ)
parameter_size(t)
eval_map(t,ones(parameter_size(t)))

model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
model = Model2Fit_Transformed_Parameters(model,Varying_σ_Map(3,map_μ_to_σ))

parameter_size(model)
θ_init = ones(8)
θ_init[2]=4
θ_init[4]=8
θ_init[6]=25

nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = NLS_Solver.LevenbergMarquardt_Conf()
result = NLS_Solver.solve(nls,θ_init,conf)

@assert NLS_Solver.converged(result)
θ_fit = NLS_Solver.solution(result)
get_model(model)
get_model_θ(model,θ_fit)

Y_init = eval_y(model,X,θ_init)
Y_fit = eval_y(model,X,θ_fit)

writedlm("test.txt",hcat(X,Y,Y_init,Y_fit))

