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


# # Model
# #


# # σ(μ)
# map_μ_to_σ = Map_Affine(1.0=>1.0,30.0=>5.0)

# # map (μ,σ(μ))
# struct Varying_σ_Map <: NLS_Fit.Abstract_Map
#     _n_peak::Int
#     _map_μ_to_σ::NLS_Fit.Abstract_Map
# end

# NLS_Fit.parameter_size(varying_σ_map::Varying_σ_Map) =
#     2*varying_σ_map._n_peak + NLS_Fit.parameter_size(varying_σ_map._map_μ_to_σ)

# function NLS_Fit.eval_map(varying_σ_map::Varying_σ_Map,hat_θ::AbstractVector)
#     @assert parameter_size(varying_σ_map) == length(hat_θ)

#     n_peak = varying_σ_map._n_peak
#     map_μ_to_σ = varying_σ_map._map_μ_to_σ

#     input_θ =  @view hat_θ[1:(2*varying_σ_map._n_peak)]
#     map_μ_to_σ_θ = @view hat_θ[(2*varying_σ_map._n_peak + 1):end]

#     θ = Vector{eltype(hat_θ)}(undef,3*varying_σ_map._n_peak)

#     # create input views
#     view_h = @view input_θ[1:2:2*n_peak]
#     view_μ = @view input_θ[2:2:2*n_peak]
    
#     # copy height, μ and σ(μ)
#     θ[1:3:3*n_peak] .= view_h
#     θ[2:3:3*n_peak] .= view_μ
#     θ[3:3:3*n_peak] .= eval_map(map_μ_to_σ,view_μ,map_μ_to_σ_θ)

#     θ
# end


# t = Varying_σ_Map(4,map_μ_to_σ)
# parameter_size(t)
# eval_map(t,ones(parameter_size(t)))

# model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
# model = Model2Fit_Transformed_Parameters(model,Varying_σ_Map(3,map_μ_to_σ))

# parameter_size(model)
# θ_init = ones(8)
# θ_init[2]=4
# θ_init[4]=8
# θ_init[6]=25

# nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
# conf = NLS_Solver.LevenbergMarquardt_Conf()
# result = NLS_Solver.solve(nls,θ_init,conf)

# @assert NLS_Solver.converged(result)
# θ_fit = NLS_Solver.solution(result)
# get_model(model)
# get_model_θ(model,θ_fit)

# Y_init = eval_y(model,X,θ_init)
# Y_fit = eval_y(model,X,θ_fit)

# writedlm("test.txt",hcat(X,Y,Y_init,Y_fit))

#
source=collect(1:3:12)
dest=collect(2:3:12)

function src_indices_after_dest_removed(source_dest::Pair{SOURCE,DEST}) where {SOURCE<:AbstractVector{Int},DEST<:AbstractVector{Int}}
    (src,dest) = source_dest

    @assert length(src)==length(dest)

    n = length(src)

    count =  zeros(Int,n)
    i_count = 1

    # note this works even if src, dest are not sorted but this has
    # n^2 complexity
    for src_i in src
        for dest_i in dest
            if dest_i < src_i
                count[i_count] += 1
            end
        end 
        i_count += 1
    end
    src - count
end
src_indices_after_dest_removed(dest=>dest)

src_indices_after_dest_removed(source=>dest)

# An helper that  maps
#
# [ θ̂_1, θ̂_2, ... , θ̂_n, θ̂f_1, θ̂f_2, ..., θ̂f_nf ]
#
# to
#
# [ θ_1, θ_2, ... , θ_{n+m} ] (m is dim of src = dim of dest)
#
# where θ̂_i for i ∈ src have be taken, and f_{θf}(θ̂_i) inserted at pos dest
#
struct Transform_Src_Insert_Dest_Map <: NLS_Fit.Abstract_Map
    _map::NLS_Fit.Abstract_Map
    _src_dest::Matrix{Int} # src=col(1), dest=col(2)
end

function Transform_Src_Insert_Dest_Map(map::NLS_Fit.Abstract_Map, source_dest::Pair{SOURCE,DEST}) where {SOURCE<:AbstractVector{Int},DEST<:AbstractVector{Int}}
    @assert length(first(source_dest))==length(last(source_dest))
    
    src = src_indices_after_dest_removed(source_dest)
    dest = last(source_dest)

    Transform_Src_Insert_Dest_Map(map,hcat(src,dest))
end

NLS_Fit.parameter_size(tm::Transform_Src_Insert_Dest_Map) = NLS_Fit.parameter_size(tm._map)

function NLS_Fit.eval_map(tm::Transform_Src_Insert_Dest_Map,hat_model_θ::AbstractVector,hat_map_θ::AbstractVector)
    @assert parameter_size(tm) == length(hat_map_θ)

    src  = @view tm._src_dest[:,1]
    dest = @view tm._src_dest[:,2]
    
    θ_mapped = eval_map(tm._map,(@view hat_model_θ[src]),hat_map_θ)
   
    θ = insert_some_elements(hat_model_θ, dest, θ_mapped)

    θ
end

# Demo ****************
#
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
n = NLS_Fit.parameter_size(model)

map_μ_to_σ = Map_Affine(1.0=>1.0,30.0=>5.0)
source=collect(2:3:n)
dest=collect(3:3:n)

map = Transform_Src_Insert_Dest_Map(map_μ_to_σ,source=>dest)
model = Model2Fit_Transformed_Parameters(model,n-length(dest),map)

NLS_Fit.parameter_size(model)
NLS_Fit.get_model_θ(model,Float64[1:8;])

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

# parameter_size(t)
# hat_model_θ = Float64[1:12;]
# hat_model_θ = deleteat!(hat_model_θ,dest)
# hat_map_θ = ones(2)
# NLS_Fit.eval_map(t,hat_model_θ,hat_map_θ)

