export Model2Fit_Transformed_Parameters
export get_model_θ, get_model
export create_model_transform_src_insert_dest

@doc raw"""

Create a new model, where the  parameter vector ``\theta`` is computed using an
[`Abstract_Map`](@ref):
```math
 \hat{\theta}\mapsto \theta = f_{\hat{\theta}_f}(\hat{\theta}_m)
```

The resulting model is thus:
```math
\hat{m}(X,\hat{\theta}=[\hat{\theta}_m,\hat{\theta}_f]) = m(X,\theta = f_{\hat{\theta}_f}(\hat{\theta}_m))
```

# Constructor

## General case

```julia
Model2Fit_Transformed_Parameters(model::Abstract_Model2Fit,
                                 map_domain_size::Int,
                                 map::Abstract_Map)
```

- `map_domain_size` is the expected ``\hat{\theta}_m`` size. This
  quantity is used to implement model `parameter_size()` method. The
  reason is that the current [`Abstract_Map`](@ref) interface only
  provides ``\hat{\theta}_f`` size (its `parameter_size()` method) but
  tells nothing about ``\hat{\theta}_m`` size.
- Also, please note that ``\theta`` size is not necessary equal to
  ``\hat{\theta}_m`` size

## Using `Transform_Src_Insert_Dest_Map`

Very often this structure is used with a map of type
[`Transform_Src_Insert_Dest_Map`](@ref). In that case it is
simpler to use a dedicated function:
[`create_model_transform_src_insert_dest`](@ref) (follow this link to
see an example).

# Example

TODO: add an example using the general constructor

# Extra methods

Beside the regular [`Abstract_Model2Fit`](@ref) interface, you can use these extra methods:

- [`get_model(model::Model2Fit_Transformed_Parameters)`](@ref) 
- [`get_model_θ(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)`](@ref) 

## Not exported

These methods are not exported but may be useful:
- [`get_model_hat_θ_view(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)`](@ref) 
- [`get_map_hat_θ_view(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)`](@ref) 
"""
struct Model2Fit_Transformed_Parameters{MODEL <: Abstract_Model2Fit,
                                        MAP <: Abstract_Map} <: Abstract_Model2Fit
    _model::MODEL
    _map_domain_size::Int # length(X) where map(X,θ̂f), mandatory to define parameter_size
    # _map_codomain_size must be equal to the wrapped _model parameter_size
    _map::MAP

    function Model2Fit_Transformed_Parameters(model::MODEL,
                                              map_domain_size::Int,
                                              map::MAP) where {MODEL <: Abstract_Model2Fit,
                                                               MAP <: Abstract_Map}
        new{MODEL,MAP}(model, map_domain_size, map)
    end

end

@doc raw"""

This function creates an [`Model2Fit_Transformed_Parameters`](@ref)
instance from a map of type
[`Transform_Src_Insert_Dest_Map`](@ref).

```julia
function create_model_transform_src_insert_dest(model, map, src=>dest) -> model
```

# Example
   
```jldoctest
using NLS_Fit

# Two Gaussian peaks with σ(μ) an affine function
#
model = Gaussian_Peak() + Gaussian_Peak() 

# model θ is:
# idx:  1   2  3   4   5   6
#  θ : h1, μ1, σ1, h2, μ2, σ2
#
# hence:
#
src  = [2,5]
dest = [3,6]

# The σ(μ) relation: σ(μ=1)=1 and σ(μ=100)=10 for [θ̂A, θ̂B]=[1, 1]
#
f_σ_μ = Map_Affine(1.0=>1.0,100.0=>10.0)

# Now the model_σ_μ parameters are:
#
# [ h1, μ1, h2, μ2, θ̂A, θ̂B ]
#
# the initial model is called with:
#
# [ h1, μ1, σ1 = affine(μ1, θ̂A, θ̂B), h2, μ2, , σ2 = affine(μ2, θ̂A, θ̂B) ]
#
model_σ_μ = create_model_transform_src_insert_dest(model,f_σ_μ,src=>dest)

# Examine that the compute θ takes into account the σ(μ) dependance:
#
θ̂ = Float64[1,5, 2, 95, 1, 1]
θ = get_model_θ(model_σ_μ, θ̂)

# output
6-element Vector{Float64}:
  1.0
  5.0
  1.3636363636363638
  2.0
 95.0
  9.545454545454547

```

"""
function create_model_transform_src_insert_dest(model::Abstract_Model2Fit,
                                                f_map::Abstract_Map,
                                                src_dest::Pair{SOURCE,DEST}) where {SOURCE<:AbstractVector{Int},DEST<:AbstractVector{Int}}
    @assert length(first(src_dest)) == length(last(src_dest))
    
    # Initialize g_map (see Transform_Src_Insert_Dest_Map
    #
    g_map = NLS_Fit.Transform_Src_Insert_Dest_Map(f_map, src_dest)
    
    # compute map_domain_size:
    # -> this is model θ length minus the number of inserted element
    #
    map_domain_size = parameter_size(model) - length(first(src_dest))
    @assert map_domain_size ≥ 0
        
    Model2Fit_Transformed_Parameters(model, map_domain_size, g_map)
end 

# Specific methods  ================
#
@doc raw"""
```julia
get_model_hat_θ_view(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)
```

Return a view on ``\hat{\theta}_m`` knowing ``\hat{\theta}=[\hat{\theta}_m,\hat{\theta}_f]``
"""
function get_model_hat_θ_view(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)
    @assert parameter_size(hat_model) == length(hat_θ)
    
    @view hat_θ[1:hat_model._map_domain_size]
end

@doc raw"""
```julia
get_model_hat_θ_view(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)
```

Return a view on ``\hat{\theta}_f`` knowing ``\hat{\theta}=[\hat{\theta}_m,\hat{\theta}_f]``
"""
function get_map_hat_θ_view(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)
    @assert parameter_size(hat_model) == length(hat_θ)

    @view hat_θ[(hat_model._map_domain_size+1):end]
end

@doc raw"""
```julia
get_model(mp::Model2Fit_Transformed_Parameters)::Abstract_Model2Fit
```

Get back the wrapped model

See [`Model2Fit_Transformed_Parameters`](@ref) 
"""
get_model(model::Model2Fit_Transformed_Parameters) = model._model

@doc raw"""
```julia
get_model_θ(mp::Model2Fit_Transformed_Parameters, hat_θ::AbstractVector) -> θ::AbstractVector
```

Retrieve the parameter vector θ associated to the wrapped model [`get_model`](@ref).
"""
function get_model_θ(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)
    model_hat_θ = get_model_hat_θ_view(hat_model,hat_θ)
    map_hat_θ = get_map_hat_θ_view(hat_model,hat_θ)

    θ = eval_map(hat_model._map,model_hat_θ,map_hat_θ)

    # check that map codomain = model parameter size
    @assert parameter_size(get_model(hat_model)) == length(θ)

    θ
end 

# Visit  ================
#
visit_submodel_size(hat_model::Model2Fit_Transformed_Parameters) = 1
visit_get_submodel(hat_model::Model2Fit_Transformed_Parameters,submodel_idx::Int) = get_model(hat_model)
visit_get_X(hat_model::Model2Fit_Transformed_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,hat_θ::AbstractVector) = X
visit_get_Y(hat_model::Model2Fit_Transformed_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,hat_θ::AbstractVector) = Y
visit_get_θ(hat_model::Model2Fit_Transformed_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,hat_θ::AbstractVector) = get_model_θ(hat_model,hat_θ)


# Interface ================
#
parameter_size(hat_model::Model2Fit_Transformed_Parameters) = hat_model._map_domain_size + parameter_size(hat_model._map)

function accumulate_y!(hat_model::Model2Fit_Transformed_Parameters,Y::AbstractVector,X::AbstractVector,hat_θ::AbstractVector)
    model = get_model(hat_model)
    θ = get_model_θ(hat_model,hat_θ)

    accumulate_y!(model,Y,X,θ)
end
