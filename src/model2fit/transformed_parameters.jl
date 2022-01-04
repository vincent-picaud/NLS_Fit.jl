export Model2Fit_Transformed_Parameters
export get_model_θ, get_model

@doc raw"""

Create a new model, where ``\theta`` parameters are computed using an
[`Abstract_Map`](@ref) to define ``\hat{\theta}_m\mapsto
f_{\hat{\theta}_f}(\hat{\theta}_m)``.

```math
\hat{m}(X,\hat{\theta}=[\hat{\theta}_m,\hat{\theta}_f]) = m(X,\theta = f_{\hat{\theta}_f}(\hat{\theta}_m))
```

# Constructor

```julia
Model2Fit_Transformed_Parameters(model::Abstract_Model2Fit,
                                 map_domain_size::Int,
                                 map::Abstract_Map)
```

- `map_domain_size` is the expected ``\hat{\theta}_m`` size. This
  quantity is used to implement model `parameter_size()` method. The
  reason is that current [`Abstract_Map`](@ref) interface only
  provides ``\hat{\theta}_f`` size (its `parameter_size()` method) but
  tells nothing about ``\hat{\theta}_m`` size.
- Also, please note that ``\theta`` size is not necessary equal to
  ``\hat{\theta}_m`` size

# See
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
    _map_domain_size::Int # length(X) where map(X,θ̂_map), mandatory to define parameter_size
    # _map_codomain_size must be = to parameter_size(_model)
    _map::MAP

    function Model2Fit_Transformed_Parameters(model::MODEL,
                                              map_domain_size::Int,
                                              map::MAP) where {MODEL <: Abstract_Model2Fit,
                                                               MAP <: Abstract_Map}
        new{MODEL,MAP}(model, map_domain_size, map)
    end
    
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
