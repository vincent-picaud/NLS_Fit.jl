export Model2Fit_Transformed_Parameters
export get_model_θ, get_model

@doc raw"""
Create a new model, where ``\theta`` parameters are computed using an
[`Abstract_Map`](@ref).

```math
\hat{m}(X,\hat{\theta)) = m(X,\theta = f(\hat{\theta))
```

# Construtor

```julia
Model2Fit_Transformed_Parameters(model:: Abstract_Model2Fit,
                                 map::Abstract_Map)
```


"""
struct Model2Fit_Transformed_Parameters{MODEL <: Abstract_Model2Fit,
                                        MAP <: Abstract_Map} <: Abstract_Model2Fit
    _model::MODEL
    _map::MAP

    function Model2Fit_Transformed_Parameters(model::MODEL,
                                              map::MAP) where {MODEL <: Abstract_Model2Fit,
                                                               MAP <: Abstract_Map}
        
        new{MODEL,MAP}(model, map)
    end
    
end
# Specific methods  ================
#
@doc raw"""
```julia
get_model(mp::Model2Fit_Transformed_Parameters)::Abstract_Model2Fit
```

Get back the wrapped model
"""
get_model(model::Model2Fit_Transformed_Parameters) = model._model

@doc raw"""
```julia
get_model_θ(mp::Model2Fit_Transformed_Parameters,θ::AbstractVector) -> θ::AbstractVector
```

Retrieve the parameter vector θ associated to the wrapped model [`get_model`](@ref).
"""
function get_model_θ(hat_model::Model2Fit_Transformed_Parameters,hat_θ::AbstractVector)
    θ = eval_map(hat_model._map,hat_θ)

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
parameter_size(hat_model::Model2Fit_Transformed_Parameters) = parameter_size(hat_model._map)

function accumulate_y!(hat_model::Model2Fit_Transformed_Parameters,Y::AbstractVector,X::AbstractVector,hat_θ::AbstractVector)
    model = get_model(hat_model)
    θ = get_model_θ(hat_model,hat_θ)

    accumulate_y!(model,Y,X,θ)
end
