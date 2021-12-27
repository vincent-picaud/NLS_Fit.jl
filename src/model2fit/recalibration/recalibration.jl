# ****************************************************************
# Implement a generic Recalibration procedure using an Abstract_Map
# ****************************************************************
#
# TODO: change name to Model2Fit_Recalibration
#
export Model2Fit_Recalibration
export get_calibration_map,get_calibration_map_θ
export get_calibrated_model,get_calibrated_model_θ
export eval_calibrated_x

@doc raw"""
Create a calibrable model

# Constructors

```julia
calibrable_model = Model2Fit_Recalibration(model_to_calibrate, calibration_map)
```

where `model_to_calibrate` is an [`Abstract_Model2Fit`](@ref) and
`calibration_map` is an [`Abstract_Map`](@ref)

# Explanation

The call

```julia
eval_y(calibrable_model,X_hat,[θ, θ_map])
```

is equivalent to 

```julia
eval_y(model_to_calibrate,eval_calibrated_x(X_hat, θ_map), [θ])
```

# Also see

- [`get_calibration_map`](@ref) 
- [`get_calibration_map_θ`](@ref) 
- [`get_calibrated_model`](@ref) 
- [`get_calibrated_model_θ`](@ref) 
- [`eval_calibrated_x`](@ref) 

"""
struct Model2Fit_Recalibration{MODEL2FIT_TYPE <: Abstract_Model2Fit, MAP_TYPE <: Abstract_Map}  <: Abstract_Model2Fit
    _model2calibrate::MODEL2FIT_TYPE
    _map::MAP_TYPE
end

# Specific methods  ================
#

@doc raw"""
```julia
get_calibration_map(model::Model2Fit_Recalibration)::Abstract_Map
```

Retrieve the calibration map, of type [`Abstract_Map`](@ref).

See: [`Model2Fit_Recalibration`](@ref) 
"""
get_calibration_map(model::Model2Fit_Recalibration) = model._map

@doc raw"""
```julia
get_calibration_map_θ(model::Model2Fit_Recalibration,
                      θ::AbstractVector)::AbstractVector
```

Retrieve the calibration map parameters

See: [`Model2Fit_Recalibration`](@ref) 
"""
function get_calibration_map_θ(model::Model2Fit_Recalibration,θ::AbstractVector)
    @assert length(θ) == parameter_size(model)

    submodel = get_calibrated_model(model)
    s=parameter_size(submodel)

    @view θ[(s+1):end]
end


@doc raw"""
```julia
eval_calibrated_x(m::Model2Fit_Recalibration,X_hat::AbstractVector,θ::AbstractVector) -> AbstractVector
```

Compute the calibrated ``X`` from the reference domain ``\hat{X}`` for
the given transformation parameters ``θ``.

See: [`Model2Fit_Recalibration`](@ref) 
"""
function eval_calibrated_x(model::Model2Fit_Recalibration,X_hat::AbstractVector,θ::AbstractVector)
    map = get_calibration_map(model)
    map_θ = get_calibration_map_θ(model,θ)
    
    eval_map(map, X_hat, map_θ)
end

@doc raw"""
```julia
get_calibrated_model(model::Model2Fit_Recalibration)::Abstract_Model2Fit
```

Retrieve the underlying model to calibrate

See: [`Model2Fit_Recalibration`](@ref) 
"""
get_calibrated_model(model::Model2Fit_Recalibration) = model._model2calibrate

@doc raw"""
```julia
get_calibrated_model_θ(model::Model2Fit_Recalibration,
                       θ::AbstractVector)::AbstractVector
```
Retrieve the parameter vector of the underlying model to calibrate

See: [`Model2Fit_Recalibration`](@ref) 
"""
function get_calibrated_model_θ(model::Model2Fit_Recalibration,θ::AbstractVector)
    @assert length(θ) == parameter_size(model)

    submodel = get_calibrated_model(model)
    s=parameter_size(submodel)

    @view θ[1:s]
end

# Visit  ================
#
visit_submodel_size(model::Model2Fit_Recalibration) = 1

visit_get_submodel(model::Model2Fit_Recalibration,submodel_idx::Int) = get_calibrated_model(model)
visit_get_Y(model::Model2Fit_Recalibration,submodel_idx::Int,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector) = Y
visit_get_X(model::Model2Fit_Recalibration,submodel_idx::Int,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector) = eval_calibrated_x(model,X_hat,θ)
visit_get_θ(model::Model2Fit_Recalibration,submodel_idx::Int,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector) = get_calibrated_model_θ(model,θ)


# Interface ================
#
parameter_size(m::Model2Fit_Recalibration) = parameter_size(get_calibrated_model(m))+parameter_size(get_calibration_map(m))

function accumulate_y!(m::Model2Fit_Recalibration,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector)
  
    calibrated_X = eval_calibrated_x(m,X_hat,θ)
    calibrated_model = get_calibrated_model(m)
    calibrated_model_θ = get_calibrated_model_θ(m,θ)

    accumulate_y!(calibrated_model,Y,calibrated_X,calibrated_model_θ)
end

