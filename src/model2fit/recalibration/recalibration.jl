export Recalibration
export eval_calibrated_x

# ****************************************************************
# Implement affine recalibration
# ****************************************************************
#
@doc raw"""
Create a calibrable model

# Constructors

```julia
calibrable_model = Recalibration_Affine(model_to_calibrate, calibration_map)
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

"""
struct Recalibration{MODEL2FIT_TYPE <: Abstract_Model2Fit, MAP_TYPE <: Abstract_Map}  <: Abstract_Model2Fit
    _model2calibrate::MODEL2FIT_TYPE
    _map::MAP_TYPE
end

# Specific methods  ================
#
@doc raw"""
```julia
eval_calibrated_x(m::Recalibration,X_hat::AbstractVector,θ::AbstractVector) -> AbstractVector
```

Compute the calibrated ``X`` from the reference domain ``\hat{X}`` for
the given transformation parameters ``θ``.

"""
function eval_calibrated_x(model::Recalibration,X_hat::AbstractVector,θ::AbstractVector)
    submodel = model._model2calibrate
    s=parameter_size(submodel)
    θ_map = @view θ[(s+1):end]
    
    eval_x(model._map, X_hat, θ_map)
end

get_submodel(model::Recalibration) = model._model2calibrate
get_map(model::Recalibration) = model._map

function get_submodel_θ(model::Recalibration,θ::AbstractVector)
    submodel = get_submodel(model)
    s=parameter_size(submodel)
    θ_model = @view θ[1:s]
end

# Visit  ================
#
visit_submodel_size(model::Recalibration) = 1

visit_get_submodel(model::Recalibration,submodel_idx::Int) = get_submodel(model)
visit_get_Y(model::Recalibration,submodel_idx::Int,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector) = Y
visit_get_X(model::Recalibration,submodel_idx::Int,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector) = eval_calibrated_x(model,X_hat,θ)
visit_get_θ(model::Recalibration,submodel_idx::Int,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector) = get_submodel_θ(model,θ)


# Interface ================
#
parameter_size(m::Recalibration) = parameter_size(get_submodel(m))+parameter_size(get_map(m))

function accumulate_y!(m::Recalibration,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector)
  
    X_calibrated = eval_calibrated_x(m,X_hat,θ)
    submodel = get_submodel(m)
    submodel_θ = get_submodel_θ(m,θ)

    accumulate_y!(submodel,Y,X_calibrated,submodel_θ)
end

