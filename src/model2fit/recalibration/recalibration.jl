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

# Visit  ================
#
visit_submodel_size(model::Recalibration) = 1

_visit_get_submodel(model::Recalibration) = model._model2calibrate

function visit_get_submodel(model::Recalibration,submodel_idx::Int)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)

    _visit_get_submodel(model)
end

function _visit_get_X(model::Recalibration,X_hat::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(model)

    submodel = _visit_get_submodel(model)
    s=parameter_size(submodel)
    θ_map = @view θ[(s+1):end]
    
    eval_x(model._map, X_hat, θ_map)
end 
function visit_get_X(model::Recalibration,submodel_idx::Int,X_hat::AbstractVector,θ::AbstractVector)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)

    _visit_get_X(model,X_hat,θ)
end

function _visit_get_θ(model::Recalibration,θ::AbstractVector)
    submodel = _visit_get_submodel(model)
    s=parameter_size(submodel)
    θ_model = @view θ[1:s]
end
function visit_get_θ(model::Recalibration,submodel_idx::Int,X::AbstractVector,θ::AbstractVector)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)

    _visit_get_θ(model,θ)
end

# Interface ================
#
parameter_size(m::Recalibration) = parameter_size(m._model2calibrate)+parameter_size(m._map)

@doc raw"""
```julia
eval_calibrated_x(m::Recalibration,X_hat::AbstractVector,θ::AbstractVector) -> AbstractVector
```

Compute the calibrated ``X`` from the reference domain ``\hat{X}`` for
the given transformation parameters ``θ``.

"""
function eval_calibrated_x(m::Recalibration,X_hat::AbstractVector,θ::AbstractVector)
    _visit_get_X(m,X_hat,θ)
end

function accumulate_y!(m::Recalibration,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector)
  
    X_calibrated = _visit_get_X(m,X_hat,θ)
    θ_model = _visit_get_θ(m,θ)

    accumulate_y!(m._model2calibrate,Y,X_calibrated,θ_model)
end

