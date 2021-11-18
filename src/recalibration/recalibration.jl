export Recalibration

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

where `model_to_calibrate` is an [`AbstractModel`](@ref) and
`calibration_map` is an [`Abstract_Map`](@ref)

# Explanation

The call

```julia
eval_y(calibrable_model,X_hat,[θ, θ_map])
```

is equivalent to 

```julia
eval_y(model_to_calibrate,eval_x(X_hat, θ_map), [θ])
```

"""
struct Recalibration{MODEL2FIT_TYPE <: Abstract_Model2Fit, MAP_TYPE <: Abstract_Map}  <: Abstract_Model2Fit
    _model2calibrate::MODEL2FIT_TYPE
    _map::MAP_TYPE
end

parameter_size(m::Recalibration) = parameter_size(m._model2calibrate)+parameter_size(m._map)

@doc raw"""
```julia
eval_x(m::Recalibration,X_hat::AbstractVector,θ::AbstractVector) -> AbstractVector
```

Compute the calibrated ``X`` from the reference domain ``\hat{X}`` for
the given transformation parameters ``θ``.

"""
function eval_x(m::Recalibration,X_hat::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)

    s=parameter_size(m._model2calibrate)
    θ_map = @view θ[(s+1):end]
 
    eval_x(m._map, X_hat, θ_map)

end

function eval_y!(m::Recalibration,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector)
  
    X = eval_x(m._map, X_hat, θ)
    
    s=parameter_size(m._model2calibrate)
    θ_model = @view θ[1:s]

    eval_y!(m._model2calibrate,Y,X,θ_model)
end

