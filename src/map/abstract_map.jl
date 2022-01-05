@doc raw"""
```julia
abstract type Abstract_Map end
```

This type is an abstraction of a 
```math
f:X ↦ Y=f_{\hat{\theta}}(X)
```
where ``\hat{\theta}`` is map parameter vector.

Such transformations can be used in various contexts: 

- calibration task: we have a reference ``\hat{X}`` to a calibrated
  ``X``, see [`Model2Fit_Recalibration`](@ref).

- parameter transformation: we map a parameter vector ``\hat{θ}`` to a new one ``θ``,
  see [`Model2Fit_Transformed_Parameters`](@ref).

# Interface

- [`parameter_size`](@ref)
- [`eval_map`](@ref) 

"""
abstract type Abstract_Map end

@doc raw"""
```julia
parameter_size(::Abstract_Map) -> Int
```

Return ``\hat{θ}`` parameter vector length

# See also
- [`Abstract_Map`](@ref) 

"""
parameter_size(::Abstract_Map) = @assert(false,"To implement")

@doc raw"""
```julia
eval_map(m::Abstract_Map,X_hat::AbstractVector,hat_θ::AbstractVector) -> X::AbstractVector
```

Compute ``X=X(\hat{X})``.

# See also
- [`Abstract_Map`](@ref) 

"""
eval_map(m::Abstract_Map,X_hat::AbstractVector,hat_θ::AbstractVector) = @assert(false,"To implement")   
