@doc raw"""

Map base type. 

This type is an abstraction of a ``f:X ↦ Y=f(X)``.

Such transformations can be used in various contexts:
- calibration task: we have a reference ``\hat{X}`` to a calibrated ``X``
- parameter transformation: we map a parameter vector ``\theta\prime` to a new one ``\theta``.

# Interface

- [`parameter_size`](@ref)
- [`eval_map`](@ref) 

"""
abstract type Abstract_Map end

@doc raw"""
```julia
parameter_size(::Abstract_Map) -> Int
```

Return θ parameter vector length
"""
parameter_size(::Abstract_Map) = @assert(false,"To implement")

@doc raw"""
```julia
eval_map(m::Abstract_Map,X_hat::AbstractVector,θ::AbstractVector) -> X::AbstractVector
```

Compute ``X=X(\hat{X})``.
"""
eval_map(m::Abstract_Map,X_hat::AbstractVector,θ::AbstractVector) = @assert(false,"To implement")   
