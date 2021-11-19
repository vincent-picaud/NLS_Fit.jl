@doc raw"""

Map base type. 

This type is an abstraction of ``\hat{X}↦X``.

Such transformations are used for recalibration or parameters varying
with the space variable ``X``.

# Interface

- parameter_size
- eval_x

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
eval_x(m::Abstract_Map,X_hat::AbstractVector,θ::AbstractVector) -> X::AbstractVector
```

Compute ``X=X(\hat{X})``.
"""
eval_x(m::Abstract_Map,X_hat::AbstractVector,θ::AbstractVector) = @assert(false,"To implement")   
