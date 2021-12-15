export Map_Translate

# ****************************************************************

@doc raw"""
```julia
Map_Translate(;scale = 1)
```
Defines a translate transformation:

```math
X(\hat{X}) = \hat{X} + scale (θ - 1)
```

In peculiar, for ``(θ) = (1)`` the transformation is identity. 
This is due to the presence of the ``θ- 1`` factor. This factor is introduced to stay consistent with other transformations where we
try to have identity transform for parameter vector of ones.
"""
struct Map_Translate{T} <: Abstract_Map
    _scale::T
end

function Map_Translate(;scale=1)
    Map_Translate{typeof(scale)}(scale)
end

parameter_size(m::Map_Translate) = 1

function eval_map(m::Map_Translate{T},
                  X_hat::AbstractVector{T_HAT},
                  θ::AbstractVector{T_θ}) where {T_HAT,T,T_θ}

    @assert length(θ) == parameter_size(m)

    θ = θ[1]

    X_T = promote_type(T_HAT,T,T_θ)
    X=similar(X_hat,X_T)

    @. X = X_hat + m._scale*(θ-1)

    X
end
