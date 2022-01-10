export Map_Quadratic
export eval_map

# ****************************************************************

@doc raw"""
```julia
Map_Quadratic(X_hat_A => X_A, X_hat_B => X_B, X_hat_C => X_C)

Map_Quadratic(X_hat_A, X_hat_B, X_hat_C)
```

Same semantic as [`Map_Affine`](@ref) 
"""
struct Map_Quadratic{T} <: Abstract_Map
    _hat_to_A::Pair{T,T}
    _hat_to_B::Pair{T,T}
    _hat_to_C::Pair{T,T}
end

function Map_Quadratic(x_hat_A::T,x_hat_B::T,x_hat_C::T) where {T}
    Map_Quadratic(x_hat_A=>one(T),
                  x_hat_B=>one(T),
                  x_hat_C=>one(T))
end

parameter_size(m::Map_Quadratic) = 3

function eval_map(m::Map_Quadratic{T},
                  X_hat::AbstractVector{T_HAT},
                  θ::AbstractVector{T_θ}) where {T_HAT,T,T_θ}

    @assert length(θ) == parameter_size(m)
    
    X_T = promote_type(T_HAT,T,T_θ)

    θ_A = θ[1]
    θ_B = θ[2]
    θ_C = θ[3]

    X_hat_A, X_A = m._hat_to_A
    X_hat_B, X_B = m._hat_to_B
    X_hat_C, X_C = m._hat_to_C


    basis_A(x) = (X_hat_B - x)*(X_hat_C - x)/((X_hat_B - X_hat_A)*(X_hat_C - X_hat_A))
    basis_B(x) = (x - X_hat_A)*(X_hat_C - x)/((X_hat_B - X_hat_A)*(X_hat_C - X_hat_B))
    basis_C(x) = (x - X_hat_A)*(x - X_hat_B)/((X_hat_C - X_hat_A)*(X_hat_C - X_hat_B))

    transform(x) = basis_A(x)*X_A*θ_A + basis_B(x)*X_B*θ_B + basis_C(x)*X_C*θ_C

    X=similar(X_hat,X_T)
    @. X = transform(X_hat)

    X
end

