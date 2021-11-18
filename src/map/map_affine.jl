@doc raw"""
```julia
Map_Affine(X_hat_A => X_A, X_hat_B => X_B)

Map_Affine(X_hat_A, X_hat_B)
```

The second constructor is a shortcut for
```julia
Map_Affine(X_hat_A => X_hat_A, X_hat_B => X_hat_B)
```

Define an affine transformation:

```math
X(\hat{X}) = L_A(\hat{X}) X_A θ_A +  L_B(\hat{X}) X_B θ_B
```

where `L_A` and `L_B` are the Lagrange basis.

In peculiar, for `(θ_A,θ_B) = (1,1)` we have an affine map such that
`X(X_hat_A) = X_A` and `X(X_hat_B) = X_B`.
"""
struct Map_Affine{T} <: Abstract_Map
    _hat_to_A::Pair{T,T}
    _hat_to_B::Pair{T,T}
end

function Map_Affine(x_hat_A::T,x_hat_B::T) where {T}
    Map_Affine(x_hat_A=>x_hat_A,x_hat_B=>x_hat_B)
end

parameter_size(m::Map_Affine) = 2

function eval_x(m::Map_Affine{T},
                X_hat::AbstractVector{T_HAT},
                θ::AbstractVector{T_θ}) where {T_HAT,T,T_θ}
    
    X_T = promote_type(T_HAT,T,T_θ)

    θ_A = θ[1]
    θ_B = θ[2]
    
    Δ = 1/(first(m._hat_to_B)-first(m._hat_to_A))

    basis_A(x) = (first(m._hat_to_B)-x)*Δ*last(m._hat_to_A)
    basis_B(x) = (x-first(m._hat_to_A))*Δ*last(m._hat_to_B)
    transform(x) = basis_A(x)*θ_A + basis_B(x)*θ_B

    X=similar(X_hat,X_T)
    @. X = transform(X_hat)

    X
end
