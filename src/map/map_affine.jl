export Map_Affine
export eval_map

# ****************************************************************

@doc raw"""
Defines an affine transformation:

```math
X(\hat{X}) = L_A(\hat{X}) X_A θ_A +  L_B(\hat{X}) X_B θ_B
```
where ``L_A`` and ``L_B`` are the Lagrange basis.

# Constructors

## Scaled interpolation 

```julia
Map_Affine(X_hat_A => X_A, X_hat_B => X_B)
```

for ``(θ_A,θ_B) = (1,1)`` we have an affine map such that
``X(\hat{X}_A) = X_A`` and ``X(\hat{X}_B) = X_B``.

## Usual interpolation

The second constructor 
```julia
Map_Affine(X_hat_A, X_hat_B)
```

is a shortcut for

```julia
Map_Affine(X_hat_A => 1, X_hat_B => 1)
```
that defines an **usual** interpolation (without scaling):
``(θ_A,θ_B)`` are simply the values at ``(X_A,X_B)``
"""
struct Map_Affine{T} <: Abstract_Map
    _hat_to_A::Pair{T,T}
    _hat_to_B::Pair{T,T}
end

function Map_Affine(x_hat_A::T,x_hat_B::T) where {T}
    Map_Affine(x_hat_A=>one(T),x_hat_B=>one(T))
end

parameter_size(m::Map_Affine) = 2

function eval_map(m::Map_Affine{T},
                X_hat::AbstractVector{T_HAT},
                θ::AbstractVector{T_θ}) where {T_HAT,T,T_θ}

    @assert length(θ) == parameter_size(m)
    
    X_T = promote_type(T_HAT,T,T_θ)

    θ_A = θ[1]
    θ_B = θ[2]

    X_hat_A = first(m._hat_to_A)
    X_A = last(m._hat_to_A)
    X_hat_B = first(m._hat_to_B)
    X_B = last(m._hat_to_B)

    Δ = 1/(X_hat_B-X_hat_A)

    basis_A(x) = (X_hat_B-x)*Δ
    basis_B(x) = (x-X_hat_A)*Δ
    transform(x) = basis_A(x)*X_A*θ_A + basis_B(x)*X_B*θ_B

    X=similar(X_hat,X_T)
    @. X = transform(X_hat)

    X
end



# note: maybe better to use:
# @article{article,
# author = {Murray, Kevin and Müller, Samuel and Turlach, Berwin},
# year = {2016},
# month = {02},
# pages = {1-21},
# title = {Fast and flexible methods for monotone polynomial fitting},
# volume = {86},
# journal = {Journal of Statistical Computation and Simulation},
# doi = {10.1080/00949655.2016.1139582}
# }                               
