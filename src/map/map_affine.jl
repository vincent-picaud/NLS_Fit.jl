export Map_Affine, Map_Affine_Monotonic
export eval_x

# ****************************************************************

@doc raw"""
```julia
Map_Affine(X_hat_A => X_A, X_hat_B => X_B)

Map_Affine(X_hat_A, X_hat_B)
```

The second constructor is a shortcut for
```julia
Map_Affine(X_hat_A => X_hat_A, X_hat_B => X_hat_B)
```

Defines an affine transformation:

```math
X(\hat{X}) = L_A(\hat{X}) X_A θ_A +  L_B(\hat{X}) X_B θ_B
```

where ``L_A`` and ``L_B`` are the Lagrange basis.

In peculiar, for ``(θ_A,θ_B) = (1,1)`` we have an affine map such that
``X(\hat{X}_A) = X_A`` and ``X(\hat{X}_B) = X_B``.
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

# ****************************************************************

# The only difference is that θ_B is remplaced by θ_A + θ_B

@doc raw"""
```julia
Map_Affine_Monotonic(X_hat_A => X_A, X_hat_B => X_B)

Map_Affine_Monotonic(X_hat_A, X_hat_B)
```

The second constructor is a shortcut for
```julia
Map_Affine_Monotonic(X_hat_A => X_hat_A, X_hat_B => X_hat_B)
```
Same as [`Map_Affine`](@ref) but uses another parametrization
that allows to insure monotonic map using simple bound constraints.

The `Map_Affine` map is defined as follows:
```math
X(\hat{X}) = L_A(\hat{X}) X_A θ_A +  L_B(\hat{X}) X_B θ_B
```

With this parametrization you have to add the ``θ_B ≥ θ_A`` to insure
that the map is increasing.

On the other hand, with the `Map_Affine_Monotonic`
parametrization:

```math
X(\hat{X}) = L_A(\hat{X}) X_A θ_A +  L_B(\hat{X}) ( (X_B-X_A) θ_B + X_A θ_A )
```

with 
```math
\frac{d}{d\hat{X}} X = \frac{X_B-X_A}{\hat{X}_B-\hat{X}_A} θ_B
```

we see that the ``θ_B`` parameter now plays the role a slope factor
and you simply have to impose ``θ_B ≥ 0`` to preserve the increasing
or decreasing character given by the sign of ``
\frac{X_B-X_A}{\hat{X}_B-\hat{X}_A}`` (which is constant once the
structure has been initialized).

In peculiar, for ``(θ_A,θ_B) = (1,1)`` we have an affine map such that
``X(\hat{X}_A) = X_A`` and ``X(\hat{X}_B) = X_B``. Whereas for
``(θ_A,θ_B) = (1,0)`` we have a constant one.

"""
struct Map_Affine_Monotonic{T} <: Abstract_Map
    _hat_to_A::Pair{T,T}
    _hat_to_B::Pair{T,T}
end

function Map_Affine_Monotonic(x_hat_A::T,x_hat_B::T) where {T}
    Map_Affine_Monotonic(x_hat_A=>x_hat_A,x_hat_B=>x_hat_B)
end

parameter_size(m::Map_Affine_Monotonic) = 2

function eval_x(m::Map_Affine_Monotonic{T},
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
    transform(x) = basis_A(x)*X_A*θ_A + basis_B(x)*((X_B-X_A)*θ_B+X_A*θ_A)

    X=similar(X_hat,X_T)
    @. X = transform(X_hat)

    X
end
