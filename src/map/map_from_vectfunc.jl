export Map_From_VectFunc

@doc raw"""

An [`Abstract_Map`](@ref) instance that can be easily created by
providing a Julia function

# Constructor
```julia
Map_From_VectFunc(n_θ,f)
```
where `(X,θ)->f(X,θ)` is the wrapped function.

# Example

```jldoctest
using NLS_Fit

f(X,θ) = θ[1]*X

map = Map_From_VectFunc(1,f)

X = Float64[1,2,3]

θ = Float64[4]

eval_map(map,X,θ)

# output
3-element Vector{Float64}:
  4.0
  8.0
 12.0

```

"""
struct Map_From_VectFunc{F<:Function} <: Abstract_Map
    _n_θ::Int
    _f::F
end

parameter_size(m::Map_From_VectFunc) = m._n_θ

function eval_map(m::Map_From_VectFunc,
                  X::AbstractVector,
                  θ::AbstractVector)
    @assert length(θ) == parameter_size(m)

    m._f(X,θ)
end 
