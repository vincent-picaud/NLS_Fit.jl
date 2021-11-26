export parameter_size, eval_y!, eval_y, alloc_y

@doc raw"""

Abstract type, base of all model to fit.

# Interface

- [`parameter_size`](@ref) 
- [`eval_y!`](@ref) 

"""
abstract type Abstract_Model2Fit end

@doc raw"""
```julia
parameter_size(::Abstract_Model2Fit)::Int
```

Return length of the expected parameter vector `θ`


Also see : 
- [`Abstract_Model2Fit`](@ref)
"""
parameter_size(::Abstract_Model2Fit) = @assert(false,"To implement!")

@doc raw"""
```julia
eval_y!(::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)::AbstractVector
```

Accumulate model contribution into vector `Y`.

Also see : 
- [`Abstract_Model2Fit`](@ref)
- [`eval_y`](@ref) 
"""
eval_y!(::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)  = @assert(false,"To implement!")

# convenience

function alloc_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector)
    zeros(promote_type(eltype(X),eltype(θ)),length(X))
end

@doc raw"""
```julia
eval_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector)::AbstractVector
```

A convenience function that call [`eval_y!`](@ref) using a zero
initialized `Y` vector. This returned vector contains model values.

Also see : 
- [`Abstract_Model2Fit`](@ref)
"""
eval_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector) = eval_y!(m,alloc_y(m,X,θ),X,θ)
