export Abstract_Model2Fit
export parameter_size, accumulate_y!, eval_y, alloc_y

@doc raw"""

Abstract type, base of all model to fit.

# Interface

The interface is very simple

- [`parameter_size`](@ref) 
- [`accumulate_y!`](@ref) 

# Convenience functions

- [`eval_y`](@ref) 

"""
abstract type Abstract_Model2Fit end

# Visit  ================
#
visit_submodel_size(model::Abstract_Model2Fit) = @assert(false,"To implement for $(typeof(model)) type!")
visit_get_submodel(model::Abstract_Model2Fit,submodel_idx::Int) = @assert(false,"To implement for $(typeof(model)) type!")
visit_get_Y(model::Abstract_Model2Fit,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = @assert(false,"To implement for $(typeof(model)) type!")
visit_get_X(model::Abstract_Model2Fit,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = @assert(false,"To implement for $(typeof(model)) type!")
visit_get_θ(model::Abstract_Model2Fit,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = @assert(false,"To implement for $(typeof(model)) type!")

# Public functions ================
#
@doc raw"""
```julia
parameter_size(::Abstract_Model2Fit)::Int
```

Return length of the expected parameter vector `θ`

# Also see 
- [`Abstract_Model2Fit`](@ref)
"""
parameter_size(model::Abstract_Model2Fit) = @assert(false,"To implement for $(typeof(model)) type!")

@doc raw"""
```julia
accumulate_y!(::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)::AbstractVector
```

Accumulate model contribution into vector `Y`.

# Also see 
- [`Abstract_Model2Fit`](@ref)
- [`eval_y`](@ref) 
"""
accumulate_y!(model::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)  = @assert(false,"To implement for $(typeof(model)) type!")

# convenience

function alloc_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector)
    zeros(promote_type(eltype(X),eltype(θ)),length(X))
end

@doc raw"""
```julia
eval_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector)::AbstractVector
```

A convenience function that call [`accumulate_y!`](@ref) using a zero
initialized `Y` vector. This returned vector contains model values.

# Also see 
- [`Abstract_Model2Fit`](@ref)
"""
eval_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector) = accumulate_y!(m,alloc_y(m,X,θ),X,θ)
