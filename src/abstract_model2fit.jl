export parameter_size, eval_y!, eval_y, alloc_y

@doc raw"""

Abstract type, base of all model to fit.
"""
abstract type Abstract_Model2Fit end

parameter_size(::Abstract_Model2Fit) = @assert(false,"To implement!")

eval_y!(::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)  = @assert(false,"To implement!")

# convenience

function alloc_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector)
    zeros(promote_type(eltype(X),eltype(θ)),length(X))
end

eval_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector) = eval_y!(m,alloc_y(m,X,θ),X,θ)
