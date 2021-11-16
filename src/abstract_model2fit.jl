@doc raw"""

Abstract type, base of all model to fit.
"""
abstract type Abstract_Model2Fit end

parameter_size(::Abstract_Model2Fit) = @assert(false,"To implement!")

# Here eval Y_i for point X_i. When performing fit we have a collection of such X_i
#
eval_y(::Abstract_Model2Fit,X_i::Any,θ::AbstractVector)  = @assert(false,"To implement!")
