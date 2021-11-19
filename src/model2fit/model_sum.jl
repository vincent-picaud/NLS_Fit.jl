# Note:
# ----------------------------------------------------------------
#
# with:
# model = Gaussian_Peak() + Gaussian_Peak()
# θ=rand(parameter_size(model))
#
# julia> @btime eval_y($model,2.0,$θ)
#   281.779 ns (8 allocations: 208 bytes)
# 0.06722452348922268
# 
# struct Model2Fit_Sum <: Abstract_Model2Fit
#     _left::Abstract_Model2Fit
#     _right::Abstract_Model2Fit
# end
#
# versus:
#
# struct Model2Fit_Sum{LEFT_TYPE<:Abstract_Model2Fit,
#                      RIGHT_TYPE<:Abstract_Model2Fit} <: Abstract_Model2Fit
#     _left::LEFT_TYPE
#     _right::RIGHT_TYPE
# end
#
# julia> @btime eval_y($model,2.0,$θ)
#   77.596 ns (0 allocations: 0 bytes)
# 0.018749729970880733
#
# -> one must use parametrized types to limit memory allocations
#
struct Model2Fit_Sum{LEFT_TYPE<:Abstract_Model2Fit,
                     RIGHT_TYPE<:Abstract_Model2Fit} <: Abstract_Model2Fit
    _left::LEFT_TYPE
    _right::RIGHT_TYPE
end

parameter_size(m::Model2Fit_Sum) = parameter_size(m._left) + parameter_size(m._right)

function eval_y!(m::Model2Fit_Sum,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(X) == length(Y)
    
    eval_y!(m._left,Y,X,@view θ[1:parameter_size(m._left)]) 
    eval_y!(m._right,Y,X,@view θ[(parameter_size(m._left)+1):end])

    Y
end


import Base: (+)
Base. +(left::Abstract_Model2Fit,right::Abstract_Model2Fit) = Model2Fit_Sum(left,right)
Base. +(left::Model2Fit_Empty,right::Abstract_Model2Fit) = right
Base. +(left::Abstract_Model2Fit,right::Model2Fit_Empty) = left


