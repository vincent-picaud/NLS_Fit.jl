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
# TODO: use Tuple{LEFT_TYPE,RIGHT_TYPE) -> avoid if ... as we can use [1] and [2]
#
struct Model2Fit_Sum{LEFT_TYPE<:Abstract_Model2Fit,
                     RIGHT_TYPE<:Abstract_Model2Fit} <: Abstract_Model2Fit
    _left_right::Tuple{LEFT_TYPE,RIGHT_TYPE}
end

# Internal
#
get_submodel_θ_1(model::Model2Fit_Sum,θ::AbstractVector) = @view θ[1:parameter_size(model._left_right[1])]
get_submodel_θ_2(model::Model2Fit_Sum,θ::AbstractVector) = @view θ[(parameter_size(model._left_right[1])+1):end]

# Visit  ================
#
visit_submodel_size(model::Model2Fit_Sum) = 2

visit_get_submodel(model::Model2Fit_Sum,submodel_idx::Int) = model._left_right[submodel_idx]
visit_get_Y(model::Model2Fit_Sum,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = Y
visit_get_X(model::Model2Fit_Sum,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = X
function visit_get_θ(model::Model2Fit_Sum,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert submodel_idx==1 || submodel_idx==2

    if submodel_idx==1
        return get_submodel_θ_1(model,θ)
    end
    
    # submodel_idx == 2
    get_submodel_θ_2(model,θ)
end

# Interface  ================
#
parameter_size(m::Model2Fit_Sum) = parameter_size(m._left_right[1])+ parameter_size(m._left_right[2])


function accumulate_y!(m::Model2Fit_Sum,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
        @assert length(X) == length(Y)
    
    accumulate_y!(m._left_right[1],Y,X,get_submodel_θ_1(m,θ))
    accumulate_y!(m._left_right[2],Y,X,get_submodel_θ_2(m,θ))
    
    Y
end

# Some Base function overloadings  ================
#
import Base: (+)
Base. +(left::Abstract_Model2Fit,right::Abstract_Model2Fit) = Model2Fit_Sum((left,right))
Base. +(left::Model2Fit_Empty,right::Abstract_Model2Fit) = right
Base. +(left::Abstract_Model2Fit,right::Model2Fit_Empty) = left


