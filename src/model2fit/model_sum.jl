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
    _left::LEFT_TYPE
    _right::RIGHT_TYPE
end

# Visit  ================
#
visit_submodel_size(model::Model2Fit_Sum) = 2

function visit_get_submodel(model::Model2Fit_Sum,submodel_idx::Int)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)

    if submodel_idx==1
        return model._left
    end

    # submodel_idx == 2
    model._right
end

function visit_get_X(model::Model2Fit_Sum,submodel_idx::Int,X::AbstractVector,θ::AbstractVector)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)
    @assert length(θ) == parameter_size(model)
    
    X
end

function visit_get_θ(model::Model2Fit_Sum,submodel_idx::Int,X::AbstractVector,θ::AbstractVector)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)
    @assert length(θ) == parameter_size(model)
    
    if submodel_idx==1
        return @view θ[1:parameter_size(model._left)]
    end

    # submodel_idx == 2
    @view θ[(parameter_size(model._left)+1):end]
end

# Interface  ================
#
parameter_size(m::Model2Fit_Sum) =
    parameter_size(visit_get_submodel(m,1)) +
    parameter_size(visit_get_submodel(m,2))

function accumulate_y!(m::Model2Fit_Sum,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(X) == length(Y)

    accumulate_y!(visit_get_submodel(m,1),Y,X,visit_get_θ(m,1,X,θ))
    accumulate_y!(visit_get_submodel(m,2),Y,X,visit_get_θ(m,2,X,θ))

    Y
end

# Some Base function overloadings  ================
#
import Base: (+)
Base. +(left::Abstract_Model2Fit,right::Abstract_Model2Fit) = Model2Fit_Sum(left,right)
Base. +(left::Model2Fit_Empty,right::Abstract_Model2Fit) = right
Base. +(left::Abstract_Model2Fit,right::Model2Fit_Empty) = left


