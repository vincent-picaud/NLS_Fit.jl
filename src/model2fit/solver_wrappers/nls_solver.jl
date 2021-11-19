export NLS_ForwardDiff_From_Model2Fit

import NLS_Solver
#: AbstractNLS, parameter_size, residue_size, eval_r, eval_r_J

using ForwardDiff

# Important: using "template" parameters allows to have 1 alloc in
# eval_r (versus 6 if one uses X,Y::AbstractVector)
#
#
# Finally MODEL2FIT_TYPE <: Abstract_Model2Fit is parametrized 
# julia> @btime eval_r($nls,$θ)
#   2.724 μs (44 allocations: 928 bytes)
# 10-element Vector{Float64}:
#  0.8056529547193196
#  0.10344414895180448
#  0.7746398216971467
#  0.5705008291976804
#  0.758504377425209
#  0.6572756090278895
#  0.22021418677443244
#  0.9012363674393303
#  0.25785439137881183
#  0.7594027591650859
#
# julia> @btime eval_r($nls,$θ)
#   842.590 ns (1 allocation: 160 bytes)
# 10-element Vector{Float64}:
#  -0.6484862692515055
#   0.10436137221132202
#   0.3416337767340343
#  -0.5195066030780292
#   0.0635334492635572
#  -0.2671004288239941
#  -0.13064413518496942
#  -0.143200160754688
#   0.517014266738562
#   0.21068675449678864
# 
struct NLS_ForwardDiff_From_Model2Fit{MODEL2FIT_TYPE <: Abstract_Model2Fit,
                                      X_ELEMENT_TYPE,
                                      Y_ELEMENT_TYPE,
                                      X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                      Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
    _fit_model::MODEL2FIT_TYPE
    _X::X_TYPE
    _Y::Y_TYPE
 
    function NLS_ForwardDiff_From_Model2Fit(fit_model::MODEL2FIT_TYPE,
                                            X::X_TYPE,Y::Y_TYPE) where{MODEL2FIT_TYPE <: Abstract_Model2Fit,
                                                                       X_ELEMENT_TYPE,
                                                                       Y_ELEMENT_TYPE,
                                                                       X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                                                       Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} 
        @assert length(X) == length(Y)

        new{MODEL2FIT_TYPE,X_ELEMENT_TYPE,Y_ELEMENT_TYPE,X_TYPE,Y_TYPE}(fit_model,X,Y)
    end
end



NLS_Solver.parameter_size(nls::NLS_ForwardDiff_From_Model2Fit) = parameter_size(nls._fit_model)
NLS_Solver.residue_size(nls::NLS_ForwardDiff_From_Model2Fit)  = length(nls._Y)

function NLS_Solver.eval_r(nls::NLS_ForwardDiff_From_Model2Fit,θ::AbstractVector)
    Y_model = eval_y(nls._fit_model,nls._X,θ)
    @. Y_model = nls._Y - Y_model
    Y_model
end


function NLS_Solver.eval_r_J(nls::NLS_ForwardDiff_From_Model2Fit, θ::AbstractVector)
    
    # r_evaluation = (r,θ)->(r .= NLS_Solver.eval_r(nls,θ))
    
    # r = Vector{T}(undef,NLS_Solver.residue_size(nls))

    #    J = ForwardDiff.jacobian(r_evaluation, r, θ)

    r = NLS_Solver.eval_r(nls,θ)
    J = ForwardDiff.jacobian(θ->NLS_Solver.eval_r(nls,θ), θ)
    r,J
end

