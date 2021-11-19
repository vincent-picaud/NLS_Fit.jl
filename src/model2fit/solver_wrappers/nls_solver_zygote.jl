# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# An attempt of wrapping using Zygote
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Does not work at all in our case:
# julia> @btime eval_r_J($nls_zygote,$θ)
#   1.421 s (2632107 allocations: 81.61 MiB) 
#
# versus Forward diff
#
# julia> @btime eval_r_J($nls,$θ)
#   36.179 μs (9 allocations: 26.11 KiB)
#
# using Zygote

# struct NLS_Zygote_From_Fit_Model{MODEL2FIT_TYPE <: Abstract_Model2Fit,
#                                  X_ELEMENT_TYPE,
#                                  Y_ELEMENT_TYPE,
#                                  X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
#                                  Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
#     _fit_model::MODEL2FIT_TYPE
#     _X::X_TYPE
#     _Y::Y_TYPE
 
#     function NLS_Zygote_From_Fit_Model(fit_model::MODEL2FIT_TYPE,
#                                             X::X_TYPE,Y::Y_TYPE) where{MODEL2FIT_TYPE <: Abstract_Model2Fit,
#                                                                        X_ELEMENT_TYPE,
#                                                                        Y_ELEMENT_TYPE,
#                                                                        X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
#                                                                        Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} 
#         @assert length(X) == length(Y)

#         new{MODEL2FIT_TYPE,X_ELEMENT_TYPE,Y_ELEMENT_TYPE,X_TYPE,Y_TYPE}(fit_model,X,Y)
#     end
# end


# NLS_Solver.parameter_size(nls::NLS_Zygote_From_Fit_Model) = parameter_size(nls._fit_model)
# NLS_Solver.residue_size(nls::NLS_Zygote_From_Fit_Model)  = length(nls._Y)

# function NLS_Solver.eval_r(nls::NLS_Zygote_From_Fit_Model,θ::AbstractVector{T}) where T
#     # Cannot use adjoint as zip() currently not supported 
#     #    map(((X_i,Y_i);)->Y_i-eval_y(nls._fit_model,X_i,θ),zip(nls._X,nls._Y))
#     r = map(X_i->eval_y(nls._fit_model,X_i,θ),nls._X)
#     nls._Y-r
# end


# function NLS_Solver.eval_r_J(nls::NLS_Zygote_From_Fit_Model, θ::AbstractVector{T}) where {T}
    
#     r = eval_r(nls,θ)

#     J = jacobian(θ->eval_r(nls,θ), θ)

#     r,J
# end
