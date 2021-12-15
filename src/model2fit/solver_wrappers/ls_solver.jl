# Defines some function to find linear parameters
#
export solve_linear_parameters

@doc raw"""

Find optimal linear parameters
"""
function solve_linear_parameters(model::Abstract_Model2Fit,
                                 X::AbstractVector,
                                 Y::AbstractVector,
                                 θ::AbstractVector,
                                 linear_param_idx::AbstractVector{Int})
    @assert length(X) == length(Y)
    @assert length(θ) == parameter_size(model)

    # Build J of the least squares problem
    #
    # | J.θ_linear - Y |^2
    #
    θ_copy = copy(θ)
    θ_copy[linear_param_idx] .= 0

    J = zeros(promote_type(eltype(X),eltype(Y),eltype(θ)),
              length(X),
              length(linear_param_idx))
    
    for (col_idx,θ_idx) in enumerate(linear_param_idx)
        θ_copy[θ_idx] = 1
        accumulate_y!(model,(@view J[:,col_idx]),X,θ_copy)
        θ_copy[θ_idx] = 0
    end
    
    # Solve the problem
    #
    # H.θ_linear = rhs
    #
    # where H = J'J, rhs = J'Y
    #
    rhs = J'*Y 
    H = J'*J

    θ_linear = H \ rhs

    θ_copy[linear_param_idx] .= θ_linear

    θ_copy
end
