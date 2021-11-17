# Implement affine recalibration
#
export Recalibration_Affine, eval_x

@doc raw"""
Create a calibrable model

The calibrated ``X`` viewed by the wrapped model when calling

```julia
veal_y(calibrable_model,X_hat,[θ, θ_A, θ_B])
```

is 

```math
X(\hat{X}) = L_A(\hat{X}) X_A θ_A +  L_B(\hat{X}) X_B θ_B
```

In other terms, the wrapped model is called as follows:


```julia
eval_y(model,X(X_hat, θ_A, θ_B), [θ])
```

With this parametrization when θ_A=θ_B=1 
"""
struct Recalibration_Affine{T<:Real, Model2Calibrate <: Abstract_Model2Fit}  <: Abstract_Model2Fit
    _m2c::Model2Calibrate
    _hat_to_cal_A::Pair{T,T}
    _hat_to_cal_B::Pair{T,T}
end

# function Recalibration_Affine(m2c::Abstract_Model2Fit,
#                               x_hat_to_cal_A::Pair{T_HAT,T},
#                               x_hat_to_cal_B::Pair{T_HAT,T}) where {T_HAT,T}
#     @assert( first(x_hat_to_cal_A) != first(x_hat_to_cal_B) ) # singular transform

#     Recalibration_Affine(m2c,x_hat_to_cal_A,x_hat_to_cal_B)
# end

function Recalibration_Affine(m2c::Abstract_Model2Fit,x_hat_A::T,x_hat_B::T) where {T<:Real}
    Recalibration_Affine(m2c,x_hat_A=>x_hat_A,x_hat_B=>x_hat_B)
end

parameter_size(m::Recalibration_Affine) = parameter_size(m._m2c)+2

function _transform_affine(X_hat::AbstractVector{T_HAT},
                           hat_to_cal_A::Pair{T,T},
                           hat_to_cal_B::Pair{T,T},
                           θ_A::T_θ,
                           θ_B::T_θ) where {T_HAT,T,T_θ}
    X_T = promote_type(T_HAT,T,T_θ)
    
    Δ = 1/(first(hat_to_cal_B)-first(hat_to_cal_A))
    basis_A(x) = (first(hat_to_cal_B)-x)*Δ*last(hat_to_cal_A)
    basis_B(x) = (x-first(hat_to_cal_A))*Δ*last(hat_to_cal_B)
    transform(x) = basis_A(x)*θ_A + basis_B(x)*θ_B

    X=similar(X_hat,X_T)
    @. X = transform(X_hat)

    X
end

function eval_x(m::Recalibration_Affine,X_hat::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)

    s=parameter_size(m._m2c)
    θ_A=θ[s+1]
    θ_B=θ[s+2]
    
    X = _transform_affine(X_hat,
                          m._hat_to_cal_A,m._hat_to_cal_B,
                          θ_A, θ_B)
end

function eval_y!(m::Recalibration_Affine,Y::AbstractVector,X_hat::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)

    s=parameter_size(m._m2c)
    θ_A=θ[s+1]
    θ_B=θ[s+2]
 
    X = _transform_affine(X_hat,
                          m._hat_to_cal_A,m._hat_to_cal_B,
                          θ_A, θ_B)
                         
    eval_y!(m._m2c,Y,X,@view θ[1:parameter_size(m._m2c)])
end
