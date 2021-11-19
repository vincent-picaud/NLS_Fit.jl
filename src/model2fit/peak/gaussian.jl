export Gaussian_Peak

struct Gaussian_Peak <: Abstract_Model2Fit_Peak
end

parameter_size(::Gaussian_Peak) = 3

function eval_y!(m::Gaussian_Peak,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(Y) == length(X)
    
    h=θ[1]
    μ=θ[2]
    σ=θ[3]

    @. Y += h*exp(-((X-μ)/σ)^2/2)

    Y
end

# ================================================================

using StaticArrays: @SVector, SVector

struct Const_μ_Gaussian_Peak{μ_T <: Real} <: Abstract_Model2Fit_Peak
    _μ::μ_T
end

parameter_size(::Const_μ_Gaussian_Peak) = 2

function eval_y!(m::Const_μ_Gaussian_Peak,Y::AbstractVector,X::AbstractVector,θ::AbstractVector{T}) where {T}
    @assert length(θ) == parameter_size(m)
    @assert length(Y) == length(X)

    h=θ[1]
    μ=m._μ
    σ=θ[2]

    θ̂ = @SVector T[h,μ,σ]

    eval_y!(Gaussian_Peak(),Y,X,θ̂)
end