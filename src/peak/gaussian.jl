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
