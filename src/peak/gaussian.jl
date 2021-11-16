export Gaussian_Peak
export parameter_size, eval_y
struct Gaussian_Peak <: Abstract_Model2Fit_Peak
end

parameter_size(::Gaussian_Peak) = 3

function eval_y(m::Gaussian_Peak,x::Real,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    
    h=θ[1]
    μ=θ[2]
    σ=θ[3]

    t = ((x-μ)/σ)^2
    
    h*exp(-t/2)
end
