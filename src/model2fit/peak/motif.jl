export Peak_Motif

using StaticArrays: @SVector, SVector

# note:
# ----------------------------------------------------------------
# as before uses a parametrized type to limit memory allocs
#
struct Peak_Motif{P<:Abstract_Model2Fit_Peak} <: Abstract_Model2Fit_Peak
    _peak::P
    _profile::Matrix{Float64} # position, height
end

# ================================================================

# Gaussian peak specialisation
# θ = (h,σ)
#
parameter_size(pm::Peak_Motif{Gaussian_Peak}) = 2

function accumulate_y!(pm::Peak_Motif{Gaussian_Peak},Y::AbstractVector,X::AbstractVector,θ::AbstractVector{T}) where {T}
    @assert length(θ) == parameter_size(pm)

    h_glob = θ[1]
    σ_glob = θ[2]

    n = size(pm._profile,1)

    for i in 1:n
        μ_loc = pm._profile[i,1]
        h_loc = pm._profile[i,2]

        accumulate_y!(pm._peak, Y, X, @SVector T[ h_glob*h_loc,μ_loc,σ_glob])
    end

    Y
end 
