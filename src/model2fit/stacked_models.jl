# Stack models together

export Model2Fit_Stacked


@doc raw"""

Define a composite model, where each submodel is responsible of a given ROI

| (1)        | (2)       | ... 
|------------|-----------| ...
| models[1]  | model[2]  | ...

Domain (1) range is defined as 1:ROI_sizes[1]
Domain (2) range is defined as (ROI_sizes[1]+1):ROI_sizes[2]
etc...
"""
struct Model2Fit_Stacked{MODEL2FIT<:Abstract_Model2Fit} <: Abstract_Model2Fit
    _models::Vector{MODEL2FIT}
    _ROI_lengths::Vector{Int}
end 

parameter_size(m::Model2Fit_Stacked) = sum(NLS_Fit.parameter_size,m._models)

function accumulate_y!(m::Model2Fit_Stacked,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(X) == length(Y)
    @assert length(X) == sum(m._ROI_lengths)
    @assert length(m._models) == length(m._ROI_lengths)

    if isempty(m._models) return Y end
    
    XY_range_first = 1
    θ_range_first = 1

    for (model,roi_length) in zip(m._models,m._ROI_lengths)

        θ_length = parameter_size(model)

        XY_range = XY_range_first:(XY_range_first + roi_length -1)
        θ_range = θ_range_first:(θ_range_first + θ_length - 1)

        accumulate_y!(model,(@view Y[XY_range]),(@view X[XY_range]),(@view θ[θ_range]))

        XY_range_first += roi_length 
        θ_range_first += θ_length 
    end

    Y
end
