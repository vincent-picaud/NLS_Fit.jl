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
    _X_ranges::Vector{UnitRange{Int}}
    _θ_ranges::Vector{UnitRange{Int}}

    _ROI_lengths::Vector{Int}

    function Model2Fit_Stacked(models::ABSVECT_MODEL2FIT,
                               ROI_lengths::AbstractVector{Int}) where {MODEL2FIT<:Abstract_Model2Fit,
                                                                        ABSVECT_MODEL2FIT<:AbstractVector{MODEL2FIT}}

        @assert length(models) == length(ROI_lengths)
        
        # Compute θ and X ranges associated to each submodel
        #
        n_models = length(models)
        
        X_ranges = Vector{UnitRange{Int}}(undef,n_models)
        θ_ranges = Vector{UnitRange{Int}}(undef,n_models)
        
        if n_models>0

            X_range_first = 1
            θ_range_first = 1

            resize!(X_ranges,0)
            resize!(θ_ranges,0)
            
            for (model,roi_length) in zip(models,ROI_lengths)

                X_range = X_range_first:(X_range_first + roi_length -1)
                push!(X_ranges,X_range)

                θ_length = parameter_size(model)
                θ_range = θ_range_first:(θ_range_first + θ_length - 1)
                push!(θ_ranges,θ_range)
                
                X_range_first += roi_length 
                θ_range_first += θ_length 
            end
            
        end 


        new{MODEL2FIT}(models,X_ranges,θ_ranges)

        
    end
end 

# Helper  ================
#
_get_XY(model::Model2Fit_Stacked,submodel_idx::Int,XY::AbstractVector) = @view XY[model._X_ranges[submodel_idx]]
_get_θ(model::Model2Fit_Stacked,submodel_idx::Int,θ::AbstractVector) = @view θ[model._θ_ranges[submodel_idx]]

# Visit  ================
#
visit_submodel_size(model::Model2Fit_Stacked) = length(model._models)

visit_get_submodel(model::Model2Fit_Stacked,submodel_idx::Int) = model._models[submodel_idx]
visit_get_Y(model::Model2Fit_Stacked,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = _get_XY(model,submodel_idx,Y)
visit_get_X(model::Model2Fit_Stacked,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = _get_XY(model,submodel_idx,X)
visit_get_θ(model::Model2Fit_Stacked,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = _get_θ(model,submodel_idx,θ)

# Interface ================
#
parameter_size(m::Model2Fit_Stacked) = last(last(m._θ_ranges))

function accumulate_y!(model::Model2Fit_Stacked,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    n_submodel = visit_submodel_size(model)

    for i in 1:n_submodel
        submodel = visit_get_submodel(model,i)
        subY = _get_XY(model,i,Y)
        subX = _get_XY(model,i,X)
        subθ = _get_θ(model,i,θ)

        accumulate_y!(submodel,subY,subX,subθ)
    end 

    Y
end
