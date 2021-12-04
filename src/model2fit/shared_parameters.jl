# Create a model with shared parameters ****************
#
export Model2Fit_Shared_Parameters
export get_model, get_model_θ

struct Model2Fit_Shared_Parameters{MODEL <: Abstract_Model2Fit,
                                   INDICES <: AbstractVector{Int}}  <: Abstract_Model2Fit
    _model::MODEL
    _indices::INDICES

    function Model2Fit_Shared_Parameters(model::MODEL,
                                         indices::INDICES) where {MODEL <: Abstract_Model2Fit,
                                                                  INDICES <: AbstractVector{Int}}
        @assert (length(indices) == 0) || ( (1 ≤ minimum(indices)) && (maximum(indices) ≤ parameter_size(model)) )
        @assert allunique(indices)
        
        new{MODEL,INDICES}(model, indices)

    end
    
end

# Convenience methods ================
#
@doc raw"""
```julia
get_model(mp::Model2Fit_Shared_Parameters) -> Absatrct_Model2Fit
```

Get back the wrapped model
"""
get_model(mp::Model2Fit_Shared_Parameters) = mp._model

@doc raw"""
```julia
get_model_θ(mp::Model2Fit_Shared_Parameters,θ::AbstractVector) -> θ::AbstractVector
```

Retrieve the parameter vector θ associated to the wrapped model [`get_model`](@ref).
"""
function get_model_θ(model::Model2Fit_Shared_Parameters,θ::AbstractVector)
    # Count unshared parameters 
    #
    θ_not_shared_size = parameter_size(model._model) - length(model._indices)
    θ_not_shared = @view θ[1:θ_not_shared_size]

    # Retrieve shared value ================
    #
    # This is the last parameter
    #
    @assert θ_not_shared_size+1 == length(θ)
    θ_shared = θ[end]

    # Create the new θ with θ_shared set a [indices]
    #
    insert_some_elements(θ_not_shared, model._indices, θ_shared)
end

# Visit  ================
#
# note: no need to check
#   @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)
#   @assert length(X) == length(Y)
# as this is done higher in the visit() method
#
visit_submodel_size(model::Model2Fit_Shared_Parameters) = 1
visit_get_submodel(model::Model2Fit_Shared_Parameters,submodel_idx::Int) = model._model
visit_get_X(model::Model2Fit_Shared_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = X
visit_get_Y(model::Model2Fit_Shared_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = Y
visit_get_θ(model::Model2Fit_Shared_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = get_model_θ(model,θ)


# Interface ================
#
parameter_size(mp::Model2Fit_Shared_Parameters) = parameter_size(mp._model) - length(mp._indices) + 1 

function accumulate_y!(mp::Model2Fit_Shared_Parameters,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(mp)
    
    θ_model = get_model_θ(mp,θ)
    model = get_model(mp)

    accumulate_y!(model,Y,X,θ_model)
end

